#include <algorithm>

#include "../partitions.hpp"

#include "symmprod_dag.hpp"

//------------------------------------------------
//                           _score_partition
//------------------------------------------------
/*
_score_partition(p) = isempty(p) ? Inf : (1e9 * length(p) + maximum(p))`
*/
double _score_partition(std::vector<int> p)
{
    if (p.size() == 0)
        return std::numeric_limits<double>::infinity();

    // todo: replace this with proper std::max-type function
    int max_p = p[0];
    for (auto p : p) {
        max_p = std::max(p, max_p);
    }
    return 1e9 * p.size() + max_p;
}

//------------------------------------------------
//                                    _get_ns
//------------------------------------------------
//  Roughly, takes a partition (set of subsets) and checks specnew_dict to see
//  whether each subset has been visited before
//  (Note that each single particle basis function is guaranteed to exist in specnew_dict)
//  If yes, stores the output of specnew_dict for each subset.
//  If no (for any subset?), returns empty vector (?)
//
/*
function _get_ns(p, specnew, specnew_dict)
   out = Vector{Int}(undef, length(p))
   for (i, kk_) in enumerate(p)
      if haskey(specnew_dict, kk_)
         out[i] = specnew_dict[kk_]
      else
         return Int[]
      end
   end
   return out
end
*/
std::vector<int> _get_ns(
    std::vector<std::vector<int>> p,
    std::vector<std::vector<int>> specnew,
    std::map<std::vector<int>,int> specnew_dict)
{
    auto out = std::vector<int>(p.size());
    for (int i=0; i<p.size(); ++i) {
        auto kk_ = p[i];
        if (specnew_dict.find(kk_) == specnew_dict.end()) {
            return std::vector<int>();
        } else {
            out[i] = specnew_dict[kk_];
        }
    }
    return out;
}

//------------------------------------------------
//                            _find_partition
//------------------------------------------------
/*
function _find_partition(kk, specnew, specnew_dict)
   worstp = _get_ns([ [k] for k in kk ], specnew, specnew_dict)
   @assert worstp == kk
   bestp = worstp
   bestscore = _score_partition(bestp)

   for ip in partitions(1:length(kk))
      p = _get_ns([ kk[i] for i in ip ], specnew, specnew_dict)
      score = _score_partition(p)
      if !isempty(p) && score < bestscore
         bestp = p
         bestscore = score
      end
   end

   return bestp
end
*/
std::vector<int> _find_partition(
    std::vector<int> kk,
    std::vector<std::vector<int>> specnew,
    std::map<std::vector<int>,int> specnew_dict)
{
    std::vector<std::vector<int>> partitioned_kk;
    for (auto k : kk)
        partitioned_kk.push_back({k});
    auto worstp = _get_ns(partitioned_kk, specnew, specnew_dict);
    //@assert worstp == kk
    auto bestp = worstp;
    auto bestscore = _score_partition(bestp);

    std::vector<int> range;
    for (int i=0; i<kk.size(); ++i)
        range.push_back(i);
    for (auto partition : partitions(range)) {
        std::vector<std::vector<int>> partitioned_kk;
        for (auto subset : partition) {
            partitioned_kk.push_back(std::vector<int>(subset.size()));
            for (int i=0; i<subset.size(); ++i)
                partitioned_kk[partitioned_kk.size()-1][i] = kk[subset[i]];
        }
        auto p = _get_ns(partitioned_kk, specnew, specnew_dict);
        auto score = _score_partition(p);
        if (p.size()>0 && score<bestscore) {
            bestp = p;
            bestscore = score;
        }
    }

    return bestp;
}

//------------------------------------------------
//                          _insert_partition
//------------------------------------------------
//  TODO: ikk is not used in either function
/*
# return value is the number of fake nodes added to the dag
function _insert_partition!(nodes, specnew, specnew_dict,
                            kk, p,
                            ikk, specN)
   if length(p) == 2
      newnode = BinDagNode((p[1], p[2]))
      push!(nodes, newnode)
      push!(specnew, kk)
      specnew_dict[kk] = length(specnew)
      return 0
   else
      # @show kk, p
      # @infiltrate
      # reduce the partition by pushing a new node
      push!(nodes, BinDagNode((p[1], p[2])))
      kk1 = sort(vcat(specnew[p[1]], specnew[p[2]]))
      push!(specnew, kk1)
      specnew_dict[kk1] = length(specnew)
      # and now recurse with the reduced partition
      return 1 + _insert_partition!(nodes, specnew, specnew_dict,
                         kk, vcat( [length(nodes)], p[3:end] ),
                         ikk, specN)
   end
end
*/
int _insert_partition(
    std::vector<std::pair<int,int>> &nodes,
    std::vector<std::vector<int>> &specnew,
    std::map<std::vector<int>,int> &specnew_dict,
    std::vector<int> kk,
    std::vector<int> p,
    int ikk,
    std::vector<std::vector<int>> specN)
{
    if (p.size() == 2) {
        auto newnode = std::pair<int,int>{p[0], p[1]};
        nodes.push_back(newnode);
        specnew.push_back(kk);
        specnew_dict[kk] = specnew.size()-1;  // TODO: note -1 for C++ indexing
        return 0;
    } else {
        // reduce the partition by pushing a new node
        nodes.push_back({p[0], p[1]});
        auto kk1 = specnew[p[0]];
        kk1.insert(kk1.end(), specnew[p[1]].begin(), specnew[p[1]].end());
        std::sort(kk1.begin(), kk1.end());
        specnew.push_back(kk1);
        specnew_dict[kk1] = specnew.size()-1;
        // and now recurse with the reduced partition
        p.erase(p.begin());
        p.erase(p.begin());
        p.insert(p.begin(), nodes.size()-1);  // TODO: note this uses 0-indexing
        return 1 + _insert_partition(
            nodes, specnew, specnew_dict, kk, p, ikk, specN);
    }
}

//------------------------------------------------
//                     BuildSparseSymmProdDAG
//------------------------------------------------
/*
"""
Construct the DAG used to evaluate an AA basis and returns it as a `SparseSymmProdDAG`

Arguments
* `spec` : AA basis specification, list of vectors of integers / indices pointing into A 

Kwargs: 
* `filter = _-> true` : 
* `verbose = false` : print some information about the 
"""
function SparseSymmProdDAG(spec::AbstractVector; 
                           filter = _->true, 
                           verbose = false, 
                           T = Float64)
   @assert issorted(length.(spec))
   @assert all(issorted, spec)
   # we need to separate them into 1-p and many-p
   spec1 = spec[ length.(spec) .== 1 ]
   IN = (length(spec1)+1):length(spec)
   specN = spec[IN]

   # start assembling the dag
   nodes = BinDagNode[]
   sizehint!(nodes, length(spec))
   specnew = Vector{Int}[]
   specnew_dict = Dict{Vector{Int}, Int}()
   sizehint!(specnew, length(spec))

   # add the full 1-particle basis (N=1) into the dag
   num1 = maximum( maximum(vv) for vv in spec )
   for i = 1:num1
      push!(nodes, BinDagNode((i, 0)))
      push!(specnew, [i])
      specnew_dict[ [i] ] = length(specnew)
   end

   # now we can construct the rest
   extranodes = 0
   for (ikk, kk) in enumerate(specN)
      # find a good partition of kk
      p = _find_partition(kk, specnew, specnew_dict)
      extranodes += _insert_partition!(nodes, specnew, specnew_dict,
                                       kk, p, ikk, specN)
   end

   verbose && @info("Extra nodes inserted into the dag: $extranodes")
   numstore = length(nodes)
   num1old = num1

   projection = [ specnew_dict[vv] for vv in spec ]

   # re-organise the dag layout to minimise numstore
   # nodesfinal, num1, numstore = _reorder_dag!(nodes)

   return SparseSymmProdDAG{T}(nodes, num1, numstore, projection)
end
*/
SparseSymmProdDAG BuildSparseSymmProdDAG(std::vector<std::vector<int>> spec)
{
    //@assert issorted(length.(spec))
    //@assert all(issorted, spec)
    // we need to separate them into 1-p and many-p
    std::vector<std::vector<int>> spec1;
    std::vector<std::vector<int>> specN;
    for (auto s : spec) {
        if (s.size() == 1) {
            spec1.push_back(s);
        } else {
            specN.push_back(s);
        }
    }
    
    // start assembling the dag
    std::vector<std::pair<int,int>> nodes;
    nodes.reserve(spec.size());
    std::vector<std::vector<int>> specnew;
    specnew.reserve(spec.size());
    // strictly, unordered_map is a better analog than map to a Julia Dict()
    // but that doesn't seem to easily allow vectors as keys
    // see: https://jimmy-shen.medium.com/stl-map-unordered-map-with-a-vector-for-the-key-f30e5f670bae
    std::map<std::vector<int>,int> specnew_dict;

    // add the full 1-particle basis (N=1) into the dag
    int num1 = 0;
    for (auto vv : spec) {
        for (auto v : vv) {
            num1 = std::max(v+1, num1); // TODO: the +1 is because C++ uses zero indexing
        }
    }
    for (int i=0; i<num1; ++i) {
        nodes.push_back({i, 0});
        specnew.push_back({i});
        specnew_dict[{i}] = specnew.size()-1;
    }

    // now we can construct the rest
    int extranodes = 0;
    for (int ikk=0; ikk<specN.size(); ++ikk) {
        auto kk = specN[ikk];
        // find a good partition of kk
        auto p = _find_partition(kk, specnew, specnew_dict);
        extranodes += _insert_partition(
            nodes, specnew, specnew_dict, kk, p, ikk, specN);
    }

    int numstore = nodes.size();
    
    std::vector<int> projection(spec.size());
    for (int i=0; i<projection.size(); ++i)
        projection[i] = specnew_dict[spec[i]];

    return {nodes, num1, numstore, projection};
}
