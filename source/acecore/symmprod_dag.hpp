#pragma once

#include <map>
#include <vector>

struct SparseSymmProdDAG
{
    std::vector<std::pair<int,int>> nodes;
    int num1;
    int numstore;
    std::vector<int> projection;

    // note factory function below
};

double _score_partition(std::vector<int> p);

std::vector<int> _get_ns(
    std::vector<std::vector<int>> p,
    std::vector<std::vector<int>> specnew,
    std::map<std::vector<int>,int> specnew_dict);

std::vector<int> _find_partition(
    std::vector<int> kk,
    std::vector<std::vector<int>> specnew,
    std::map<std::vector<int>,int> specnew_dict);

int _insert_partition(
    std::vector<std::pair<int,int>> &nodes,
    std::vector<std::vector<int>> &specnew,
    std::map<std::vector<int>,int> &specnew_dict,
    std::vector<int> kk,
    std::vector<int> p,
    int ikk,
    std::vector<std::vector<int>> specN);

SparseSymmProdDAG BuildSparseSymmProdDAG(std::vector<std::vector<int>> spec);
