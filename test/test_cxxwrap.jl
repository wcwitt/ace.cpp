module LibACE
    using CxxWrap
    @wrapmodule(joinpath("../build/","libace-cxxwrap"))
    function __init__()
        @initcxx
    end
end

@show LibACE.naive_sph_harm(1,0,0.5,0.5)
@show LibACE.naive_sph_harm_xyz(1,0,0.5,0.5,0.5)
@show LibACE.spherical_bessel_radial(1,0,0.5,1.0)
@show LibACE.determine_basis(5.0,1.0)
