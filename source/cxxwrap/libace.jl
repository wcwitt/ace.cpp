module LibACE
    using CxxWrap
    @wrapmodule(joinpath(@__DIR__,"../../build/libace-cxxwrap"))
    function __init__()
        @initcxx
    end
end
