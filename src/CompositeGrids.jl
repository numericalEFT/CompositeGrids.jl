module CompositeGrids
using StaticArrays

include("old/grid.jl")
export Grid

include("grid/simple.jl")
export SimpleGrid

include("grid/composite.jl")
export CompositeGrid

include("grid/interpolate.jl")
export Interp

end # module
