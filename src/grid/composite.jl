module CompositeGrid

using StaticArrays, FastGaussQuadrature
include("simple.jl")
using .SimpleGrid

struct Composite{T<:AbstractFloat} <: AbstractGrid{T}
end


end
