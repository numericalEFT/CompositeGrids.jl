module SimpleGrid

using StaticArrays, FastGaussQuadrature

include("chebyshev.jl")

abstract type AbstractGrid{T} end

struct Arbitrary{T<:AbstractFloat} <: AbstractGrid{T}
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}

    function Arbitrary{T}(grid) where {T<:AbstractFloat}
        bound = [grid[1],grid[end]]
        size = length(grid)
        return new{T}(bound, size, grid)
    end
end

function Base.floor(grid::Arbitrary{T}, x) where {T}
    if x <= grid.grid[1]
        return 1
    elseif x >= grid.grid[end]
        return grid.size-1
    end

    result = searchsortedfirst(grid.grid, x)-1
    return result
end

Base.getindex(grid::Arbitrary, i) = grid.grid[i]
Base.firstindex(grid::Arbitrary) = 1
Base.lastindex(grid::Arbitrary) = grid.size

struct Uniform{T<:AbstractFloat} <: AbstractGrid{T}
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}

    function Uniform{T}(bound, N) where {T<:AbstractFloat}
        Ntot = N - 1
        interval = (bound[2]-bound[1])/Ntot
        grid = bound[1] .+ Vector(1:N) .* interval .- ( interval )

        return new{T}(bound, N, grid)
    end
end

function Base.floor(grid::Uniform{T}, x) where {T}
    result = (x-grid.grid[1])/(grid.grid[end]-grid.grid[1])*(grid.size-1)+1
    if result <=0
        return 1
    elseif result >= grid.size
        return grid.size-1
    else
        return Base.floor(Int, result)
    end

end

Base.getindex(grid::Uniform, i) = grid.grid[i]
Base.firstindex(grid::Uniform) = 1
Base.lastindex(grid::Uniform) = grid.size

struct BaryCheb{T<:AbstractFloat} <: AbstractGrid{T}
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}
    weight::Vector{T}

    function BaryCheb{T}(bound, N) where {T<:AbstractFloat}
        order = N
        x, w =barychebinit(order)
        grid = zeros(T, N)
        a, b = bound[1], bound[2]
        weight = (b - a) / 2  .* w
        grid = (a + b) / 2 .+ (b - a) / 2 .* x

        return new{T}(bound, N, grid, weight)
    end
end

function Base.floor(grid::BaryCheb{T}, x) where {T}
    if x <= grid.grid[1]
        return 1
    elseif x >= grid.grid[end]
        return grid.size-1
    end

    result = searchsortedfirst(grid.grid, x)-1
    return result
end

Base.getindex(grid::BaryCheb, i) = grid.grid[i]
Base.firstindex(grid::BaryCheb) = 1
Base.lastindex(grid::BaryCheb) = grid.size

struct GaussLegendre{T<:AbstractFloat} <: AbstractGrid{T}
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}
    weight::Vector{T}

    function GaussLegendre{T}(bound, N) where {T<:AbstractFloat}
        order = N
        x, w = gausslegendre(order)
        grid = zeros(T, N)
        a, b = bound[1], bound[2]
        weight = (b - a) / 2  .* w
        grid = (a + b) / 2 .+ (b - a) / 2 .* x

        return new{T}(bound, N, grid, weight)
    end
end

function Base.floor(grid::GaussLegendre{T}, x) where {T}
    if x <= grid.grid[1]
        return 1
    elseif x >= grid.grid[end]
        return grid.size-1
    end

    result = searchsortedfirst(grid.grid, x)-1
    return result
end

Base.getindex(grid::GaussLegendre, i) = grid.grid[i]
Base.firstindex(grid::GaussLegendre) = 1
Base.lastindex(grid::GaussLegendre) = grid.size


end
