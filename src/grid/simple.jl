"""
Basic grids including common grids like arbitrary grids, uniform grids, log grids,
and optimized grids like barycheb for interpolation and gausslegendre for integration.

"""
module SimpleG

export AbstractGrid, Uniform, BaryCheb, GaussLegendre, Arbitrary, Log, denseindex
export GridBoundType, OpenBound, ClosedBound, PeriodicBound
using StaticArrays, FastGaussQuadrature

using ..BaryChebTools
const barychebinit = BaryChebTools.barychebinit
const vandermonde = BaryChebTools.vandermonde
const invvandermonde = BaryChebTools.invvandermonde
# include("chebyshev.jl")
# export barychebinit, barycheb

"""
All Grids are derived from AbstractGrid; ClosedGrid has bound[1], bound[2] == grid[1], grid[end],
while OpenGrid has bound[1]<grid[1]<grid[end]<bound[2]
"""
abstract type AbstractGrid{T} <: AbstractArray{T,1} end

abstract type GridBoundType end
struct OpenBound <: GridBoundType end
struct ClosedBound <: GridBoundType end
struct PeriodicBound <: GridBoundType end
const OPENBOUND = OpenBound()
const CLOSEDBOUND = ClosedBound()
const PERIODICBOUND = PeriodicBound()

# decide by relation between bound and [grid[1], grid[end]] for open and close
function BoundType(::Type, bound, gpbound)
    if bound[1] == gpbound[1] && bound[2] == gpbound[2]
        return ClosedBound()
    else
        return OpenBound()
    end
end

# periodic if claimed periodic
BoundType(::Type{<:PeriodicBound}, bound, gpbound) = PeriodicBound()

# for grids

# abstract type OpenGrid{T} <: AbstractGrid{T} end
# abstract type ClosedGrid{T} <: AbstractGrid{T} end
# abstract type PeriodicGrid{T} <: AbstractGrid{T} end

"""

    struct Arbitrary{T<:AbstractFloat} <: AbstractGrid

Arbitrary grid generated from given sorted grid.

#Members:
- `bound` : boundary of the grid
- `size` : number of grid points
- `grid` : grid points
- `weight` : integration weight

#Constructor:
-    function Arbitrary{T}(grid) where {T<:AbstractFloat}
"""
struct Arbitrary{T<:Real,BT} <: AbstractGrid{T}
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}
    weight::Vector{Float64}

    """
        function Arbitrary{T}(grid; bound=[grid[1], grid[end]]) where {T<:AbstractFloat}

    create Arbitrary from grid.
    """
    function Arbitrary{T,BTIN}(grid;
        bound=[grid[1], grid[end]]) where {T<:Real,BTIN}
        # allow customized bound that's different from [grid[1], grid[end]]
        @assert bound[1] <= grid[1]
        @assert bound[2] >= grid[end]
        @assert allunique(grid) "Grid points should be unique!"
        if !(issorted(grid))
            @warn "Input grid is not sorted, using sorted instead."
            grid = sort(grid)
        end
        BT = typeof(BoundType(BTIN, bound, [grid[1], grid[end]]))
        if BT == PeriodicBound && (bound[1] == grid[1] && bound[2] == grid[end])
            size = length(grid) - 1
            weight = zeros(Float64, size)
            for i in 1:size
                if i == 1
                    if size != 1
                        weight[1] = 0.5 * (grid[2] - grid[1]) + 0.5 * (bound[2] - grid[end-1])
                    else
                        weight[1] = 0 # allow arbitrary grid for 1 gridpoint, but integrate undefined
                    end
                elseif i == size
                    weight[end] = 0.5 * (grid[end] - grid[end-1])
                else
                    weight[i] = 0.5 * (grid[i+1] - grid[i-1])
                end
            end
            return new{T,BT}(bound, size, grid[1:end-1], weight)
        else
            size = length(grid)
            weight = zeros(Float64, size)
            for i in 1:size
                if i == 1
                    if size != 1
                        weight[1] = 0.5 * (grid[2] - grid[1])
                    else
                        weight[1] = 0 # allow arbitrary grid for 1 gridpoint, but integrate undefined
                    end
                elseif i == size
                    weight[end] = 0.5 * (grid[end] - grid[end-1])
                else
                    weight[i] = 0.5 * (grid[i+1] - grid[i-1])
                end
            end
            return new{T,BT}(bound, size, grid, weight)
        end
    end
end
Arbitrary{T}(grid;
    bound=[grid[1], grid[end]],
    boundtype::BTIN=CLOSEDBOUND) where {T<:Real,BTIN} = Arbitrary{T,BTIN}(grid; bound=bound)
function Arbitrary(grid;
    bound=[grid[1], grid[end]],
    isperiodic=false)
    T = eltype(grid)
    BTIN = isperiodic ? PeriodicBound : ClosedBound
    Arbitrary{T,BTIN}(grid; bound=bound)
end


"""
    function Base.floor(grid::AbstractGrid, x) #where {T}

use basic searchsorted function to find the index of largest
grid point smaller than x.

return 1 for x<grid[1] and grid.size-1 for x>grid[end].
"""
function Base.floor(grid::AbstractGrid, x) #where {T}
    if x <= grid.grid[1]
        return 1
    elseif x >= grid.grid[end]
        if grid.size != 1
            return grid.size - 1
        else
            return 1
        end
    end

    result = searchsortedfirst(grid.grid, x) - 1
    return Base.floor(Int, result)
end

function Base.floor(grid::Arbitrary{T,PeriodicBound}, x) where {T}
    fracx = (x - grid.grid[1]) / (grid.bound[2] - grid.bound[1])
    fracx = (fracx % 1 + 1) % 1
    x = grid.grid[1] + fracx * (grid.bound[2] - grid.bound[1])
    if x < grid.grid[1] || x > grid.grid[end]
        return grid.size
    else
        result = searchsortedfirst(grid.grid, x) - 1
        return Base.floor(Int, result)
    end
end

Base.length(grid::AbstractGrid) = grid.size
Base.size(grid::AbstractGrid) = (grid.size,)
Base.size(grid::AbstractGrid, I::Int) = grid.size

Base.view(grid::AbstractGrid, inds...) = Base.view(grid.grid, inds...)
# set is not allowed for grids
Base.getindex(grid::AbstractGrid, i) = grid.grid[i]
Base.firstindex(grid::AbstractGrid) = 1
Base.lastindex(grid::AbstractGrid) = grid.size
# iterator
Base.iterate(grid::AbstractGrid) = (grid.grid[1], 1)
Base.iterate(grid::AbstractGrid, state) = (state >= grid.size) ? nothing : (grid.grid[state+1], state + 1)

# Base.IteratorSize
Base.IteratorSize(::Type{AbstractGrid{T}}) where {T} = Base.HasLength()
Base.IteratorEltype(::Type{AbstractGrid{T}}) where {T} = Base.HasEltype()
Base.eltype(::Type{AbstractGrid{T}}) where {T} = eltype(T)

"""
    show(io::IO, grid::AbstractGrid)

Write a text representation of the AbstractGrid 
`grid` to the output stream `io`.
"""
function Base.show(io::IO, grid::AbstractGrid; isSimplified=false)
    if isSimplified
        print(io, "$(typeof(grid)): 1D Grid with $(grid.size) grid points.\n")
    else
        print(io,
            "$(typeof(grid)): 1D Grid with:\n"
            * "- bound: $(grid.bound)\n"
            * "- size: $(grid.size)\n"
            * "- grid: $(grid.grid)\n"
        )
    end
end

"""
    struct Uniform{T<:AbstractFloat} <: AbstractGrid

Uniform grid generated on [bound[1], bound[2]] with N points

#Members:
- `bound` : boundary of the grid
- `size` : number of grid points
- `grid` : grid points
- `weight` : integration weight

#Constructor:
-    function Uniform{T}(bound, size) where {T<:AbstractFloat}
"""
struct Uniform{T<:AbstractFloat,BT} <: AbstractGrid{T}
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}
    weight::Vector{T}
end
"""
    function Uniform{T}(bound, N) where {T<:AbstractFloat}

create Uniform grid.
"""
function Uniform{T,BTIN}(bound, N;
    gpbound=bound) where {T<:AbstractFloat,BTIN}
    @assert bound[1] < bound[2]
    @assert bound[1] <= gpbound[1]
    @assert bound[2] >= gpbound[2]
    @assert N > 0
    BT = typeof(BoundType(BTIN, bound, gpbound))
    if BT == PeriodicBound && (bound[1] == gpbound[1] && bound[2] == gpbound[2])
        Ntot = N
        interval = (gpbound[2] - gpbound[1]) / Ntot
        grid = gpbound[1] .+ Vector(1:N) .* interval .- (interval)
        weight = interval .* ones(N)
        return Uniform{T,BT}(bound, N, grid, weight)
    else
        Ntot = N - 1
        interval = (gpbound[2] - gpbound[1]) / Ntot
        grid = gpbound[1] .+ Vector(1:N) .* interval .- (interval)
        weight = similar(grid)
        for i in 1:N
            if i == 1
                weight[1] = 0.5 * (grid[2] - bound[1])
            elseif i == N
                weight[end] = 0.5 * (bound[2] - grid[end-1])
            else
                weight[i] = 0.5 * (grid[i+1] - grid[i-1])
            end
        end
        return Uniform{T,BT}(bound, N, grid, weight)
    end
end
function Uniform{T}(bound, N;
    gpbound=bound,
    boundtype::BTIN=CLOSEDBOUND) where {T,BTIN}
    return Uniform{T,BTIN}(bound, N; gpbound=gpbound)
end
function Uniform(bound, N;
    gpbound=bound,
    isperiodic=false)
    T = typeof(bound[1])
    BTIN = isperiodic ? PeriodicBound : ClosedBound
    return Uniform{T,BTIN}(bound, N; gpbound=gpbound)
end

"""
    function Base.floor(grid::Uniform{T}, x) where {T}

find the index of largest
grid point smaller than x.

return 1 for x<grid[1] and grid.size-1 for x>grid[end].
"""
function Base.floor(grid::Uniform{T}, x) where {T}
    result = (x - grid.grid[1]) / (grid.grid[end] - grid.grid[1]) * (grid.size - 1) + 1
    if result < 1
        return 1
    elseif result >= grid.size
        return grid.size - 1
    else
        return Base.floor(Int, result)
    end

end
function Base.floor(grid::Uniform{T,PeriodicBound}, x) where {T}
    result = (x - grid.grid[1]) / (grid.grid[end] - grid.grid[1]) * (grid.size - 1)
    result = (result % grid.size + grid.size) % grid.size + 1
    if result < 1
        return grid.size
    elseif result >= grid.size
        return grid.size
    else
        return Base.floor(Int, result)
    end
end

"""
    struct BaryCheb{T<:AbstractFloat} <: AbstractGrid

BaryCheb grid generated on [bound[1], bound[2]] with order N.

#Members:
- `bound` : boundary of the grid
- `size` : number of grid points
- `grid` : grid points
- `weight` : interpolation weight

#Constructor:
-    function BaryCheb{T}(bound, size) where {T<:AbstractFloat}
"""
struct BaryCheb{T<:AbstractFloat} <: AbstractGrid{T}
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}
    weight::Vector{T}
    invVandermonde::Matrix{T}
    """
        function BaryCheb{T}(bound, N) where {T<:AbstractFloat}

    create BaryCheb grid.
    """
    function BaryCheb{T}(bound, N) where {T<:AbstractFloat}
        @assert bound[1] < bound[2]
        @assert N > 1
        order = N
        x, w = barychebinit(order)
        grid = zeros(T, N)
        a, b = bound[1], bound[2]
        weight = (b - a) / 2 .* w
        grid = (a + b) / 2 .+ (b - a) / 2 .* x
        invVandermonde = inv(transpose(vandermonde(x)))

        return new{T}(bound, N, grid, weight, invVandermonde)
    end
    function BaryCheb{T}(bound, N, invVandermonde) where {T<:AbstractFloat}
        # use given Vandermonde matrix, useful for composite grid that has many BaryCheb subgrids with same order
        order = N
        x, w = barychebinit(order)
        grid = zeros(T, N)
        a, b = bound[1], bound[2]
        weight = (b - a) / 2 .* w
        grid = (a + b) / 2 .+ (b - a) / 2 .* x
        return new{T}(bound, N, grid, weight, invVandermonde)
    end
end

"""
    struct GaussLegendre{T<:AbstractFloat} <: AbstractGrid

GaussLegendre grid generated on [bound[1], bound[2]] with order N.

#Members:
- `bound` : boundary of the grid
- `size` : number of grid points
- `grid` : grid points
- `weight` : integration weight

#Constructor:
-    function GaussLegendre{T}(bound, size) where {T<:AbstractFloat}
"""
struct GaussLegendre{T<:AbstractFloat} <: AbstractGrid{T}
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}
    weight::Vector{T}

    """
        function GaussLegendre{T}(bound, N) where {T<:AbstractFloat}

    create GaussLegendre grid.
    """
    function GaussLegendre{T}(bound, N) where {T<:AbstractFloat}
        @assert bound[1] < bound[2]
        @assert N > 1
        order = N
        x, w = gausslegendre(order)
        grid = zeros(T, N)
        a, b = bound[1], bound[2]
        weight = (b - a) / 2 .* w
        grid = (a + b) / 2 .+ (b - a) / 2 .* x

        return new{T}(bound, N, grid, weight)
    end
end

"""
    struct Log{T<:AbstractFloat} <: AbstractGrid

Log grid generated on [bound[1], bound[2]] with N grid points.
Minimal interval is set to be minterval. Dense to sparse if d2s, vice versa.

On [0, 1], a typical d2s Log grid looks like
[0, λ^(N-1), ..., λ^2, λ, 1].

#Members:
- `bound` : boundary of the grid
- `size` : number of grid points
- `grid` : grid points
- `weight` : integration weight

- `λ` : scale parameter
- `d2s` : dense to sparse or not

#Constructor:
-    function Log{T}(bound, size, minterval, d2s) where {T<:AbstractFloat}
"""
struct Log{T<:AbstractFloat} <: AbstractGrid{T}
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}
    weight::Vector{T}

    λ::T
    d2s::Bool

    """
        function Log{T}(bound, N, minterval, d2s) where {T<:AbstractFloat}

    create Log grid.
    """
    function Log{T}(bound, N, minterval, d2s) where {T<:AbstractFloat}
        @assert bound[1] < bound[2]
        @assert N > 1
        @assert 0 < minterval < bound[2] - bound[1]
        grid = zeros(T, N)
        M = N - 2
        λ = (minterval / (bound[2] - bound[1]))^(1.0 / M)

        if d2s
            for i in 1:M
                grid[i+1] = bound[1] + (bound[2] - bound[1]) * λ^(M + 1 - i)
            end
        else
            for i in 2:M+1
                grid[i] = bound[2] - (bound[2] - bound[1]) * λ^(i - 1)
            end
        end
        grid[1] = bound[1]
        grid[end] = bound[2]
        weight = similar(grid)
        for i in 1:N
            if i == 1
                weight[1] = 0.5 * (grid[2] - grid[1])
            elseif i == N
                weight[end] = 0.5 * (grid[end] - grid[end-1])
            else
                weight[i] = 0.5 * (grid[i+1] - grid[i-1])
            end
        end

        return new{T}(bound, N, grid, weight, λ, d2s)
    end
end

Log{T}(; bound, N, minterval, d2s) where {T<:AbstractFloat} = Log{T}(bound, N, minterval, d2s)

function Base.show(io::IO, grid::Log; isSimplified=false)
    if isSimplified
        print(io,
            "$(typeof(grid)): 1D " * (grid.d2s ? "dense to sparse" : "sparse to dense")
            * " Log Grid with $(grid.size) grid points.\n"
        )
    else
        print(io,
            "$(typeof(grid)): 1D " * (grid.d2s ? "dense to sparse" : "sparse to dense")
            * " Log Grid\n"
            * "- bound: $(grid.bound)\n"
            * "- size: $(grid.size)\n"
            * "- grid: $(grid.grid)\n"
        )
    end
end


"""
    function Base.floor(grid::Log{T}, x) where {T}

find the index of largest
grid point smaller than x.

return 1 for x<grid[1] and grid.size-1 for x>grid[end].
"""
function Base.floor(grid::Log{T}, x) where {T}
    if x <= grid.grid[1]
        return 1
    elseif x >= grid.grid[end]
        return grid.size - 1
    end

    a, b = grid.bound[1], grid.bound[2]
    if grid.d2s
        i = Base.floor(log((x - a) / (b - a)) / log(grid.λ))
        if i > grid.size - 2
            result = 1
        else
            result = grid.size - i - 1
        end
    else
        i = Base.floor(log((b - x) / (b - a)) / log(grid.λ))
        if i > grid.size - 2
            result = grid.size - 1
        else
            result = i + 1
        end
    end

    return Base.floor(Int, result)
end

@inline function denseindex(grid::Log)
    return [(grid.d2s) ? 1 : grid.size,]
end

BoundType(::Type{<:Arbitrary{T,BT}}) where {T,BT} = BT()
BoundType(::Type{<:Uniform{T,BT}}) where {T,BT} = BT()
BoundType(::Type{<:Union{BaryCheb,GaussLegendre}}) = OPENBOUND
BoundType(::Type{<:Log}) = CLOSEDBOUND

end
