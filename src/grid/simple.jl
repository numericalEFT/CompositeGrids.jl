module SimpleGrid

export AbstractGrid, OpenGrid, ClosedGrid, Uniform, BaryCheb, GaussLegendre, Arbitrary, Log

using StaticArrays, FastGaussQuadrature

include("chebyshev.jl")

abstract type AbstractGrid end
abstract type OpenGrid <: AbstractGrid end
abstract type ClosedGrid <: AbstractGrid end

struct Arbitrary{T<:AbstractFloat} <: ClosedGrid
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}
    weight::Vector{T}

    function Arbitrary{T}(grid) where {T<:AbstractFloat}
        bound = [grid[1],grid[end]]
        size = length(grid)
        weight = similar(grid)
        for i in 1:size
            if i==1
                weight[1] = 0.5*(grid[2]-grid[1])
            elseif i==size
                weight[end] = 0.5*(grid[end]-grid[end-1])
            else
                weight[i] = 0.5*(grid[i+1]-grid[i-1])
            end
        end
        return new{T}(bound, size, grid, weight)
    end
end

function Base.floor(grid::AbstractGrid, x) #where {T}
    if x <= grid.grid[1]
        return 1
    elseif x >= grid.grid[end]
        return grid.size-1
    end

    result = searchsortedfirst(grid.grid, x)-1
    return Base.floor(Int, result)
end

Base.getindex(grid::Arbitrary, i) = grid.grid[i]
Base.firstindex(grid::Arbitrary) = 1
Base.lastindex(grid::Arbitrary) = grid.size

struct Uniform{T<:AbstractFloat} <: ClosedGrid
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}
    weight::Vector{T}

    function Uniform{T}(bound, N) where {T<:AbstractFloat}
        Ntot = N - 1
        interval = (bound[2]-bound[1])/Ntot
        grid = bound[1] .+ Vector(1:N) .* interval .- ( interval )
        weight = similar(grid)
        for i in 1:N
            if i==1
                weight[1] = 0.5*(grid[2]-grid[1])
            elseif i==N
                weight[end] = 0.5*(grid[end]-grid[end-1])
            else
                weight[i] = 0.5*(grid[i+1]-grid[i-1])
            end
        end
        return new{T}(bound, N, grid, weight)
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

struct BaryCheb{T<:AbstractFloat} <: OpenGrid
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

Base.getindex(grid::BaryCheb, i) = grid.grid[i]
Base.firstindex(grid::BaryCheb) = 1
Base.lastindex(grid::BaryCheb) = grid.size

struct GaussLegendre{T<:AbstractFloat} <: OpenGrid
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

Base.getindex(grid::GaussLegendre, i) = grid.grid[i]
Base.firstindex(grid::GaussLegendre) = 1
Base.lastindex(grid::GaussLegendre) = grid.size

struct Log{T<:AbstractFloat} <: ClosedGrid
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}
    weight::Vector{T}

    λ::T
    d2s::Bool

    function Log{T}(bound, N, minterval, d2s) where {T<:AbstractFloat}
        grid = zeros(T, N)
        M = N-2
        λ = (minterval/(bound[2]-bound[1]))^(1.0/M)

        if d2s
            for i in 1:M
                grid[i+1] = bound[1] + (bound[2]-bound[1])*λ^(M+1-i)
            end
        else
            for i in 2:M+1
                grid[i] = bound[2] - (bound[2]-bound[1])*λ^(i-1)
            end
        end
        grid[1] = bound[1]
        grid[end] = bound[2]
        weight = similar(grid)
        for i in 1:N
            if i==1
                weight[1] = 0.5*(grid[2]-grid[1])
            elseif i==N
                weight[end] = 0.5*(grid[end]-grid[end-1])
            else
                weight[i] = 0.5*(grid[i+1]-grid[i-1])
            end
        end

        return new{T}(bound, N, grid,weight, λ, d2s)
    end
end

function Base.floor(grid::Log{T}, x) where {T}
    if x <= grid.grid[1]
        return 1
    elseif x >= grid.grid[end]
        return grid.size-1
    end

    a,b=grid.bound[1],grid.bound[2]
    if grid.d2s
        i = Base.floor( log((x-a)/(b-a))/log(grid.λ)  )
        if i > grid.size-2
            result = 1
        else
            result = grid.size-i-1
        end
    else
        i = Base.floor( log((b-x)/(b-a))/log(grid.λ)  )
        if i > grid.size-2
            result = grid.size-1
        else
            result = i+1
        end
    end

    return Base.floor(Int, result)
end

Base.getindex(grid::Log, i) = grid.grid[i]
Base.firstindex(grid::Log) = 1
Base.lastindex(grid::Log) = grid.size


end
