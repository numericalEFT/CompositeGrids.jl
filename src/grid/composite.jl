module CompositeGrid

using StaticArrays, FastGaussQuadrature
include("simple.jl")
using .SimpleGrid

struct Composite{T<:AbstractFloat,PG,SG} <: ClosedGrid
    bound::SVector{2,T}
    size::Int
    grid::Vector{T}


    panel::PG
    subgrids::Vector{SG}
    inits::Vector{Int}

    function Composite{T,PG,SG}(panel, subgrids) where {T<:AbstractFloat,PG,SG}
        bound = [panel[1], panel[end]]
        @assert panel.size-1==length(subgrids)
        inits = zeros(Int, length(subgrids))
        grid = Vector{T}([])
        for i in 1:length(subgrids)
            @assert panel[i]==subgrids[i].bound[1] "$(panel[i])!=$(subgrids[i].bound[1])"
            @assert panel[i+1]==subgrids[i].bound[2] "$(panel[i+1])!=$(subgrids[i].bound[2])"
            if i == 1
                inits[i] = 1
                append!(grid, subgrids[i].grid)
            else
                if abs(grid[end] - subgrids[i].grid[1]) < eps(T)*10
                    inits[i] = length(grid)
                    append!(grid, subgrids[i].grid[2:end])
                else
                    inits[i] = length(grid)+1
                    append!(grid, subgrids[i].grid)
                end
            end

        end
        size = length(grid)
        
        return new{T,PG,SG}(bound, size, grid, panel, subgrids,inits)
    end

end

function Base.floor(grid::Composite{T,PG,SG}, x) where {T,PG,SG}
    if SG<:ClosedGrid
        i = floor(grid.panel, x)
        return grid.inits[i]-1+floor(grid.subgrids[i],x)
    end
    
    if x <= grid.grid[1]
        return 1
    elseif x >= grid.grid[end]
        return grid.size-1
    end

    result = searchsortedfirst(grid.grid, x)-1
    return result
end

Base.getindex(grid::Composite, i) = grid.grid[i]
Base.firstindex(grid::Composite) = 1
Base.lastindex(grid::Composite) = grid.size

function CompositeLogGrid(type, bound, N, minterval, d2s, order, T=Float64)
    if type == :cheb
        SubGridType = BaryCheb{T}
    elseif type == :gauss
        SubGridType = GaussLegendre{T}
    elseif type == :uniform
        SubGridType = Uniform{T}
    else
        error("$type not implemented!")
    end

    panel = Log{T}(bound, N, minterval, d2s)
    println("logpanel:",panel.grid)
    subgrids = Vector{SubGridType}([])

    for i in 1:N-1
        _bound = [panel[i],panel[i+1]]
        push!(subgrids, SubGridType(_bound,order))
    end

    return Composite{T, Log{T},SubGridType}(panel,subgrids)

end

function LogDensedGrid(type, bound, dense_on, N, minterval, order, T=Float64)
    if type == :cheb
        SubGridType = BaryCheb{T}
    elseif type == :gauss
        SubGridType = GaussLegendre{T}
    elseif type == :uniform
        SubGridType = Uniform{T}
    else
        error("$type not implemented!")
    end

    dense_on = sort(dense_on)
    @assert bound[1]<dense_on[1]<dense_on[end]<bound[2]
    dp = Vector{T}([])
    for i in 1:length(dense_on)
        if i==1
            if abs(dense_on[i]-bound[1])<minterval
                push!(dp, bound[1])
            else
                push!(dp, dense_on[i])
            end
        elseif i != length(dense_on)
            if abs(dense_on[i]-dp[end])<minterval
                if dp[end] != bound[1]
                    dp[end] = (dense_on[i]+dense_on[i-1])/2.0
                end
            else
                push!(dp, dense_on[i])
            end
        else
            if abs(dense_on[i]-bound[2])<minterval
                if abs(dp[end]-bound[2])<minterval
                    dp[end]=bound[2]
                else
                    push!(dp, bound[2])
                end
            elseif abs(dense_on[i]-dp[end])<minterval
                if dp[end] != bound[1]
                    dp[end] = (dense_on[i]+dense_on[i-1])/2.0
                end
            else
                push!(dp, dense_on[i])
            end
        end
    end

    panelgrid = Vector{T}([])
    d2slist = Vector{Bool}([])
    for i in 1:length(dp)
        if i==1
            push!(panelgrid, bound[1])
            if dp[1] != bound[1]
                push!(panelgrid, dp[1])
                push!(d2slist, false)
            end
        else
            push!(panelgrid, (dp[i]+dp[i-1])/2.0)
            push!(d2slist, !d2slist[end])
            push!(panelgrid, dp[i])
            push!(d2slist, !d2slist[end])
        end
    end
    if dp[end] !=  bound[2]
        push!(panelgrid, bound[2])
        @assert !d2slist[end] == true
        push!(d2slist, true)
    end

    panel = Arbitrary{T}(panelgrid)
    println("panel:",panel.grid)
    subgrids = Vector{Composite{T, Log{T}, SubGridType}}([])

    for i in 1:length(panel.grid)-1
        push!(subgrids, CompositeLogGrid(type, [panel[i],panel[i+1]], N, minterval, d2slist[i], order))
    end

    return Composite{T, Arbitrary{T},Composite{T, Log{T}, SubGridType}}(panel,subgrids)

end

end
