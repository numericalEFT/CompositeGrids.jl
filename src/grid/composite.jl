module CompositeGrid

using StaticArrays, FastGaussQuadrature
include("simple.jl")
using .SimpleGrid

struct Composite{T<:AbstractFloat,PG,SG} <: AbstractGrid{T}
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
            @assert panel[i]==subgrids[i].bound[1]
            @assert panel[i+1]==subgrids[i].bound[2]
            if i == 1
                inits[i] = 1
                append!(grid, subgrids[i].grid)
            else
                if grid[end] == subgrids[i].grid[1]
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


end
