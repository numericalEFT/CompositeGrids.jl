"""
Provide interpolation and integration.
"""

module Interp

using StaticArrays, FastGaussQuadrature, CompositeGrids

#include("chebyshev.jl")

# include("simple.jl")
# using .SimpleGrid
# include("composite.jl")
# using .CompositeGrid

abstract type InterpStyle end
struct FloorInterp <: InterpStyle end
struct ChebInterp <: InterpStyle end
struct CompositeInterp <: InterpStyle end

InterpStyle(::Type) = FloorInterp()
InterpStyle(::Type{<:SimpleGrid.BaryCheb}) = ChebInterp()
InterpStyle(::Type{<:CompositeGrid.Composite}) = CompositeInterp()

"""
    function linear1D(data, xgrid, x)

linear interpolation of data(x)

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
- x: x
"""

@inline function linear1D(data, xgrid, x)

    xarray = xgrid.grid

    xi0,xi1 = 0,0
    if(x<=xarray[firstindex(xgrid)])
        xi0=1
        xi1=2
    elseif(x>=xarray[lastindex(xgrid)])
        xi0=lastindex(xgrid)-1
        xi1=xi0+1
    else
        xi0=floor(xgrid,x)
        xi1=xi0+1
    end

    dx0, dx1 = x - xarray[xi0], xarray[xi1] - x

    d0, d1 = data[xi0], data[xi1]

    g = d0 * dx1 + d1 * dx0

    gx = g / (dx0 + dx1) 
    return gx
end

"""
    function interp1D(data, xgrid, x)

linear interpolation of data(x)

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
- x: x
"""

function interp1D(data, xgrid::T, x) where {T}
    interp1D(InterpStyle(T), data, xgrid, x)
end

"""
    function interp1D(::FloorInterp,data, xgrid, x)

linear interpolation of data(x), use floor and linear1D

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
- x: x
"""

function interp1D(::FloorInterp,data, xgrid, x)
    return linear1D(data, xgrid, x)
end

"""
    function interp1D(::ChebInterp, data, xgrid, x)

linear interpolation of data(x), barycheb for BaryCheb grid

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
- x: x
"""

function interp1D(::ChebInterp, data, xgrid, x)
    return SimpleGrid.barycheb(xgrid.size, x, data, xgrid.weight, xgrid.grid)
end

"""
    function interp1D(::CompositeInterp,data, xgrid, x)

linear interpolation of data(x),
first floor on panel to find subgrid, then call interp1D on subgrid 

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
- x: x
"""

function interp1D(::CompositeInterp,data, xgrid, x)
    i = floor(xgrid.panel, x)
    head, tail = xgrid.inits[i], xgrid.inits[i]+xgrid.subgrids[i].size-1
    return interp1D(data[head:tail], xgrid.subgrids[i], x)
end


"""
    function interpGrid(data, xgrid, grid)

linear interpolation of data(grid[1:end]), return a Vector

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
- grid: points to be interpolated on
"""

function interpGrid(data, xgrid::T, grid) where {T}
    interpGrid(InterpStyle(T), data, xgrid, grid)
end

"""
    function interpGrid(::Union{FloorInterp,ChebInterp}, data, xgrid, grid)

linear interpolation of data(grid[1:end]), return a Vector
simply call interp1D on each points

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
- grid: points to be interpolated on
"""

function interpGrid(::Union{FloorInterp,ChebInterp}, data, xgrid, grid)
    ff = zeros(eltype(data), length(grid))
    for (xi, x) in enumerate(grid)
        ff[xi] = interp1D(data, xgrid, x)
    end
    return ff
end

"""
    function interpGrid(::CompositeInterp, data, xgrid, grid)

linear interpolation of data(grid[1:end]), return a Vector
grid should be sorted.

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
- grid: points to be interpolated on
"""

function interpGrid(::CompositeInterp, data, xgrid, grid)
    ff = zeros(eltype(data), length(grid))

    init, curr = 1, 1
    for pi in 1:xgrid.panel.size-1
        if grid[init]< xgrid.panel[pi+1]
            head, tail = xgrid.inits[pi], xgrid.inits[pi]+xgrid.subgrids[pi].size-1
            while grid[curr]<xgrid.panel[pi+1] && curr<length(grid)
                curr += 1
            end
            if grid[curr]<xgrid.panel[pi+1] && curr==length(grid)
                ff[init:curr] = interpGrid(data[head:tail], xgrid.subgrids[pi], grid[init:curr])
            else
                ff[init:curr-1] = interpGrid(data[head:tail], xgrid.subgrids[pi], grid[init:curr-1])
            end
            # println(data[head:tail])
            # println(xgrid.subgrids[pi].grid)
            # println(grid[init:curr-1])
            # println(ff[init:curr-1])
            init = curr
        end
    end
    return ff
end

abstract type IntegrateStyle end
struct WeightIntegrate <: IntegrateStyle end
struct NoIntegrate <: IntegrateStyle end
struct CompositeIntegrate <: IntegrateStyle end

IntegrateStyle(::Type) = NoIntegrate()
IntegrateStyle(::Type{<:SimpleGrid.GaussLegendre}) = WeightIntegrate()
IntegrateStyle(::Type{<:SimpleGrid.Uniform}) = WeightIntegrate()
IntegrateStyle(::Type{<:SimpleGrid.Arbitrary}) = WeightIntegrate()
IntegrateStyle(::Type{<:SimpleGrid.Log}) = WeightIntegrate()
IntegrateStyle(::Type{<:CompositeGrid.Composite}) = CompositeIntegrate()


"""
    function integrate1D(data, xgrid)

calculate integration of data[i] on xgrid

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
"""

function integrate1D(data, xgrid::T) where {T}
    return integrate1D(IntegrateStyle(T), data, xgrid)
end

"""
    function integrate1D(::NoIntegrate, data, xgrid)

calculate integration of data[i] on xgrid
works for grids that do not have integration weight stored

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
"""

function integrate1D(::NoIntegrate, data, xgrid)
    return 0.0
    result = eltype(data)(0.0)

    grid = xgrid.grid
    for i in 1:xgrid.size
        if i==1
            weight = 0.5*(grid[2]-grid[1])
        elseif i==xgrid.size
            weight = 0.5*(grid[end]-grid[end-1])
        else
            weight = 0.5*(grid[i+1]-grid[i-1])
        end
        result += data[i]*weight
    end
    return result
end

"""
    function integrate1D(::WeightIntegrate, data, xgrid)

calculate integration of data[i] on xgrid
works for grids that have integration weight stored

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
"""

function integrate1D(::WeightIntegrate, data, xgrid)
    result = eltype(data)(0.0)

    for i in 1:xgrid.size
        result += data[i]*xgrid.weight[i]
    end
    return result
end

"""
    function integrate1D(::CompositeIntegrate, data, xgrid)

calculate integration of data[i] on xgrid
call integrate1D for each subgrid and return the sum.

#Arguments:
- xgrid: one-dimensional grid of x
- data: one-dimensional array of data
"""

function integrate1D(::CompositeIntegrate, data, xgrid)
    result = eltype(data)(0.0)

    for pi in 1:xgrid.panel.size-1
        head, tail = xgrid.inits[pi], xgrid.inits[pi]+xgrid.subgrids[pi].size-1
        result += integrate1D( data[head:tail],xgrid.subgrids[pi])
        currgrid = xgrid.subgrids[pi]
    end
    return result

end

end
