module Interp

using StaticArrays, FastGaussQuadrature

include("chebyshev.jl")

include("simple.jl")
using .SimpleGrid
include("composite.jl")
using .CompositeGrid

"""
   linear1D(data,xgrid, x) 

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

function interp1D(data, xgrid, x)
    return linear1D(data, xgrid, x)
end

function interp1D(data, xgrid::BaryCheb, x)
    return barycheb(xgrid.size, x, data, xgrid.weight, xgrid.grid)
end

function interp1D(data, xgrid::Composite, x)
    i = floor(xgrid.panel, x)
    head, tail = xgrid.inits[i], xgrid.inits[i]+xgrid.subgrids[i].size-1
    return interp1D(data[head:tail], xgrid.subgrids[i], x)
end

function interpGrid(data, xgrid, grid)
    ff = zeros(eltype(data), length(grid))
    for (xi, x) in enumerate(grid)
        ff[xi] = interp1D(data, xgrid, x)
    end
    return ff
end


function interpGrid(data, xgrid::Composite, grid)
    ff = zeros(eltype(data), length(grid))

    init, curr = 1, 1
    for pi in xgrid.panel.size-1
        if grid[init]< xgrid.panel[pi+1]
            head, tail = xgrid.inits[pi], xgrid.inits[pi]+xgrid.subgrids[pi].size-1
            while grid[curr]<xgrid.panel[pi+1]
                curr += 1
            end
            ff[init:curr-1] = interpGrid(data[head:tail], xgrid.subgrids[pi], grid[init:curr-1])
            init = curr
        end
    end
    return ff
end

end
