[![img](https://img.shields.io/badge/docs-dev-blue.svg)](https://numericaleft.github.io/CompositeGrids.jl/dev/)
[![img](https://github.com/numericaleft/CompositeGrids.jl/workflows/CI/badge.svg)](https://github.com/numericaleft/CompositeGrids.jl/actions)
[![img](https://codecov.io/gh/numericalEFT/CompositeGrids.jl/branch/main/graph/badge.svg?token=WN6HO1XASY)](https://codecov.io/gh/numericaleft/CompositeGrids.jl)


# Introduction

``CompositeGrids.jl`` is a powerful Julia package that provides a unified interface for generating a wide range of common 1D grids. In addition to basic grids, this package allows you to create composite grids by combining multiple fundamental grids. These composite grids are enriched with essential functionalities, including floor function, interpolation function, and integration function, all of which are optimized to enhance their performance on specific grids.

With ``CompositeGrids.jl``, you can effortlessly construct complex grids and efficiently handle various numerical tasks with ease. Whether you are working on scientific simulations, data analysis, or any other domain that involves grid-based calculations, this package will be your go-to tool for managing grids effectively.


# Installation
To install ``CompositeGrids.jl``, use Julia's package manager. Open the Julia REPL, type `]` to enter the package mode, and then:
```
pkg> add CompositeGrids.jl
```

# Quick Start

In this quick start example, we will demonstrate how to generate a grid from 0 to 1, log-densed at 0 and 1, and optimized for integration using the ``CompositeGrids.jl`` package. We will provide descriptive comments in the code to guide you through the process.

```julia
    using CompositeGrids
    
    # Generating a log densed composite grid with LogDensedGrid()
    tgrid = CompositeGrid.LogDensedGrid(
        type=:gauss,# The top layer grid is :gauss, optimized for integration. For interpolation use :cheb
        bound=[0.0, 1],# The grid is defined on [0.0, β]
        dense_at=[0.0, 1],# and is densed at 0.0 and β, as given by 2nd and 3rd parameter.
        N=5,# N of log grid
        minterval=0.005, # minimum interval length of log grid
        order=5 # N of bottom layer
    )
    # The grid has 3 layers.
    # The top layer is defined by the boundary and densed points. In this case its:
    println("Top layer:",tgrid.panel.grid)
    # The middle layer is a log grid with 4 points and minimum interval length 0.001:
    println("First subgrid of middle layer:",tgrid.subgrids[1].panel.grid)
    # The bottom layer is a Gauss-Legendre grid with 5 points:
    println("First subgrid of bottom layer:",tgrid.subgrids[1].subgrids[1].grid)
    
    # function to be integrated:
    f(t) = exp(t)+exp(1-t)
    # numerical value on grid points:
    data = [f(t) for (ti, t) in enumerate(tgrid.grid)]
    
    # integrate with integrate1D():
    int_result = Interp.integrate1D(data, tgrid)
    
    println("result=",int_result)
    println("comparing to:",2*(exp(1)-1))
```

```
Top layer:[0.0, 0.5, 1.0]
First subgrid of middle layer:[0.0, 0.005000000000000001, 0.023207944168063897, 0.1077217345015942, 0.5]
First subgrid of bottom layer:[0.00023455038515334025, 0.0011538267247357924, 0.0025000000000000005, 0.0038461732752642086, 0.004765449614846661]
result=3.43656365691809
comparing to:3.43656365691809
```

# Manual

## Basics

The ``CompositeGrids.jl`` package offers two modules for working with 1D grids: ``SimpleGrid`` and ``CompositeGrid``. These modules provide a collection of common 1D grids with straightforward definitions and simple structures. Additionally, ``CompositeGrid`` defines a more general type of grid, composed of a panel grid and a set of subgrids, allowing for more flexibility in grid construction.

The common interface for grids includes the following properties and methods:

 -   ``g.bound``: This property gives the boundary of the interval of the grid. It provides a clear indication of the range covered by the grid.

 -   ``g.size``: This property gives the total number of grid points in the grid. It helps to determine the grid's resolution and granularity.

 -   ``g.grid``: This property gives an array of grid points. The array contains the coordinates of all the grid points within the specified boundary.

 -   ``g[i]``: This method returns the i-th grid point, which is the same as ``g.grid[i]``. It allows for direct access to specific grid points.

 -   ``floor(g, x)``: This method returns the largest index of the grid point where ``g[i] < x``. For values of x below the first grid point, it returns 1, and for values greater than the last grid point, it returns ``(grid.size - 1)``. This ensures that both ``floor()`` and ``(floor() + 1)`` are valid grid indices for any value of x.

The ``CompositeGrids.jl`` package also provides interpolation and integration functionalities for the grids. Different implementations are available for different types of grids, allowing for efficient numerical calculations tailored to each grid type.

## Simple Grids

The `SimpleGrid` module in `CompositeGrids.jl` offers various basic grids that serve as standalone grids and components of composite grids. The available basic grids include:

- **Arbitrary Grid:** The most general basic grid, which takes an array and converts it into a grid. It provides an efficient O(ln(N)) floor function based on `searchsortedfirst()`.

- **Uniform Grid:** Defined by the boundary and the number of grid points. It offers an O(1) floor function for rapid point location.

- **Log Grid:** Defined by the boundary, number of grid points, minimum interval, and direction. It generates a log-dense grid based on the provided parameters. An O(1) floor function is provided. For example:
```julia
    using CompositeGrids
    loggrid = SimpleGrid.Log{Float64}([0.0,1.0], 6, 0.0001, true)
    println(loggrid.grid)
```
```
    [0.0, 0.00010000000000000005, 0.0010000000000000002, 0.010000000000000002, 0.1, 1.0]
```

- **BaryCheb Grid:** Specifically designed for interpolation, it is defined by the boundary and number of grid points. The grid points are distributed according to Chebyshev nodes. The floor function is not optimized, so the O(ln(N)) function will be used, but the interpolation is based on an high precision algorithm with O(N).

- **GaussLegendre Grid:** Tailored for integration purposes, it is defined by the boundary and number of grid points. The grid points are distributed according to Gauss-Legendre quadrature. The floor function is not optimized, so the O(ln(N)) function will be used. The 1D integration is optimized.

It's important to note that grids can be categorized into open grids and closed grids. Closed grids indicate that the boundary points are also included as grid points, while open grids exclude the boundary points. The ``BaryCheb`` and `GaussLegendre` grids are examples of open grids.

A detailed manual can be found [here](https://numericaleft.github.io/CompositeGrids.jl/dev/lib/simple/).

## Composite Grids

The `CompositeGrid` module in `CompositeGrids.jl` provides a general type of grid where the entire interval is first divided by a panel grid, and then each interval of the panel grid is further divided by smaller grids called subgrids. Notably, subgrids can also be composite grids themselves, allowing for hierarchical grid structures.

The `LogDensedGrid` is a particularly useful generator of composite grids that offers a general solution when a log-dense 1D grid is needed around specific points within an interval. For example, grids like &tau; grids, which require densification around 0 and &beta;, or momentum grids, which need densification around the Fermi momentum, can be efficiently generated using `LogDensedGrid`.

The `LogDensedGrid` is defined as a three-layer composite grid:

1. **Top Layer (Arbitrary Grid):** This layer is defined by the boundary and the points where the grid needs to be dense. It allows for high flexibility in defining the grid structure.

2. **Middle Layer (Log Grid):** This layer is log-dense at the specified points, ensuring finer resolution in regions of interest.

3. **Bottom Layer (Grid of Options):** The bottom layer can be one of three options: `:cheb` for BaryCheb grid used for interpolation, `:gauss` for GaussLegendre grid used for integration, and `:uniform` for a uniform grid that serves general purposes.

The floor function of the composite grid is defined recursively. When locating a grid point, the floor function of the panel grid is called first to find the corresponding subgrid, and then the floor function of the subgrid is called to determine the final result. Since subgrids can themselves be composite grids, this recursive process continues until the lowest level of subgrids is reached.

The hierarchical nature of composite grids allows for the creation of sophisticated grid structures tailored to specific needs, making them a powerful tool for various scientific and computational applications.

A detailed manual can be found [here](https://numericaleft.github.io/CompositeGrids.jl/dev/lib/composite/).


## Interpolation and Integration

### Interpolation

Interpolation in `CompositeGrids.jl` provides an estimate of the function value at a given point `x` using the provided grid and function values on the grid points. For most of the simple grids, linear interpolation is used in conjunction with the floor function to locate the corresponding grid points. Notably, the BaryCheb grid employs an optimized algorithm for interpolation, leveraging information from all grid points to yield more precise results with the same number of grid points. This enhanced interpolation is subject to the condition that the function itself is smooth enough. When working with composite grids, the interpolation process is performed recursively, and the final result depends on the type of the lowest-level grid. For higher dimensions where data is defined on a list of grids, linear interpolation is provided, even when some of the grids are of the BaryCheb type.

### Integration

Integration over 1D grids is supported in `CompositeGrids.jl`. For most of the simple grids, linear integration is employed. However, for ``GaussLegendre`` grids and ``BaryCheb`` grids, an optimized integration method is used. Similar to interpolation, integration for composite grids is also carried out recursively, and the chosen method depends on the type of the lowest-level grids.

### Differentiation

The `CompositeGrids.jl` package offers differentiation for 1D grids. Specifically, a high-precision algorithm is implemented for ``BaryCheb`` grids, resulting in accurate differentiation results.

A detailed manual can be found [here](https://numericaleft.github.io/CompositeGrids.jl/dev/lib/interpolate/).

An example of interpolation and differenciation is shown below:
```julia
using CompositeGrids
β = π

# Generating a log densed composite grid with LogDensedGrid()
tgrid = CompositeGrid.LogDensedGrid(
    type=:cheb,# The top layer grid is :cheb
    bound=[0.0, β],# The grid is defined on [0.0, β]
    dense_at=[0.0, β],# and is densed at 0.0 and β, as given by 2nd and 3rd parameter.
    N=5,# N of log grid
    minterval=0.005, # minimum interval length of log grid
    order=5 # N of bottom layer
)

# function to be represented:
f(t) = sin(t)
# numerical value on grid points:
data = [f(t) for t in tgrid]

# integrate with integrate1D():
sin1 = Interp.interp1D(data, tgrid, 1.0)
dsin1 = Interp.differentiate1D(data, tgrid, 1.0)

println("result=", (sin1, dsin1))
println("comparing to:", (sin(1.0), cos(1.0)))
```
```
result=(0.8414425112056995, 0.5400742649805592)
comparing to:(0.8414709848078965, 0.5403023058681398)
```

## Complexity of Operations

| Grid             | Floor | Interpolate | Integrate | Differentiate |
|------------------|------------------------|------------------------|----------------------|---------------------------|
| Arbitrary        | O(ln(N))  | O(ln(N))      | O(N)   | O(ln(N))  |
| Uniform          | O(1)   | O(1)       | O(N)  | O(1)   |
| Log              | O(1)   | O(1)       | O(N)     | O(1)          |
| BaryCheb         | O(ln(N))    | O(N)              | O(N)    | O(N)                |
| GaussLegendre    | O(ln(N))      | O(ln(N))          | O(N)      | O(ln(N))     |
| CompositeGrid    | O(floor(panel)*floor(sub)) | O(floor(panel)*interp(sub)) | O(N) | O(floor(panel)*diff(sub))  |

**Note 1**: For `CompositeGrid`, the complexity depends on the type of the low-level grid in the hierarchy.
**Note 2**: Interpolation and differentiation of ``BaryCheb`` are implemented with high-precision algorithm, thus much less grid points are needed despite complexity.
