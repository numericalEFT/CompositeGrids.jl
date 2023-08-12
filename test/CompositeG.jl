
@testset "CompositeGrids" begin

    # with a shift to the grid element, check if produce the correct floored index
    function check(grid, range, shift, idx_shift)
        for i in range
            @test(floor(grid, grid[i] + shift) == i + idx_shift)
            # if floor(grid, grid[i] + shift) != i + idx_shift
            #     return false
            # end
        end
        return true
    end

    @testset "Composite" begin
        uniform = SimpleGrid.Uniform{Float64}([0.0, 1.0], 3)
        gauss1 = SimpleGrid.GaussLegendre{Float64}([0.0, 0.5], 4)
        gauss2 = SimpleGrid.GaussLegendre{Float64}([0.5, 1.0], 4)
        comp = CompositeGrid.Composite{
            Float64,
            typeof(uniform),
            SimpleGrid.GaussLegendre{Float64}
        }(uniform, [gauss1, gauss2])

        println(comp)
        # println(comp.inits)

        @test floor(comp, 0.0) == 1
        @test floor(comp, comp[1]) == 1

        δ = 1.0e-12
        check(comp, 2:comp.size-1, δ, 0)
        check(comp, 2:comp.size-1, -δ, -1)

        @test floor(comp, comp[end]) == comp.size - 1
        @test floor(comp, 1.0) == comp.size - 1

    end

    @testset "CompositeLog" begin

        comp = CompositeGrid.CompositeLogGrid(:cheb, [0.0, 1.0], 4, 0.001, true, 4)
        println(comp)
        # println(comp.inits)

        @test floor(comp, 0.0) == 1
        @test floor(comp, comp[1]) == 1

        δ = 1.0e-12
        check(comp, 2:comp.size-1, δ, 0)
        check(comp, 2:comp.size-1, -δ, -1)

        @test floor(comp, comp[end]) == comp.size - 1
        @test floor(comp, 1.0) == comp.size - 1

        # comp = CompositeGrid.LogDensedGrid(:uniform, [0.0, 10.0], [0.0, 1.0, 1.0, 2.0, 2.000001], 4, 0.001, 4)
        comp = CompositeGrid.LogDensedGrid(type=:uniform,
            bound=[0.0, 10.0],
            dense_at=[0.0, 1.0, 1.0, 2.0, 2.000001],
            N=4,
            minterval=0.001,
            order=4)
        # println(comp.grid)
        # println(comp.inits)
        # println(CompositeGrid.denseindex(comp))
        # println([comp[i] for i in CompositeGrid.denseindex(comp)])

        @test floor(comp, 0.0) == 1
        @test floor(comp, comp[1]) == 1

        δ = 1.0e-12
        check(comp, 2:comp.size-1, δ, 0)
        check(comp, 2:comp.size-1, -δ, -1)

        @test floor(comp, comp[end]) == comp.size - 1
        @test floor(comp, 10.0) == comp.size - 1

        # test grid generation
        comp = CompositeGrid.LogDensedGrid(:uniform, [0.0, 10.0], [0.0,], 4, 0.001, 4)
        # println(comp.grid)
        comp = CompositeGrid.LogDensedGrid(:uniform, [0.0, 10.0], [0.0, 1.0], 4, 0.001, 4)
        # println(comp.grid)
        comp = CompositeGrid.LogDensedGrid(:uniform, [0.0, 10.0], [0.5, 1.0], 4, 0.001, 4)
        # println(comp.grid)
        comp = CompositeGrid.LogDensedGrid(:uniform, [0.0, 10.0], [0.5, 10.0], 4, 0.001, 4)
        # println(comp.grid)

        # test edge case when two panel points are exactly minterval away
        # raw panel will be [0.0, 5.0,5.001,5.002,10.0], 
        # after construction it will be [0.0, 5.001, 10.0]
        comp = CompositeGrid.LogDensedGrid(:uniform, [0.0, 10.0], [5.0, 5.002], 4, 0.001, 4)
        @test length(comp.panel) == 3
        @test isapprox(comp.panel[2], 5.001, atol=0.001)
    end
    @testset "Issue #47" begin
        function no_non_increasing(vec::Vector{T}) where {T}
            for i in 2:length(vec)
                if vec[i] <= vec[i-1]
                    # Print the neighbors
                    if i == 2
                        println("Non-increasing point: ", vec[i], " with left neighbor: ", vec[i-1])
                        return false
                    elseif i == length(vec)
                        println("Non-increasing point: ", vec[i-1], " with right neighbor: ", vec[i])
                        return false
                    else
                        println("Non-increasing point: ", vec[i-1], " with neighbors: ", vec[i-2], " and ", vec[i])
                        return false
                    end
                end
            end
            return true
        end
        kF = 1.9191582926775128
        kgrid = CompositeGrid.LogDensedGrid(:uniform, [1.0 * kF, 100 * kF], [kF,], 8, 0.1 * kF, 8)
        @test no_non_increasing(kgrid.grid)
        # println(kgrid.panel.grid)
        # loggrid = kgrid.subgrids[1].panel
        # println(loggrid.grid)
        # ugrid1, ugrid2 = kgrid.subgrids[1].subgrids[6], kgrid.subgrids[1].subgrids[7]
        # println(ugrid1.grid)
        # println(ugrid2.grid)
        # println(ugrid1.bound)
        # println(ugrid2.bound)
    end
end

