@testset "SimpleGrids" begin

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

    @testset "ArbitraryGrid" begin

        # resulted grid should be sorted
        println("Testing Arbitrary constructor, a warning should appear.")
        arbitrary = SimpleGrid.Arbitrary{Float64}([0.0, 0.2, 0.6, 0.8, 0.4, 1.0])
        @test arbitrary.grid == [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

        uniform = SimpleGrid.Uniform{Float64}([0.0, 1.0], 10)
        arbitrary = SimpleGrid.Arbitrary{Float64}(uniform.grid)
        println(arbitrary)
        @test floor(arbitrary, 0.0) == 1
        @test floor(arbitrary, arbitrary[1]) == 1

        δ = 1.0e-12
        check(arbitrary, 2:arbitrary.size-1, δ, 0)
        check(arbitrary, 2:arbitrary.size-1, -δ, -1)

        @test floor(arbitrary, arbitrary[end]) == arbitrary.size - 1
        @test floor(arbitrary, 1.0) == arbitrary.size - 1

        uniform = SimpleGrid.Uniform{Float64}([0.0, 1.0], 2)
        arbitrary = SimpleGrid.Arbitrary{Float64}(uniform.grid)
        @test floor(arbitrary, 0.0) == 1
        @test floor(arbitrary, arbitrary[1]) == 1

        δ = 1.0e-12
        check(arbitrary, 2:arbitrary.size-1, δ, 0)
        check(arbitrary, 2:arbitrary.size-1, -δ, -1)

        @test floor(arbitrary, arbitrary[end]) == arbitrary.size - 1
        @test floor(arbitrary, 1.0) == arbitrary.size - 1

        arb2 = SimpleGrid.Arbitrary{Float64}([1.0,])
        @test floor(arb2, 0.0) == 1
        @test floor(arb2, 1.0) == 1
        @test floor(arb2, 2.0) == 1
    end

    @testset "UniformGrid" begin
        uniform = SimpleGrid.Uniform{Float64}([0.0, 1.0], 10)
        @test floor(uniform, 0.0) == 1
        @test floor(uniform, uniform[1]) == 1

        δ = 1.0e-12
        check(uniform, 2:uniform.size-1, δ, 0)
        check(uniform, 2:uniform.size-1, -δ, -1)

        @test floor(uniform, uniform[end]) == uniform.size - 1
        @test floor(uniform, 1.0) == uniform.size - 1

        #uniform = Grid.Uniform{Float64,2}(0.0, 1.0, (true, true))
        uniform = SimpleGrid.Uniform{Float64}([0.0, 1.0], 2)
        @test floor(uniform, 0.0) == 1
        @test floor(uniform, uniform[1]) == 1

        δ = 1.0e-12
        check(uniform, 2:uniform.size-1, δ, 0)
        check(uniform, 2:uniform.size-1, -δ, -1)

        @test floor(uniform, uniform[end]) == uniform.size - 1
        @test floor(uniform, 1.0) == uniform.size - 1
    end

    @testset "BaryChebGrid" begin
        cheb = SimpleGrid.BaryCheb{Float64}([0.0, 1.0], 10)
        println(cheb.grid)
        @test floor(cheb, 0.0) == 1
        @test floor(cheb, cheb[1]) == 1

        δ = 1.0e-12
        check(cheb, 2:cheb.size-1, δ, 0)
        check(cheb, 2:cheb.size-1, -δ, -1)

        @test floor(cheb, cheb[end]) == cheb.size - 1
        @test floor(cheb, 1.0) == cheb.size - 1

        #cheb = Grid.BaryCheb{Float64,2}(0.0, 1.0, (true, true))
        cheb = SimpleGrid.BaryCheb{Float64}([0.0, 1.0], 2)
        @test floor(cheb, 0.0) == 1
        @test floor(cheb, cheb[1]) == 1

        δ = 1.0e-12
        check(cheb, 2:cheb.size-1, δ, 0)
        check(cheb, 2:cheb.size-1, -δ, -1)

        @test floor(cheb, cheb[end]) == cheb.size - 1
        @test floor(cheb, 1.0) == cheb.size - 1
    end

    @testset "GaussLegendreGrid" begin
        gauss = SimpleGrid.GaussLegendre{Float64}([0.0, 1.0], 10)
        # println(gauss.grid)
        @test floor(gauss, 0.0) == 1
        @test floor(gauss, gauss[1]) == 1

        δ = 1.0e-12
        check(gauss, 2:gauss.size-1, δ, 0)
        check(gauss, 2:gauss.size-1, -δ, -1)

        @test floor(gauss, gauss[end]) == gauss.size - 1
        @test floor(gauss, 1.0) == gauss.size - 1

        #gauss = Grid.GaussLegendre{Float64,2}(0.0, 1.0, (true, true))
        gauss = SimpleGrid.GaussLegendre{Float64}([0.0, 1.0], 2)
        @test floor(gauss, 0.0) == 1
        @test floor(gauss, gauss[1]) == 1

        δ = 1.0e-12
        check(gauss, 2:gauss.size-1, δ, 0)
        check(gauss, 2:gauss.size-1, -δ, -1)

        @test floor(gauss, gauss[end]) == gauss.size - 1
        @test floor(gauss, 1.0) == gauss.size - 1
    end

    @testset "LogGrid" begin
        # loggrid = SimpleGrid.Log{Float64}([0.0, 1.0], 6, 0.001, true)
        loggrid = SimpleGrid.Log{Float64}(bound=[0.0, 1.0], N=6, minterval=0.001, d2s=true)
        println(loggrid)
        # println(SimpleGrid.denseindex(loggrid))
        @test floor(loggrid, 0.0) == 1
        @test floor(loggrid, loggrid[1]) == 1

        δ = 1.0e-12
        check(loggrid, 2:loggrid.size-1, δ, 0)
        check(loggrid, 2:loggrid.size-1, -δ, -1)

        @test floor(loggrid, loggrid[end]) == loggrid.size - 1
        @test floor(loggrid, 1.0) == loggrid.size - 1

        loggrid = SimpleGrid.Log{Float64}([0.0, 1.0], 6, 0.001, false)
        println(loggrid)
        # println(SimpleGrid.denseindex(loggrid))
        @test floor(loggrid, 0.0) == 1
        @test floor(loggrid, loggrid[1]) == 1

        δ = 1.0e-12
        check(loggrid, 2:loggrid.size-1, δ, 0)
        check(loggrid, 2:loggrid.size-1, -δ, -1)

        @test floor(loggrid, loggrid[end]) == loggrid.size - 1
        @test floor(loggrid, 1.0) == loggrid.size - 1
    end

end
