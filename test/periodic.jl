@testset "Periodic Bound Grids" begin
    rng = MersenneTwister(1453)

    @testset "SimpleGrids with PeriodicBound" begin
        @testset "Arbitrary" begin
            grid = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 3.13, π, 3.15, 4.0, 5.0, 6.0, 6.28, 2π]
            ag = SimpleG.Arbitrary(grid; isperiodic=true)
            println(ag.grid)
            @test ag.grid == grid[1:end-1]

            # test floor
            δ = 1e-6
            for (i, x) in enumerate(ag)
                @test floor(ag, x + δ) == i
                @test floor(ag, x - δ) == ((i + length(ag) - 2) % length(ag)) + 1
            end

            # interp
            f(x) = π - abs(x - π)
            data = f.(ag.grid)

            @test f(δ) ≈ Interp.interp1D(data, ag, δ)
            @test f(2π - δ) ≈ Interp.interp1D(data, ag, 2π - δ)
            testx = rand(rng, 10) * 2π
            for x in testx
                @test f(x) ≈ Interp.interp1D(data, ag, x)
            end

            # @test Interp.integrate1D(data, ag) ≈ π^2
            @test isapprox(Interp.integrate1D(data, ag), π^2, atol=1e-3)
        end

        @testset "Uniform" begin
            bound = [0, 2π]
            N = 10
            ag = SimpleG.Uniform(bound, N; isperiodic=true)
            println(ag.grid)

            # test floor
            δ = 1e-6
            for (i, x) in enumerate(ag)
                @test floor(ag, x + δ) == i
                @test floor(ag, x - δ) == ((i + length(ag) - 2) % length(ag)) + 1
            end

            # interp
            f(x) = π - abs(x - π)
            data = f.(ag.grid)

            @test f(δ) ≈ Interp.interp1D(data, ag, δ)
            @test f(2π - δ) ≈ Interp.interp1D(data, ag, 2π - δ)
            testx = rand(rng, 10) * 2π
            for x in testx
                @test f(x) ≈ Interp.interp1D(data, ag, x)
            end

            @test Interp.integrate1D(data, ag) ≈ π^2
        end

        @testset "Composite" begin
            ag = CompositeGrid.LogDensedGrid(:uniform, [0.0, 2π], [π,], 4, 0.01, 4; isperiodic=true)
            println(ag.grid)
            println(ag.panel.grid)

            # test floor
            δ = 1e-6
            for (i, x) in enumerate(ag)
                @test floor(ag, x + δ) == i
                @test floor(ag, x - δ) == ((i + length(ag) - 2) % length(ag)) + 1
            end

            # interp
            f(x) = π - abs(x - π)
            data = f.(ag.grid)

            @test f(δ) ≈ Interp.interp1D(data, ag, δ)
            @test f(2π - δ) ≈ Interp.interp1D(data, ag, 2π - δ)
            testx = rand(rng, 10) * 2π
            for x in testx
                @test f(x) ≈ Interp.interp1D(data, ag, x)
            end

            @test Interp.integrate1D(data, ag) ≈ π^2
        end
    end


end