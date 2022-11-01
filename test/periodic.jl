@testset "Periodic Bound Grids" begin
    rng = MersenneTwister(1453)

    @testset "SimpleGrids with PeriodicBound" begin
        @testset "Arbitrary" begin
            grid = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0, π, 6.0, 2π]
            ag = SimpleG.Arbitrary{Float64}(grid; boundtype=SimpleG.PERIODICBOUND)
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

            testx = rand(rng, 10) * 2π
            for x in testx
                @test f(x) ≈ Interp.interp1D(data, ag, x)
            end
        end

        @testset "Uniform" begin
            bound = [0, 2π]
            N = 10
            ag = SimpleG.Uniform{Float64}(bound, N; boundtype=SimpleG.PERIODICBOUND)
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

            testx = rand(rng, 10) * 2π
            for x in testx
                @test f(x) ≈ Interp.interp1D(data, ag, x)
            end
        end
    end

end