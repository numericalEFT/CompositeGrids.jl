@testset "MC" begin
    # testing locate and volumn functions provided for monte carlo
    rng = MersenneTwister(1234)

    function test_locate_volumn(grid; δ = 1e-6)
        vol = 0.0
        for (ti, t) in enumerate(grid.grid)
            # test locate
            @test ti == Interp.locate(grid, t)
            @test ti == Interp.locate(grid, (ti == 1) ? t : t - δ)
            @test ti == Interp.locate(grid, (ti == length(grid)) ? t : t + δ)
            vol += Interp.volumn(grid, ti)
        end

        @test vol ≈ Interp.volumn(grid)
    end

    @testset "Locate and Volumn" begin
        Ng, Ns = 10, 4
        β = π

        tgrid = SimpleGrid.Uniform{Float64}([0.0, β], Ng)
        test_locate_volumn(tgrid)

        tgrid = SimpleGrid.BaryCheb{Float64}([0.0, β], Ng)
        test_locate_volumn(tgrid)

        tgrid = SimpleGrid.Log{Float64}([0.0, β], Ng, 0.001, true)
        test_locate_volumn(tgrid)

        tgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, β],[β/2,], Ng, 0.001, Ns)
        test_locate_volumn(tgrid)

    end

    function test_mc_histogram(grid; Npp=1e6, f = (x -> x))
        # println("testing:$(typeof(grid))")
        Ng = length(grid)
        a, b = grid.bound[1], grid.bound[2]
        Nmc = Npp * Ng

        hist = zeros(Ng)
        for i in 1:Nmc
            x = rand(rng) * (b - a) + a
            ind = Interp.locate(grid, x)
            vol = Interp.volumn(grid, ind)
            hist[ind] += f(x) / vol * Interp.volumn(grid)
        end

        for (ti, t) in enumerate(grid.grid)
            vi = Interp.volumn(grid, ti)
            vall = Interp.volumn(grid)
            if ti != 1 && ti != Ng
                @test isapprox(hist[ti] / Nmc, f(t),
                               rtol=5 / sqrt(Nmc*vi/vall), atol = 5*vi/vall)
            else
                # edge points has extra error because grid point is not at center of interval
                @test isapprox(hist[ti] / Nmc, f(t),
                               rtol=5 / sqrt(Nmc*vi/vall), atol = 5*vi/vall)
            end
        end
    end

    @testset "MC histogram" begin
        Ng, Nl, Ns = 16, 8, 4
        β = π

        tgrid = SimpleGrid.Uniform{Float64}([0.0, β], Ng)
        test_mc_histogram(tgrid)

        tgrid = SimpleGrid.BaryCheb{Float64}([0.0, β], Nl)
        test_mc_histogram(tgrid)

        tgrid = SimpleGrid.Log{Float64}([0.0, β], Nl, 0.1, true)
        test_mc_histogram(tgrid)

        tgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, β],[β/2,], Nl, 0.1, Ns)
        test_mc_histogram(tgrid)

    end

end
