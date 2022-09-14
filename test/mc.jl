@testset "MC" begin
    # testing locate and volumn functions provided for monte carlo
    rng = MersenneTwister(1234)

    @testset "Locate and Volumn" begin
        Ng = 33
        β = π
        tgrid = SimpleGrid.Uniform{Float64}([0.0, β], Ng)
        δ = 1e-6

        vol = 0.0
        for (ti, t) in enumerate(tgrid.grid)
            # test locate
            @test ti == Interp.locate(tgrid, t)
            @test ti == Interp.locate(tgrid, (ti == 1) ? t : t - δ)
            @test ti == Interp.locate(tgrid, (ti == length(tgrid)) ? t : t + δ)
            vol += Interp.volumn(tgrid, ti)
        end

        @test vol == Interp.volumn(tgrid)

    end

    @testset "MC histogram" begin
        Nmc = 1e6
        Ng = 33
        β = π
        tgrid = SimpleGrid.Uniform{Float64}([0.0, β], Ng)
        # tugrid = Grid.Uniform{Float64,33}(0.0, β, (true, true))
        # kugrid = Grid.Uniform{Float64,33}(0.0, maxK, (true, true))
        f(t) = t

        hist = zeros(tgrid.size)
        for i in 1:Nmc
            x = rand(rng) * β
            ind = Interp.locate(tgrid, x)
            vol = Interp.volumn(tgrid, ind)
            hist[ind] += f(x) / vol * Interp.volumn(tgrid)
        end

        for (ti, t) in enumerate(tgrid.grid)
            if ti != 1 && ti != Ng
                @test isapprox(hist[ti] / Nmc, f(t), rtol=3 / sqrt(Nmc / Ng))
            else
                # edge points has extra error because grid point is not at center of interval
                @test isapprox(hist[ti] / Nmc, f(t), rtol=3 / sqrt(Nmc / Ng), atol = 3/Ng)
            end
        end

    end

end
