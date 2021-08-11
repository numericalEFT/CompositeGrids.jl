@testset "Interpolate" begin
    @testset "Linear1D" begin
        β = π
        tgrid = SimpleGrid.Uniform{Float64}([0.0, β], 33)
        # tugrid = Grid.Uniform{Float64,33}(0.0, β, (true, true))
        # kugrid = Grid.Uniform{Float64,33}(0.0, maxK, (true, true))
        f(t) = t
        data = zeros(tgrid.size)

        for (ti, t) in enumerate(tgrid.grid)
            data[ti] = f(t)
        end

        for ti = 1:tgrid.size - 1
            t = tgrid[ti] + 1.e-6
            fbar = Interp.interp1D(data, tgrid, t)
            @test abs(f(tgrid[ti]) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt
            @test f(tgrid[ti]) < fbar
            @test f(tgrid[ti + 1]) > fbar
        end
        for ti = 2:tgrid.size
            t = tgrid[ti] - 1.e-6
            fbar = Interp.interp1D(data, tgrid, t)
            @test abs(f(tgrid[ti]) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt
            @test f(tgrid[ti]) > fbar
            @test f(tgrid[ti - 1]) < fbar
        end

        t = tgrid[1] + eps(Float64)*1e3
        fbar = Interp.interp1D(data, tgrid, t)
        @test abs(f(t) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt

        t = tgrid[tgrid.size] - eps(Float64)*1e3
        fbar = Interp.interp1D(data, tgrid, t)
        @test abs(f(t) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt

        tlist = rand(10) * β
        # println(tlist)

        for (ti, t) in enumerate(tlist)
            fbar = Interp.interp1D(data, tgrid, t)
            # println("$k, $t, $fbar, ", f(k, t))
            @test abs(f(t) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt
        end
    end

    @testset "BaryCheb" begin
        β = π
        tgrid = SimpleGrid.BaryCheb{Float64}([0.0, β], 16)
        # tugrid = Grid.Uniform{Float64,33}(0.0, β, (true, true))
        # kugrid = Grid.Uniform{Float64,33}(0.0, maxK, (true, true))
        f(t) = t
        data = zeros(tgrid.size)

        for (ti, t) in enumerate(tgrid.grid)
            data[ti] = f(t)
        end

        for ti = 1:tgrid.size - 1
            t = tgrid[ti] + 1.e-6
            fbar = Interp.interp1D(data, tgrid, t)
            @test abs(f(tgrid[ti]) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt
            @test f(tgrid[ti]) < fbar
            @test f(tgrid[ti + 1]) > fbar
        end
        for ti = 2:tgrid.size
            t = tgrid[ti] - 1.e-6
            fbar = Interp.interp1D(data, tgrid, t)
            @test abs(f(tgrid[ti]) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt
            @test f(tgrid[ti]) > fbar
            @test f(tgrid[ti - 1]) < fbar
        end

        t = tgrid[1] + eps(Float64)*1e3
        fbar = Interp.interp1D(data, tgrid, t)
        @test abs(f(t) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt

        t = tgrid[tgrid.size] - eps(Float64)*1e3
        fbar = Interp.interp1D(data, tgrid, t)
        @test abs(f(t) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt

        tlist = rand(10) * β
        # println(tlist)

        for (ti, t) in enumerate(tlist)
            fbar = Interp.interp1D(data, tgrid, t)
            # println("$k, $t, $fbar, ", f(k, t))
            @test abs(f(t) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt
        end
    end

    @testset "DensedLog" begin
        β = 4
        tgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, β], [0.0, 0.5β, β], 4, 0.001, 4)
        # tugrid = Grid.Uniform{Float64,33}(0.0, β, (true, true))
        # kugrid = Grid.Uniform{Float64,33}(0.0, maxK, (true, true))
        f(t) = t
        data = zeros(tgrid.size)

        for (ti, t) in enumerate(tgrid.grid)
            data[ti] = f(t)
        end

        for ti = 1:tgrid.size - 1
            t = tgrid[ti] + 1.e-6
            fbar = Interp.interp1D(data, tgrid, t)
            @test abs(f(tgrid[ti]) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt
            @test f(tgrid[ti]) < fbar
            @test f(tgrid[ti + 1]) > fbar
        end
        for ti = 2:tgrid.size
            t = tgrid[ti] - 1.e-6
            fbar = Interp.interp1D(data, tgrid, t)
            @test abs(f(tgrid[ti]) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt
            @test f(tgrid[ti]) > fbar
            @test f(tgrid[ti - 1]) < fbar
        end

        t = tgrid[1] + eps(Float64)*1e3
        fbar = Interp.interp1D(data, tgrid, t)
        @test abs(f(t) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt

        t = tgrid[tgrid.size] - eps(Float64)*1e3
        fbar = Interp.interp1D(data, tgrid, t)
        @test abs(f(t) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt

        tlist = rand(10) * β
        tlist = sort(tlist)
        println(tlist)
        ff = Interp.interpGrid(data, tgrid, tlist)
        

        for (ti, t) in enumerate(tlist)
            fbar = Interp.interp1D(data, tgrid, t)
            # println("$k, $t, $fbar, ", f(k, t))
            @test abs(f(t) - fbar) < 3.e-6 # linear interpolation, so error is δK+δt
            @test abs(f(t) - ff[ti]) < 3.e-6 # linear interpolation, so error is δK+δt
        end
    end

    @testset "Integrate" begin
        β = 1.0
        tgrid = SimpleGrid.GaussLegendre{Float64}([0.0, β], 4)
        # tugrid = Grid.Uniform{Float64,33}(0.0, β, (true, true))
        # kugrid = Grid.Uniform{Float64,33}(0.0, maxK, (true, true))
        f(t) = t
        data = zeros(tgrid.size)
        for (ti, t) in enumerate(tgrid.grid)
            data[ti] = f(t)
        end
        println(tgrid.grid)
        println(data)
        println(tgrid.weight)
        println(sum(data.*tgrid.weight))
        int_result = Interp.integrate1D(data, tgrid)
        @test abs(int_result - 0.5) < 3.e-6

        β = 1.0
        tgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, β], [0.0, 0.5β, β], 2, 0.001, 3)
        # tugrid = Grid.Uniform{Float64,33}(0.0, β, (true, true))
        # kugrid = Grid.Uniform{Float64,33}(0.0, maxK, (true, true))
        data = zeros(tgrid.size)
        for (ti, t) in enumerate(tgrid.grid)
            data[ti] = f(t)
        end
        println(tgrid.grid)
        println(data)
        int_result = Interp.integrate1D(data, tgrid)
        @test abs(int_result - 0.5) < 3.e-6
    end
end

