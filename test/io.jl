@testset "IO" begin
	  @testset "JLD2" begin
        function deepequal(a, b)
            typeof(a) == typeof(b) || return false
            N = fieldcount(typeof(a))
            for i in 1:N
                deepequal(getfield(a, i) , getfield(b, i)) || return false
            end
            return true
        end

	      β = 10

        tgrid = CompositeGrid.LogDensedGrid(
        :gauss,# The top layer grid is :gauss, optimized for integration. For interpolation use :cheb
        [0.0, β],# The grid is defined on [0.0, β]
        [0.0, β],# and is densed at 0.0 and β, as given by 2nd and 3rd parameter.
        5,# N of log grid
        0.005, # niminum interval length of log grid
        5 # N of bottom layer
        )

        ############# FileIO API #################
        save("example.jld2", Dict("grid" => tgrid))

        d = load("example.jld2")
        dg = d["grid"]
        println(typeof(dg))
        println(dg.grid)
        @test deepequal(dg, tgrid)

        ############# naive API ##################
        jldopen("example.jld2", "w") do file
            file["test"] = tgrid
        end

        jldopen("example.jld2", "r") do file
            g = file["test"]
            println(typeof(g))
            println(g.grid)
            @test deepequal(g, tgrid)
        end

        # remove the file
        rm("example.jld2")
    end
end
