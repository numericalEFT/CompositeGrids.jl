using CompositeGrids, Test, StaticArrays, LinearAlgebra, Printf, Random, Statistics# , JLD2, FileIO
# import Test: @test, @testset

if isempty(ARGS)
    include("BaryChebTools.jl")
    include("SimpleG.jl")
    include("CompositeG.jl")
    include("Interp.jl")
    # include("io.jl")
    include("mc.jl")
    include("periodic.jl")
else
    include(ARGS[1])
end
