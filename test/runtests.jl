include("../src/EasyStocks.jl")
using Test

@testset "EasyStocks.jl" begin
    # Write your tests here.
    @test 1+1==2
end

@testset "Make sure all functions run" begin
    @test try
        mp = ModPar()
        np = NumPar(mp)
        MS = ModelSolution(mp,np)
        SolveModel!(MS)
        true
    catch
        false
    end
end
