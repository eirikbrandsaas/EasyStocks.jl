using EasyStocks
using Test

@testset "EasyStocks.jl" begin
    # Write your tests here.
    @test 1+1==2
end

@testset "Make sure all functions run" begin
    @test try
        SolveModel()
        mp = ModPar()
        np = NumPar()
        ModelSolution(mp,np)
        true
    catch
        false
    end
end
