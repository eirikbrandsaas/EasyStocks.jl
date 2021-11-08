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

@testset "Sales Costs" begin
    @test SalesCost(1,2,0.0)==0.0 # Test that cost = 0  when sizes are different, but cost is 0
    @test SalesCost(1,1,0.1)==0.0 # Test that cost = 0  when sizes are the same
    @test SalesCost(1,2,0.1)>0.0 # Test that cost > 0  when sizes are not the same
end
