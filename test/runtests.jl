using EasyStocks
using Test

@testset "EasyStocks.jl" begin
    # Write your tests here.
    @test minimum(SolveModel())<1.0
end
