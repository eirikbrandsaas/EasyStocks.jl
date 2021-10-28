## Load everything
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("src/EasyStocks.jl")
include("test/runtests.jl")

## Setup
using CairoMakie # For plotting
mp = ModPar(q=0.005,ψ=0.015,xstar=2.0)
np = NumPar(mp)
MS = ModelSolution(mp,np)
## Run-code

@time SolveModel!(MS)

fig = Figure()
lines(fig[1,1],MS.np.xgrd, MS.α[:,1])
ylims!(-0.05,1.05)
lines(fig[1,2],MS.np.xgrd, MS.V[:,1])
ylims!(-6,1)
lines(fig[2,2],MS.np.xgrd, MS.V[:,2])
ylims!(-1.5,0)
lines(fig[2,1],MS.np.xgrd, MS.b[:,1] + MS.s[:,1])
ylims!(-0.05,5)
linkxaxes!(fig.content...)
fig
