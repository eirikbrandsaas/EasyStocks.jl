## Load everything
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("src/EasyStocks.jl")

## Setup
using CairoMakie # For plotting
## Run-code

@time MS = SolveModel()

fig = Figure()
lines(fig[1,1],MS.np.xgrd, MS.α[:,1])
ylims!(-0.05,1.05)
lines(fig[1,2],MS.np.xgrd, MS.V[:,1])
ylims!(-6,1)
lines(fig[2,1],MS.np.xgrd, MS.s[:,1])
lines(fig[2,2],MS.np.xgrd, MS.V[:,2])
lines(fig[3,1],MS.np.xgrd, MS.b[:,1])
lines(fig[3,2],MS.np.xgrd, MS.b[:,1] + MS.s[:,1])
ylims!(-0.05,5)
fig
