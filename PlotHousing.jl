##
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("src/EasyStocks.jl")

## Setup
using CairoMakie # For plotting

mp1 = ModPar(q=0.00,ψ=0.015,xstar=.0,η=0.5)
mp1.γ = 2 + 1/(1-mp1.η)
np1 = NumPar(mp1,nx=100,nsav=150,nq=31,nα=101,nh=1,xmax=10)
MS1 = ModelSolution(mp1,np1)
@time SolveModel!(MS1)

## With housing or jump?
# Jump:
mp2 = ModPar(q=mp1.q,ψ=0.015,xstar=mp1.xstar,η =mp1.η,γ = mp1.γ)
np2 = NumPar(mp1,nx=np1.nx,nsav=np1.nsav,nq=np1.nq,nα=np1.nα,nh=1,xmax=np1.xmax)
# Housing:
mp2 = ModPar(q=mp1.q,ψ=mp1.ψ,xstar=mp1.xstar,η =mp1.η,γ = mp1.γ)
np2 = NumPar(mp2,nx=np1.nx,nsav=np1.nsav,nq=np1.nq,nα=np1.nα,nh=2,xmax=np1.xmax)

np2.hgrd[1] = 0.5
np2.hgrd[2] = 3.0
MS2 = ModelSolution(mp2,np2)
@time SolveModel!(MS2)

# np.πrs[1] = 0.00
# np.πrs[:] = np.πrs/sum(np.πrs)
# sum(np.πrs.*np.rsgrd)
# sum(np.πrs.*np.rsgrd.^2 .- sum(np.πrs.*np.rsgrd)^2)
# MS3 = ModelSolution(mp,np)
# @time SolveModel!(MS3)
##
fig = Figure()
g2 = fig[1,1] = GridLayout()
g1 = fig[2,1] = GridLayout()
ax1  = Axis(g1[1,1],xlabel = "", ylabel = "")
ax2  = Axis(g1[1,2],xlabel = "", ylabel = "")
ax3  = Axis(g1[1,3],xlabel = "", ylabel = "")

v2 = Axis(g2[1,1], xlabel ="", ylabel = "Utils")
h2 = Axis(g2[1,2], xlabel ="", ylabel = "Housing")
s3 = Axis(g2[1,3])

lines!(s3,MS1.s[:,1] + MS1.b[:,1],linestyle=:dash,color=:gray,)
lines!(s3,MS2.s[:,1] + MS2.b[:,1],color=:orange,)

lines!(h2,MS1.np.xgrd,MS1.h[:,2],linestyle=:dash,color=:gray,)
lines!(h2,MS2.np.xgrd,MS2.h[:,2],color=:orange,)
ylims!(h2,0.5,2.5)

lines!(v2,MS1.np.xgrd,MS1.V[:,2],linestyle=:dash,color=:gray,)
lines!(v2,MS2.np.xgrd,MS2.V[:,2],color=:orange,)

# rowsize!(g2,2,Relative(1.33))

lines!(ax1,MS1.np.xgrd, MS1.α[:,1])
ylims!(ax1,-0.05,1.05)
# text!(ax1, "Not willing\nto pay cost \n ↓", position=(0.0,0.5),align=(:left,:center))

lines!(ax2,MS2.np.xgrd, MS1.α[:,1], linestyle=:dash,color=:gray,)
lines!(ax2,MS2.np.xgrd, MS2.α[:,1], color=:orange)
ylims!(ax2,-0.05,1.05)
text!(ax2, "↑\n Gambling\n for housing", position=(3.2,0.2),align=(:right,:center))
text!(ax2, "Stay safe \nfor housing\n↓", position=(3.8,0.2),align=(:left,:center))
text!(ax2, "Gradually re-enter\n↓", position=(4.7,0.8),align=(:left,:center))

lines!(ax3,MS2.np.xgrd, MS2.α[:,1],)
# lines!(ax3,MS2.np.xgrd, MS3.α[:,1],)
ylims!(ax3,-0.05,1.05)

fig
## Run-code
