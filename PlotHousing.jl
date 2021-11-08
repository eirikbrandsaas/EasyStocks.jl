##
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("src/EasyStocks.jl")

## Setup
using CairoMakie # For plotting

## Solve
mp1 = ModPar(q=0.005,ψ=0.0,xstar=.0,η=0.5)
mp1.γ = 2 + 1/(1-mp1.η)
np1 = NumPar(mp1,nx=100,nsav=150,nq=5,nα=101,nh=1,xmax=10)
MS1 = ModelSolution(mp1,np1)
@time SolveModel!(MS1)

## With housing or jump?
# Jump:
mp2 = ModPar(q=mp1.q,ψ=0.015,xstar=mp1.xstar,η =mp1.η,γ = mp1.γ)
np2 = NumPar(mp1,nx=np1.nx,nsav=np1.nsav,nq=np1.nq,nα=np1.nα,nh=1,xmax=np1.xmax)
# Housing:
mp2 = ModPar(q=mp1.q,ψ=mp1.ψ,xstar=mp1.xstar,η =mp1.η,γ = mp1.γ)
np2 = NumPar(mp2,nx=np1.nx,nsav=np1.nsav,nq=np1.nq,nα=np1.nα,nh=2,xmax=np1.xmax)

np2.hgrd[1] = 1.
np2.hgrd[2] = 3.0
MS2 = ModelSolution(mp2,np2)
@time SolveModel!(MS2)

np2.πrs[1] += 0.005
np2.πrs[:] = np2.πrs/sum(np2.πrs)
MS3 = ModelSolution(mp2,np2)
@time SolveModel!(MS3)
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

ih = 1

lines!(s3,MS1.s[:,ih] + MS1.b[:,ih],linestyle=:dash,color=:gray,)
lines!(s3,MS2.s[:,ih] + MS2.b[:,ih],color=:orange,)

ia = 2
lines!(h2,MS1.np.xgrd,MS1.h[:,ih,ia],linestyle=:dash,color=:gray,)
lines!(h2,MS2.np.xgrd,MS2.h[:,ih,ia],color=:orange,)
lines!(h2,MS2.np.xgrd,MS3.h[:,ih,ia],color=:black,)
ylims!(h2,0.5,2.5)

ia = 2
ih = 1
lines!(v2,MS1.np.xgrd,MS1.V[:,1,ia],linestyle=:dash,color=:gray,)
lines!(v2,MS2.np.xgrd,MS2.V[:,1,ia],color=:orange,)
lines!(v2,MS3.np.xgrd,MS3.V[:,1,ia],color=:black,)

# rowsize!(g2,2,Relative(1.33))

ih = 1
lines!(ax1,MS1.np.xgrd, MS1.α[:,ih])
ylims!(ax1,-0.05,1.05)
# text!(ax1, "Not willing\nto pay cost \n ↓", position=(0.0,0.5),align=(:left,:center))

ih = 1
lines!(ax2,MS2.np.xgrd, MS1.α[:,ih], linestyle=:dash,color=:gray,)
lines!(ax2,MS2.np.xgrd, MS2.α[:,ih], color=:orange)
ylims!(ax2,-0.05,1.05)
# text!(ax2, "↑\n Gambling\n for housing", position=(3.2,0.2),align=(:right,:center))
# text!(ax2, "Stay safe \nfor housing\n↓", position=(3.8,0.2),align=(:left,:center))
# text!(ax2, "Gradually re-enter\n↓", position=(4.7,0.8),align=(:left,:center))

ih = 1
lines!(ax3,MS2.np.xgrd, MS2.α[:,ih],color=:orange)
lines!(ax3,MS2.np.xgrd, MS3.α[:,ih],color=:black)
ylims!(ax3,-0.05,1.05)

fig
## Run-code

fig = Figure()
ax1 = Axis(fig[1,1])
MS = MS1
hitp = LinearInterpolation(MS.np.xgrd,MS2.h[:,2],extrapolation_bc=Flat())

probh2 = fill(-1.0,MS.np.nx)
for ix in eachindex(MS.np.xgrd)
  xns = MS.s[ix,1]*(1.0.+MS.np.rsgrd) .+ MS.b[ix,1] *(1.0 + MS.mp.r)
  probh2[ix] = sum(MS.np.πrs.*(hitp.(xns).==MS2.np.hgrd[2]))
end

scatter!(ax1,probh2)

MS = MS2
hitp = LinearInterpolation(MS2.np.xgrd,MS.h[:,2],extrapolation_bc=Flat())

probh2 = fill(-1.0,MS.np.nx)
for ix in eachindex(MS.np.xgrd)
  xns = MS.s[ix,1]*(1.0.+MS.np.rsgrd) .+ MS.b[ix,1] *(1.0 + MS.mp.r)
  probh2[ix] = sum(MS.np.πrs.*(hitp.(xns).==MS2.np.hgrd[2]))
end

scatter!(ax1,probh2.+1.0)
fig
