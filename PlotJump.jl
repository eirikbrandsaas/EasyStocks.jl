##
Pkg.activate(".")
Pkg.instantiate()
include("src/EasyStocks.jl")

## Setup
using CairoMakie # For plotting

mp = ModPar(q=0.0,ψ=0.0,xstar=2.0,η=0.)
np = NumPar(mp,nx=100,nsav=300,nq=41,nα=71,nh=1)
MS1 = ModelSolution(mp,np)
@time SolveModel!(MS1)
mp = ModPar(q=0.00,ψ=0.015,xstar=2.0,η =0.0)
MS2 = ModelSolution(mp,np)
@time SolveModel!(MS2)

np.πrs[1] = 0.002
np.πrs[:] = np.πrs/sum(np.πrs)
MS3 = ModelSolution(mp,np)
# @time SolveModel!(MS3)
##
fig = Figure()
g2 = fig[1,1] = GridLayout()
g1 = fig[2,1] = GridLayout()
ax1  = Axis(g1[1,1],xlabel = "", ylabel = "α")
ax2  = Axis(g1[2,1],xlabel = "", ylabel = "α")

v2 = Axis(g2[1,1], xlabel ="", ylabel = "Utils")
h2 = Axis(g2[1,2], xlabel ="", ylabel = "Housing")

lines!(h2,np.xgrd,MS1.h[:,2],linestyle=:dash,color=:gray,)
lines!(h2,np.xgrd,MS2.h[:,2],color=:orange,)
ylims!(h2,0.5,2.5)

lines!(v2,np.xgrd,MS1.V[:,2],linestyle=:dash,color=:gray,)
lines!(v2,np.xgrd,MS2.V[:,2],color=:orange,)
# ax3  = Axis(g1[3,1],xlabel = "", ylabel = "α")

# rowsize!(g2,2,Relative(1.33))
# axt  = Axis(fig[4,1],xlabel = "Wealth (w2)", ylabel = "Utils")
# lines!(axt,MS1.np.xgrd, MS1.V[:,2])
# Label(fig[4, 1, Top()], "Portfolio weight on stocks without housing", valign = :bottom)

lines!(ax1,MS1.np.xgrd, MS1.α[:,1])
ylims!(ax1,-0.05,1.05)
# text!(ax1, "Not willing\nto pay cost \n ↓", position=(0.0,0.5),align=(:left,:center))

lines!(ax2,MS2.np.xgrd, MS1.α[:,1], linestyle=:dash,color=:gray,)
lines!(ax2,MS2.np.xgrd, MS2.α[:,1], color=:orange)
ylims!(ax2,-0.05,1.05)
text!(ax2, "↑\n Gambling\n for housing", position=(3.2,0.2),align=(:right,:center))
text!(ax2, "Stay safe \nfor housing\n↓", position=(3.8,0.2),align=(:left,:center))
text!(ax2, "Gradually re-enter\n↓", position=(4.7,0.8),align=(:left,:center))

# lines!(ax3,MS2.np.xgrd, MS2.α[:,1],)
# lines!(ax3,MS2.np.xgrd, MS3.α[:,1],)
# ylims!(ax3,-0.05,1.05)

fig
## Run-code
