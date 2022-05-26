##
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("src/EasyStocks.jl")

## Setup
using CairoMakie # For plotting

mp = ModPar(q=0.0,ψ=0.0,xstar=2.0,η=0.)
np = NumPar(mp,nx=60,nsav=120,nq=15,nα=81,nh=1,ygrd=[1.0,1.0])
MS1 = ModelSolution(mp,np)
@time SolveModel!(MS1)
mp = ModPar(q=mp.q,ψ=0.015,xstar=mp.xstar,η =mp.η)
MS2 = ModelSolution(mp,np)
@time SolveModel!(MS2)

# np.πrs[1] = 0.002
# np.πrs[:] = np.πrs/sum(np.πrs)
# MS3 = ModelSolution(mp,np)
# @time SolveModel!(MS3)
##
fig = Figure()
g2 = fig[1,1] = GridLayout()
g1 = fig[2,1] = GridLayout()

ax2  = Axis(g1[1,1],xlabel = "Wealth", ylabel = "α")
ax3  = Axis(g1[2,1],xlabel = "Wealth", ylabel = "s")
v2 = Axis(g2[1,1], xlabel ="Wealth", ylabel = "Utils")
v1 = Axis(g2[1,2], xlabel ="Wealth", ylabel = "Utils")
lines!(v1,np.xgrd,MS1.V[:,1,1],linestyle=:dash,color=:gray,label="Benchmark")
lines!(v1,np.xgrd,MS2.V[:,1,1],color=:orange,label="Housing")
ylims!(v1,(MS2.V[np.nx÷7,1,1],0))


lines!(v2,np.xgrd,MS1.V[:,1,2],linestyle=:dash,color=:gray,)
lines!(v2,np.xgrd,MS2.V[:,1,2],color=:orange,)

lines!(ax3,MS2.np.xgrd, MS1.s[:,1,1]+MS1.b[:,1,1], linestyle=:dash,color=:gray)
lines!(ax3,MS2.np.xgrd, MS2.s[:,1,1]+MS2.b[:,1,1], color=:orange)
lines!(ax3,MS2.np.xgrd, MS1.s[:,1,1], linestyle=:dash,color=:blue)
lines!(ax3,MS2.np.xgrd, MS2.s[:,1,1],color=:blue)
lines!(ax3,MS2.np.xgrd, MS1.b[:,1,1], linestyle=:dash,color=:green)
lines!(ax3,MS2.np.xgrd, MS2.b[:,1,1],color=:green)

lines!(ax2,MS2.np.xgrd, MS1.α[:,1,1], linestyle=:dash,color=:gray)
lines!(ax2,MS2.np.xgrd, MS2.α[:,1,1], color=:orange)
ylims!(ax2,(-0.05,1.05))
text!(ax2, "↑\n Gambling\n for housing", position=(3.4,0.3),align=(:right,:center))
text!(ax2, "Stay safe \nfor housing\n↓", position=(4.,0.3),align=(:left,:center))
text!(ax2, "Gradually \nre-enter\n↓", position=(4.9,0.8),align=(:left,:center))

# lines!(ax3,MS2.np.xgrd, MS2.α[:,1],)
# lines!(ax3,MS2.np.xgrd, MS3.α[:,1],)
# ylims!(ax3,-0.05,1.05)
#leg = Legend(g1[1, 2], ax2)    
axislegend(v1,position=:rb)
Label(g2[1, 1, Top()], "Second Period Value Function", valign = :bottom,)
Label(g2[1, 2, Top()], "First Period Value Function", valign = :bottom,)

fig
## Run-code
function Prob_xstar(MS)
    tmp = fill(0.0,MS.np.nx)
    for ix in eachindex(MS.np.xgrd)
        b = MS.b[ix,1,1]
        s = MS.s[ix,1,1]
        for irs in eachindex(np.πrs)
            tmp[ix] += np.πrs[irs]*(xp(b,s,mp.r,np.rsgrd[irs]) > MS.mp.xstar)
        end
    end
    return tmp
end

Prob_xstar(MS3)
scatter!(ax2,MS2.np.xgrd, Prob_xstar(MS1), linestyle=:dash,color=:gray)
scatter!(ax2,MS2.np.xgrd, Prob_xstar(MS2), color=:orange)
fig