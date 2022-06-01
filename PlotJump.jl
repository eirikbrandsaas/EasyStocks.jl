##
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("src/EasyStocks.jl")

## Setup
using CairoMakie # For plotting
using Colors # To make it easier to select colors

mp = ModPar(q=0.0,ψ=0.0,xstar=2.0,η=0.)
np = NumPar(mp,nx=400,nsav=120,nq=1001,nα=101,nh=1,ygrd=[0.0,0.0],xmax=4.0,p1cons=false)
np = NumPar(mp,nx=150,nsav=120,nq=501,nα=101,nh=1,ygrd=[0.0,0.0],xmax=5.0,p1cons=false)
np.hgrd.=0.0
MS1 = ModelSolution(mp,np)
@time SolveModel!(MS1)
ψ = 0.015
mp = ModPar(q=mp.q,ψ=ψ,xstar=mp.xstar,η =mp.η)
MS2 = ModelSolution(mp,np)
@time SolveModel!(MS2)
np.πrs[1] = 0.002
np.πrs[:] = np.πrs/sum(np.πrs)

mp.ψ=0.0
mp.q=0.01
MS1c = ModelSolution(mp,np)
@time SolveModel!(MS1c)
mp.ψ=ψ
MS2c = ModelSolution(mp,np)
@time SolveModel!(MS2c)
mp.q=0.0

mp.ψ=0.0
MS1d = ModelSolution(mp,np)
@time SolveModel!(MS1d)
mp.ψ=ψ
MS2d = ModelSolution(mp,np)
@time SolveModel!(MS2d)
##
# Plot settings
using Colors # To make it easier to select colors
blue = colorant"#2B6EB2"
fig = Figure(resolution = (700, 300))
ax1  = Axis(fig[1,2],xlabel = L"Wealth $x$", ylabel = L"\alpha")
lines!(ax1,MS1.np.xgrd, MS1.α[:,1,1], linestyle=:dash,color=blue,label="Benchmark")
ixstar = findfirst(x->x==1,(np.xgrd .-mp.xstar).>0)
lines!(ax1,MS2.np.xgrd[1:ixstar-1], MS2.α[1:ixstar-1,1,1], color=:orange,label="Housing")
lines!(ax1,MS2.np.xgrd[ixstar:end], MS2.α[ixstar:end,1,1], color=:orange)
ylims!(ax1,(-0.05,1.05))
xlims!(ax1,(1,3))
fig

ax2  = Axis(fig[1,1],xlabel = L"Wealth $x$", ylabel = L"Utility $u(x')$")
lines!(ax2,MS1.np.xgrd, MS1.V[:,1,2], linestyle=:dash,color=blue)
ixstar = findfirst(x->x==1,(np.xgrd .-mp.xstar).>0)
lines!(ax2,MS2.np.xgrd[1:ixstar-1], MS2.V[1:ixstar-1,1,2], color=:orange)
lines!(ax2,MS2.np.xgrd[ixstar:end], MS2.V[ixstar:end,1,2], color=:orange)
ylims!(ax2,(-1.035,-0.165))
xlims!(ax2,(1,3))

leg = Legend(fig[2, 1:2], ax1,orientation = :horizontal,framevisible = false)
Label(fig[1, 1, Top()], L"Terminal Value $U(x')$", valign = :bottom,padding=[0 0 5 0])
Label(fig[1, 2, Top()], L"Portfolio Weight $\alpha$", valign = :bottom,padding=[0 0 5 0])
rowsize!(fig.layout, 2, Fixed(0.1))
hidespines!(ax2,:t,:r)
hidespines!(ax1,:t,:r)
save("figs/value_weight.pdf", fig)
fig

## Plots with disaster
fig = Figure(resolution = (700, 300))
ax1  = Axis(fig[1,1],xlabel = L"Wealth $x$", ylabel = L"\alpha")
lines!(ax1,MS1d.np.xgrd, MS1d.α[:,1,1], linestyle=:dash,color=blue,label="Benchmark")
lines!(ax1,MS2d.np.xgrd[1:ixstar-1], MS2d.α[1:ixstar-1,1,1], color=:orange,label="Housing")
lines!(ax1,MS2d.np.xgrd[ixstar:end], MS2d.α[ixstar:end,1,1], color=:orange)
ylims!(ax1,(-0.05,1.05))
xlims!(ax1,(1,3))

ax2  = Axis(fig[1,2],xlabel = L"Wealth $x$", ylabel = L"\alpha")
discont=sort!([1;ixstar;findfirst(x->x>0,MS2c.α[:,1,1]);findlast(x->x==0,MS2c.α[:,1,1])+1;length(MS2c.α[:,1,1])])
lines!(ax2,MS1d.np.xgrd[1:discont[2]-1], MS1c.α[1:discont[2]-1,1,1], linestyle=:dash,color=blue,label="Benchmark")
lines!(ax2,MS1d.np.xgrd[discont[2]:end], MS1c.α[discont[2]:end,1,1], linestyle=:dash,color=blue)
for i in 1:length(discont)-1
    lines!(ax2,MS2d.np.xgrd[discont[i]:discont[i+1]-1], MS2c.α[discont[i]:discont[i+1]-1], color=:orange,label="Housing")
end
ylims!(ax2,(-0.05,1.05))
xlims!(ax2,(1,3.))
leg = Legend(fig[2, 1:2], ax1,orientation = :horizontal,framevisible = false)
Label(fig[1, 1, Top()], L"Disaster Risk $p_{tail}>0$", valign = :bottom,padding=[0 0 5 0])
Label(fig[1, 2, Top()], L"Participation cost $q>0$", valign = :bottom,padding=[0 0 5 0])
rowsize!(fig.layout, 2, Fixed(0.1))
hidespines!(ax2,:t,:r)
hidespines!(ax1,:t,:r)
save("figs/disaster_cost.pdf", fig)
fig

