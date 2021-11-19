function SolveModel!(MS::ModelSolution)

  LastPeriod!(MS)
  FirstPeriod!(MS)

  return MS

end

function LastPeriod!(MS::ModelSolution)
  np = MS.np
  mp = MS.mp
  ia = np.na
  y = np.ygrd[ia]
  vtmp = fill(-Inf64,np.nh)
  ctmp = fill(NaN64,np.nh)
  for (ix,xv) in enumerate(np.xgrd)
    for ih in eachindex(np.hgrd)
    b = 0.0
    s = 0.0
    for (ihn,hv) in enumerate(np.hgrd) # Find choices
      ms = SalesCost(ih,ihn,mp.ms)
      xtmp = xv-ms
      ctmp[ihn] = c = findc2(xtmp,y,b,s,mp.q,hv,np.hgrd[1])
      if c > 0
        vtmp[ihn] = utilh(c,hv,mp.γ,mp.η) + KKKShifter(mp.xstar,xv,mp.ψ)
      end
    end
    imax = argmax(vtmp)
    MS.h[ix,ih,ia] = np.hgrd[imax]
    MS.c[ix,ih,ia] = ctmp[imax]
    MS.V[ix,ih,ia] = vtmp[imax]
    end
  end
end

function FirstPeriod!(MS::ModelSolution)
  np = MS.np
  mp = MS.mp

  savgrd = range(0.0,stop=np.xgrd[end],length=np.nsav)
  αgrd = range(0.0,stop=1.0,length=np.nα)
  vtmp = fill(0.0,(np.nsav,np.nα))

  ia = 1
  y = np.ygrd[ia]
  for ih in eachindex(np.hgrd)
    Vnxt_intrp = CubicSplineInterpolation((np.xgrd,),MS.V[:,ih,ia+1],extrapolation_bc=Interpolations.Line())
  for (ix,xv) in enumerate(np.xgrd)
    vtmp .= -Inf64
    for (isav,sav) in enumerate(savgrd)
      for (iα,α) in enumerate(αgrd)
        b = sav*(1.0-α)
        s = sav*α
        c = findc1(xv,y,b,s,mp.q)
        if c > 0
          xpn = xp(b,s,mp.r,np.rsgrd)
          vtmp[isav,iα] = util(c,mp.γ) + mp.β*sum(np.πrs .* Vnxt_intrp.(xpn))
        end
      end
    end

    imax = argmax(vtmp)
    isav = imax[1]
    iα = imax[2]
    MS.α[ix,ih] = α = αgrd[iα]
    MS.s[ix,ih] = s = savgrd[isav]*α
    MS.b[ix,ih] = b = savgrd[isav]*(1.0-α)
    if savgrd[isav] == 0.0 # If there is no saving, portfolio weight is'nt defined
      MS.α[ix,ih] = NaN64
    end
    MS.c[ix,ih,ia] = findc1(xv,y,s,b,mp.q)
    MS.V[ix,ih,ia] = vtmp[imax]
  end
  end

end

" Find consumption as a function of choices"
function findc2(x,y,b,s,q,h,hmin) # Find todays consumption
  c = x + y - b - s - q*part(s) - (h-hmin) # First unit of housing is free :)
end

function findc1(x,y,b,s,q) # Find todays consumption
  c = x + y - b - s - q*part(s)
end

"Takes portfolio choices and returns possible next period wealth"
function xp(b,s,r,rsgrd::Vector{Float64})
  b*(1.0+r) .+ s*(1.0 .+ rsgrd)
end

"Returns true if household participates"
function part(s::AbstractFloat) # Whether you are participating or not
  if s>0.0
      part = true
  else
      part = false
  end

  return part
end

"Utility shifter as in Karlman, Kinnerud, Kragh-Soerensen (2021)"
function KKKShifter(xstar::AbstractFloat, x::AbstractFloat,ψ::AbstractFloat)
  ψ*(x>xstar)
end


function SalesCost(ih::Integer,ihn::Integer,κ::AbstractFloat)
  if ihn != ih
    cost = κ
  else
    cost = 0.0
  end

  return cost
end

