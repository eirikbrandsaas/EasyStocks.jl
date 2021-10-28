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
  for (ix,xv) in enumerate(np.xgrd)
    MS.s[ix,ia] = b = 0.0
    MS.b[ix,ia] = s = 0.0
    MS.c[ix,ia] = c = findc(xv,y,b,s,mp.q)
    MS.V[ix,ia] = util(c,mp.γ) + KKKShifter(mp.xstar,xv,mp.ψ)
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
  Vnxt_intrp = CubicSplineInterpolation((np.xgrd,),MS.V[:,ia+1],extrapolation_bc=Interpolations.Line())
  for (ix,xv) in enumerate(np.xgrd)
    vtmp .= -Inf64
    for (isav,sav) in enumerate(savgrd)
      for (iα,α) in enumerate(αgrd)
        b = sav*(1.0-α)
        s = sav*α
        c = findc(xv,y,b,s,mp.q)
        if c > 0
          xpn = xp(b,s,mp.r,np.rsgrd)
          vtmp[isav,iα] = util(c,mp.γ) + mp.β*sum(np.πrs .* Vnxt_intrp.(xpn))
        end
      end
    end

    imax = argmax(vtmp)
    isav = imax[1]
    iα = imax[2]
    MS.α[ix,ia] = α = αgrd[iα]
    MS.s[ix,ia] = s = savgrd[isav]*α
    MS.b[ix,ia] = b = savgrd[isav]*(1.0-α)
    if savgrd[isav] == 0.0 # If there is no saving, portfolio weight is'nt defined
      MS.α[ix,ia] = NaN64
    end
    MS.c[ix,ia] = findc(xv,y,s,b,mp.q)
    MS.V[ix,ia] = vtmp[imax]

  end

end

function findc(x,y,b,s,q) # Find todays consumption
  c = x + y - b - s - q*part(s)
end

function xp(b,s,r,rsgrd::Vector{Float64}) # Next period wealth
  b*(1.0+r) .+ s*(1.0 .+ rsgrd)
end

function part(s::AbstractFloat) # Whether you are participating or not
  if s>0.0
      part = true
  else
      part = false
  end

  return part
end

function KKKShifter(xstar::AbstractFloat, x::AbstractFloat,ψ::AbstractFloat)
  ψ*(x>xstar)
end


