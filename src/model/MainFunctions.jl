function SolveModel()
  mp = ModPar()
  np = NumPar(mp)
  MS = ModelSolution(mp,np)

  LastPeriod!(MS)
  FirstPeriod!(MS)

  return MS

end

function LastPeriod!(MS::ModelSolution)
  np = MS.np
  mp = MS.mp
  ia = np.na
  for (ix,xv) in enumerate(np.xgrd)
    MS.c[ix,ia] = c = xv
    MS.V[ix,ia] = util(c,mp.γ)
  end
end

function FirstPeriod!(MS::ModelSolution)
  np = MS.np
  mp = MS.mp

  nchoice = np.nx*10
  sgrd = range(0.0,stop=np.xgrd[end],length=nchoice)
  bgrd = range(0.0,stop=np.xgrd[end],length=nchoice)
  vtmp = fill(0.0,(nchoice,nchoice))

  ia = 1
  Vnxt_intrp = LinearInterpolation((np.xgrd,),MS.V[:,ia+1],extrapolation_bc=Interpolations.Flat())
  for (ix,xv) in enumerate(np.xgrd)
    vtmp .= -Inf64
    for (is,sv) in enumerate(sgrd)
      for (ib,bv) in enumerate(bgrd)
        c = findc(xv,bv,sv,mp.q)
        if c > 0
          xpn = xp(bv,sv,mp)
          vtmp[is,ib] = util(c,mp.γ) + mp.β*Vnxt_intrp(xpn)
        else
          vtmp[is,ib] = -Inf64
        end
      end
    end

    imax = argmax(vtmp)
    is = imax[1]
    ib = imax[2]
    MS.s[ix,ia] = sgrd[is]
    MS.b[ix,ia] = bgrd[ib]
    MS.V[ix,ia] = vtmp[is,ib]

  end

end

function findc(x,b,s,q) # Find todays consumption
  c = x - b - s - q*part(s)
end

function xp(b,s,mp::ModPar) # Next period wealth
  b*(1.0+mp.r) + s*(1.0 + mp.μ)
end

function part(s::AbstractFloat) # Whether you are participating or not
  if s>0.0
      part = true
  else
      part = false
  end

  return part
end
