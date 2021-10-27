mutable struct ModPar
  γ :: Float64
  β :: Float64
  q :: Float64
  r :: Float64
  σ :: Float64
  μ :: Float64

  function ModPar(;
    γ = 2.0,
    β = 1.0,
    q = 0.0,
    r = 0.0,
    σ = 0.1,
    μ = 0.0,
  )

    new(γ, β, q, r, σ, μ)
  end
end

struct NumPar
  nx :: Int
  nh :: Int
  na :: Int

  xmax :: Float64
  hmax :: Float64

  xgrd :: Vector{Float64}
  hgrd :: Vector{Float64}

  function NumPar(;
      nx=11,
      nh=11,
      na=2,
      xmax=1.0,
      hmax=1.0,)


    xgrd = range(0,stop=xmax,length=nx)
    hgrd = range(0,stop=xmax,length=nx)

    new(nx, nh, na, xmax, hmax, xgrd, hgrd)
  end
end

mutable struct ModelSolution
  V :: Array{Float64}
  α :: Array{Float64}
  b :: Array{Float64}
  s :: Array{Float64}
  h :: Array{Float64}
  c :: Array{Float64}
  mp :: ModPar
  np :: NumPar

  function ModelSolution(mp::ModPar,np::NumPar)
    V = fill(-Inf64,(np.nx,np.na))
    α = fill(-Inf64,(np.nx,np.na))
    b = fill(-Inf64,(np.nx,np.na))
    s = fill(-Inf64,(np.nx,np.na))
    h = fill(-Inf64,(np.nx,np.na))
    c = fill(-Inf64,(np.nx,np.na))

    new(V, α, b, s, h, c, mp, np)
  end
end
