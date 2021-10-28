mutable struct ModPar
  γ :: Float64
  β :: Float64
  q :: Float64
  r :: Float64
  σ :: Float64
  μ :: Float64
  ψ :: Float64


  function ModPar(;
    γ = 2.0,
    β = 1.0,
    q = 0.0,
    r = 0.0,
    σ = 0.16,
    μ = 0.02,
    ψ = 0.0,
  )

    new(γ, β, q, r, σ, μ, ψ)
  end
end

struct NumPar
  nx :: Int
  nh :: Int
  na :: Int
  nq :: Int # How many points to use in quadrature
  nα :: Int # α-choice grid density
  nsav :: Int # savings-choice grid density

  xmax :: Float64
  hmax :: Float64

  xgrd :: StepRangeLen{Float64}
  hgrd :: Vector{Float64}
  ygrd :: Vector{Float64}

  πrs :: Vector{Float64}
  rsgrd :: Vector{Float64}

  function NumPar(mp::ModPar;
      nx=51,
      nh=11,
      na=2,
      nq=11,
      nα=21,
      nsav=51,
      xmax=8.0,
      hmax=1.0,)


    xgrd = range(0,stop=xmax,length=nx)
    hgrd = range(0,stop=xmax,length=nx)
    ygrd = [1.,1.]

    # Use expectation package
    E = expectation(Normal(mp.μ,mp.σ),n=nq)
    rsgrd=nodes(E) # Return grid
    rsgrd=max.(-1.0,rsgrd) # Limited liability
    πrs=weights(E) # PMF
    @assert sum(πrs)≈1 # Check that PMF sums to one

    new(nx, nh, na, nq, nα, nsav, xmax, hmax, xgrd, hgrd, ygrd, πrs, rsgrd)
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
