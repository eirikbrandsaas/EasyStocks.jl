function util(c::AbstractFloat,γ::AbstractFloat)
  (c^(1.0-γ))/(1.0-γ)
end

function utilh(c::AbstractFloat,h::AbstractFloat,γ::AbstractFloat,η::AbstractFloat)
  ((c^(1.0-η)*h^η)^(1.0-γ))/(1.0-γ)
end
