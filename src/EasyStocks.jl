
using Interpolations # To interpolate value function
using Expectations # To easily find expected values
using Distributions # To create normal stock returns =)
using CairoMakie # For plotting
# Load files

include("model/Structs.jl")
include("model/Fundamentals.jl")
include("model/MainFunctions.jl")
