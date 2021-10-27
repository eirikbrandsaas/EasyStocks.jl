module EasyStocks

# Dependencies
using Base: AbstractFloat
using Interpolations

# Load files

include("model/Structs.jl")
include("model/Fundamentals.jl")
include("model/MainFunctions.jl")

# Export functions
export SolveModel,
       ModPar,
       NumPar,
       ModelSolution
end

