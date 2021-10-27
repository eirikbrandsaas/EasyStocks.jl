module EasyStocks

# Dependencies
using Interpolations

# Load files

include("model/MainFunctions.jl")
include("model/Structs.jl")

# Export functions
export SolveModel,
       ModPar,
       NumPar,
       ModelSolution
end

