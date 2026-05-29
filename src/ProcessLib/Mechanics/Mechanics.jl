# This file defines the main mechanics process, including the necessary functions.
# Note that the sequence of function definitions must be maintained to ensure that all necessary functions are defined before they are called.
include("ExternalForce.jl")
include("VonMisesReturnMapping.jl")
include("LocalAssemble.jl")
include("Solvers.jl")