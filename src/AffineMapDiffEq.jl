module AffineMapDiffEq

using DifferentialEquations
using KrylovKit
using StaticArrays
using LinearAlgebra
using Interpolations
using ForwardDiff
using Random: rand!

# Export core types and functions
export dynamic_problemSampled, 
       affine, 
       spectralradius, 
       spectrum, 
       LinMap,
       partialpart, 
       valuepart, 
       find_fix_pont

# Include sub-modules/files
include("types.jl")
include("interpolation.jl")
include("mapping.jl")
include("spectrum.jl")

end # module AffineMapDiffEq
