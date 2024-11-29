module DDE_mapping


using KrylovKit
using DifferentialEquations
using StaticArrays
using LinearAlgebra
using Interpolations

using ForwardDiff



export dynamic_problemSampled,affine,spectralradius,
spectrum,LinMap,
partialpart,valuepart,find_fix_pont

include("DDE_mapping_types.jl")
include("DDE_mapping_functions.jl")

test()=println("test")

end # module DDE_mapping
