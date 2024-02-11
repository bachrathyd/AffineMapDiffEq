using Profile
5+5 

using KrylovKit
using DifferentialEquations
using StaticArrays
using LinearAlgebra
using Interpolations



using ForwardDiff
#using Plots
#PythonPlot()
#plotlyjs()

#using  KrylovKit
#using  DifferentialEquations
#using  StaticArrays
#using  LinearAlgebra
#using  Interpolations

struct dynamic_problemSampled#{TCompFloat,Tint,Tfloat}
    DDEdynProblem#::DDEProblem #e.g. =DDEProblem(....)
    alg#::MethodOfSteps # e.g.: alg = MethodOfSteps(Tsit5())
    maxdelay#::Tfloat
    Tperiod#::Tfloat
    dt#::Tfloat
    StateSmaplingTime#::LinRange#Vector{Tfloat}
    eigN#::Tint # number of eigen vectors
    eigs#::Vector{TCompFloat}
    #eigsA#::Vector{Vector{TCompFloat}}
    zerofixpont#::Bool
    affineinteration::Int#::Bool
    #fixpont
end

#DODOAU
#EPSI_TODO_REMOVE
