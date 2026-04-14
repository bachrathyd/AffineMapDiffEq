

mutable struct dynamic_problemSampled#{TCompFloat,Tint,Tfloat}
    Problem#::DDEProblem #e.g. =DDEProblem(....)
    alg#::Dict{Symbol, Any}               XXXX#::MethodOfSteps # e.g.: alg = MethodOfSteps(Tsit5())
    maxdelay#::Tfloat
    Tperiod#::Tfloat # deprecated - the time span of the problem defined in the DDEProblem is used
    #dt#::Tfloat
    zerofixpont::Bool
    affineinteration::Int#::Bool
    
    
    StateSmaplingTime#::LinRange#Vector{Tfloat}
    
    Krylov_arg#::Dict{Symbol, Any}               XXXX#::MethodOfSteps # e.g.: alg = MethodOfSteps(Tsit5())
   
    #eigN#::Tint # number of eigen vectors
    #KrylovTol#::Float64
    #KrylovExtraDim::Int
    #Δu_Δλ # for continuation
end

#DODOAU
#EPSI_TODO_REMOVE
