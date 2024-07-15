

struct dynamic_problemSampled#{TCompFloat,Tint,Tfloat}
    Problem#::DDEProblem #e.g. =DDEProblem(....)
    alg#::Dict{Symbol, Any}               XXXX#::MethodOfSteps # e.g.: alg = MethodOfSteps(Tsit5())
    maxdelay#::Tfloat
    Tperiod#::Tfloat
    #dt#::Tfloat
    StateSmaplingTime#::LinRange#Vector{Tfloat}
    eigN#::Tint # number of eigen vectors
    #eigs#::Vector{TCompFloat}
    #eigsA#::Vector{Vector{TCompFloat}}
    zerofixpont::Bool
    affineinteration::Int#::Bool
    #adaptive::Bool#Variable stepsize
    KrylovTol#::Float64
    KrylovExtraDim::Int
    #fixpont
end

#DODOAU
#EPSI_TODO_REMOVE
