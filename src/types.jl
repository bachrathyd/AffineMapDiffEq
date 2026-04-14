# AffineMapDiffEq.jl Types

"""
    dynamic_problemSampled{P,A,T,ST,K,IT}

A structure containing all parameters and configurations for the sampled dynamic problem,
including the underlying DDE problem, solver arguments, and mapping settings.
"""
mutable struct dynamic_problemSampled{P,A,T,ST,K,IT}
    Problem::P             # DDEProblem (e.g., DDEProblem(....))
    alg::A                 # Solver options (e.g., Dict(:alg => MethodOfSteps(BS3())))
    maxdelay::T            # Maximum delay time
    zerofixpont::Bool      # Whether to assume the fixed point is at zero
    affineinteration::Int  # Number of affine mapping iterations
    StateSmaplingTime::ST  # Sampling points for the discretized state history
    Krylov_arg::K          # Arguments for the Krylov solver (e.g., KrylovKit.Arnoldi)
    perturbation_size::Float64 # Finite difference step size (0.0 for ForwardDiff)
    itp_type::IT           # Interpolation type (e.g., BSpline(Cubic(Line(OnGrid()))))
end
