# AffineMapDiffEq.jl Interpolation Wrappers

using Interpolations
using StaticArrays

"""
    AffineHistory{I, V}

A type-stable wrapper for interpolation objects used in DDE history functions.
"""
struct AffineHistory{I, V}
    itp::I
    is_vector_of_vectors::Bool
    N_dim::Int
end

# Multi-dispatch call to avoid branching inside the ODE loop
@inline function (ah::AffineHistory{I, V})(p, t) where {I, V}
    if ah.is_vector_of_vectors
        return [ah.itp(i, t) for i in 1:ah.N_dim]
    else
        # Force the return type to V if it's a StaticArray
        return ah.itp(t)::V
    end
end

# Support for derivative hints if needed by some solvers
@inline function (ah::AffineHistory{I, V})(p, t, ::Type{Val{1}}) where {I, V}
    if ah.is_vector_of_vectors
        return [Interpolations.gradient(ah.itp, i, t)[2] for i in 1:ah.N_dim]
    else
        return Interpolations.gradient(ah.itp, t)[1]
    end
end

"""
    construct_history(s, StateSmaplingTime)

Analyzes the input state sequence and constructs an appropriate `AffineHistory` wrapper.
"""
function construct_history(s::T, StateSmaplingTime) where {T}
    if eltype(s) <: AbstractVector && !(eltype(s) <: StaticArray)
        s_mat = reduce(hcat, s)
        itp = interpolate(s_mat, (NoInterp(), BSpline(Cubic(Line(OnGrid())))))
        sitp = scale(itp, 1:size(s_mat, 1), StateSmaplingTime)
        return AffineHistory{typeof(sitp), eltype(s)}(sitp, true, size(s_mat, 1))
    else
        itp = interpolate(s, BSpline(Cubic(Line(OnGrid()))))
        sitp = scale(itp, StateSmaplingTime)
        return AffineHistory{typeof(sitp), eltype(s)}(sitp, false, 0)
    end
end
