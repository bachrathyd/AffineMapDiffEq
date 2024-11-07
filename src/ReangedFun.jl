#ReangedFun

using LinearAlgebra
using BenchmarkTools
using StaticArrays


using FunctionWrappers
import FunctionWrappers: FunctionWrapper

#h(p, t::Float64) = ...
#h(p, t::Float64, deriv::Type{Val{1}}) = ...

#TODO: what is better? a "mutable struct" or s "struct" with f::Marry{1,FunctionWrapper{Tout,Tuple{Tin}} )
mutable struct RangeFun{Tin,Tout}#<: Function   F<:Function
    f::FunctionWrapper{Tout,Tuple{Tin}}
    range::MVector{2,Tin}
end

function Base.show(io::IO, fc::RangeFun) # I have to remove <: Function to have pretty-print
    println(io, "functinos: $(fc.f)")
    println(io, "range $(fc.range)")
end


# ------------- Initializations --------------------
#error("provide no type")
#function funcomp()
#    funcomp{Float64}(Float64[], funcomp[], [0.0, 0.0])
#end
function Base.return_types(fc::RangeFun{T,Tout}) where {T,Tout}
    return [Tout]
end

#TODO: output can be Any if F is not defined well
function RangeFun(foo::F) where {F<:Function}
    Tout = Base.return_types(foo, (Float64,))[1]
    RangeFun{Float64,Tout}(FunctionWrapper{Tout,Tuple{Float64}}(foo), MVector(0.0, 1.0))
end
function RangeFun(foo::F, range::Vector{T}) where {F<:Function} where {T}
    Tout = Base.return_types(foo, (T,))[1]
    RangeFun{T,Tout}(FunctionWrapper{Tout,Tuple{T}}(foo), MVector(range...))
end

function RangeFun(foo::F, range::Tuple{T,T}) where {F<:Function} where {T}
    Tout = Base.return_types(foo, (T,))[1]
    RangeFun{T,Tout}(FunctionWrapper{Tout,Tuple{T}}(foo), MVector(range...))
end

#TODO:: only the first elements of the Base.return_types, but is should not be a problem, because the inputstypes are well defined
function (fc::RangeFun{Tin,Tout})(t::Tin)::Tout where {Tin} where {Tout}
    @inbounds out = fc.f(t)
    return out::Tout
end


#--------------- Functions --------------------
#keep the range of the first argument! (Mostly it is assumned that the ranges are the same)
function Base.:+(v::RangeFun{Tin,Tout}, w::RangeFun{Tin,Tout})::RangeFun{Tin,Tout} where {Tin} where {Tout}
    return RangeFun(FunctionWrapper{Tout,Tuple{Tin}}((t::Tin) -> v(t) + w(t)), v.range)
end
#keep the range of the first argument! (Mostly it is assumned that the ranges are the same)
function Base.:-(v::RangeFun{Tin,Tout}, w::RangeFun{Tin,Tout})::RangeFun{Tin,Tout} where {Tin} where {Tout}
    return RangeFun(FunctionWrapper{Tout,Tuple{Tin}}((t::Tin) -> v(t) - w(t)), v.range)
end

#Base.:*(α, v): multiply v with a scalar α, which can be of a different scalar type; in particular this method is used to create vectors similar to v but with a different type of underlying scalars.
function Base.:*(α::N, v::RangeFun{Tin,Tout})::RangeFun{Tin,Tout} where {N,Tin,Tout}
    return RangeFun(FunctionWrapper{Tout,Tuple{Tin}}((t::Tin) -> α .* v(t)), v.range)
end
#Base.:*(α, v): multiply v with a scalar α, which can be of a different scalar type; in particular this method is used to create vectors similar to v but with a different type of underlying scalars.
function Base.:*(v::RangeFun{Tin,Tout}, α::N)::RangeFun{Tin,Tout} where {N,Tin,Tout}
    return RangeFun(FunctionWrapper{Tout,Tuple{Tin}}((t::Tin) -> α .* v(t)), v.range)
end
#foo = 10.0 * his
#foo([1,2], -1.3)


#Base.similar(v): a way to construct vectors which are exactly similar to v
function Base.similar(v::RangeFun{Tin,Tout})::RangeFun{Tin,Tout} where {N,Tin,Tout}
    # #TODO: return an error if v is also empty - maybe it is not a problem
    wout = 0.0 * deepcopy(v) #TODO: do we need to copy it?
    #empty!(v)
end
#TODO: this is  bug (or feature): foo([1,2], -1.3)

function Base.empty!(w::RangeFun{Tin,Tout})::RangeFun{Tin,Tout} where {Tin,Tout}
    w.f = FunctionWrapper{Tout,Tuple{Tin}}((t::Tin) -> zero(Tout))
    w.range .= MVector(zero(Tin), one(Tin))
    return w
end



#LinearAlgebra.mul!(w, v, α): out of place scalar multiplication; multiply vector v with scalar α and store the result in w
function LinearAlgebra.mul!(w::RangeFun{Tin,Tout}, v::RangeFun{Tin,Tout}, α::N)::RangeFun{Tin,Tout} where {Tin,Tout,N<:Number}
    #empty!(w)
    w.f = FunctionWrapper{Tout,Tuple{Tin}}((t::Tin) -> α .* v(t))
    w.range .= v.range
    return w
end


#LinearAlgebra.rmul!(v, α): in-place scalar multiplication of v with α; in particular with α = false, v is the corresponding zero vector
function LinearAlgebra.rmul!(v::RangeFun{Tin,Tout}, α::N)::RangeFun{Tin,Tout} where {Tin,Tout,N<:Number}
    v_foo=deepcopy(v.f)
    v.f = FunctionWrapper{Tout,Tuple{Tin}}((t::Tin) -> v_foo(t) .* α)
    return v
end

#LinearAlgebra.axpy!(α, v, w): store in w the result of α*v + w
function LinearAlgebra.axpy!(α::N, v::RangeFun{Tin,Tout}, w::RangeFun{Tin,Tout})::RangeFun{Tin,Tout} where {Tin,Tout,N<:Number}
    w_foo=deepcopy(w.f)
    w.f = FunctionWrapper{Tout,Tuple{Tin}}((t::Tin) -> α .* v(t)+ w_foo(t))
    return w
end

##TODO: ezzel lenne jó convert(::Type{funcomp}, x::F) where F<:Function= funcomp(x)::funcomp
#function LinearAlgebra.axpy!(α::Number, v::funcomp, fw::Function)::funcomp
#    LinearAlgebra.axpy!(α, v, funcomp(fw, v.range))
#end



#LinearAlgebra.axpby!(α, v, β, w): store in w the result of α*v + β*w
function LinearAlgebra.axpby!(α::N, v::RangeFun{Tin,Tout}, β::N, w::RangeFun{Tin,Tout})::RangeFun{Tin,Tout} where {Tin,Tout,N<:Number}
    w_foo=deepcopy(w.f)
    w.f = FunctionWrapper{Tout,Tuple{Tin}}((t::Tin) -> α .* v(t)+ β .*w_foo(t))
    return w
end



#####TODO: ehhez kell a struktúrába valami range (validity range-et definiálni)
####using Integrals
#####LinearAlgebra.dot(v,w): compute the inner product of two vectors
####function LinearAlgebra.dot(v::RangeFun{Tin,Tout},  w::RangeFun{Tin,Tout})::Tnorm where Tnorm
####    #TODO: nem az egészet kell összeintegrálni, 
####    #TODO: nem jók az idők, mert a függvény az csak a hsitory function. Csak a mappingben kell a sol és abból kivenni a cuccot...
####
####    tstart = maximum([v.range[1], w.range[1]])
####    tend = minimum([v.range[2], w.range[2]])
####
####    prob = IntegralProblem((t, p) -> v(t)' * w(t), tstart, tend)
####    VW = solve(prob, HCubatureJL(), reltol=1e-5, abstol=1e-5).u# 
####    #println(VW)
####    return VW::Tnorm
####end
####
#####LinearAlgebra.norm(v): compute the 2-norm of a vector
####function LinearAlgebra.norm(v::funcomp)::Float64
####    return sqrt(dot(v, v))::Float64
####end
