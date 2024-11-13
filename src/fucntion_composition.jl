#(history)function compatibility for Krylovkit - incomplete

#julia> myfun(s1::String, s2::String) = invoke(myfun, tuple(Any, Any), s1, s2) # overwrite
#WARNING: Method definition myfun(String, String) [...] overwritten
#https://discourse.julialang.org/t/overwriting-functions/404



using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Integrals


#h(p, t::Float64) = ...
#h(p, t::Float64, deriv::Type{Val{1}}) = ...
struct funcomp{T,Tout}#<: Function   F<:Function
    scaler::Vector{T}
    f::Vector{Function}
    range::Vector{T}
    #outtype=Tout::DataType
end

function Base.show(io::IO, fc::funcomp) # I have to remove <: Function to have pretty-print
    println(io, "scaler: $(fc.scaler)")
    println(io, "functinos: $(fc.f)")
    println(io, "range $(fc.range)")
end

function Base.show(io::IO, fc::funcomp) # I have to remove <: Function to have pretty-print
    println(io, typeof(fc))
    println(io, "scaler: $(fc.scaler)")
    println(io, "Num. functinos: $(size(fc.f,1))")
    println(io, "range $(fc.range)")
end


# ------------- Initializations --------------------
#error("provide no type")
#function funcomp()
#    funcomp{Float64}(Float64[], funcomp[], [0.0, 0.0])
#end
function Base.return_types(fc::funcomp{A,Tout}) where {A,Tout}
    return [Tout]
end
function funcomp(foo::F) where {F<:Function}
    Tout = Base.return_types(foo, (Float64,))[1]
    funcomp{Float64,Tout}([1.0], [foo], [0.0, 1.0])
end
function funcomp(foo::F, range::Vector{T}) where {F<:Function} where {T}
    Tout = Base.return_types(foo, (Float64,))[1]
    funcomp{Float64,Tout}([1.0], [foo], range)
end

function funcomp(scaler::Float64, foo::F) where {F<:Function}
    Tout:Base.return_types(foo, (Float64,))[1]
    funcomp{Float64,Tout}([scaler], [foo], [0.0, 1.0])
end
function funcomp(scaler::Float64, foo::F, range::Vector{T}) where {F<:Function} where {T}
    Tout = Base.return_types(foo, (Float64,))[1]
    funcomp{Float64,Tout}([scaler], [foo], range)
end

function funcomp(scaler::Vector{Float64,}, foos::Vector)# where {F<:Function}
    Tout = Base.return_types(foos[1], (Float64,))[1]
    reduce!(funcomp{Float64,Tout}(scaler, Any[foos...], [0.0, 1.0]))
end

function funcomp(scaler::Vector{Float64}, foos::Vector{F}, range::Vector{T}) where {F<:Function} where {T}
    Tout = Base.return_types(foos[1], (Float64,))[1]
    reduce!(funcomp{Float64,Tout}(scaler, Any[foos...], range))
end
function funcomp(scaler::Vector{Float64}, foos::Vector{F}, range::Vector{T}) where {F<:Any} where {T}
    Tout = Base.return_types(foos[1], (Float64,))[1]
    reduce!(funcomp{Float64,Tout}(scaler, Any[foos...], range))
end

function funcomp(foos::Vector)
    println("aaaaabbbbbaaaa")
    Tout = Base.return_types(foos[1], (Float64,))[1]
    reduce!(funcomp{Float64,Tout}(ones(Float64, size(foos, 1)), Any[foos...], [0.0, 1.0]))

end
function funcomp(foos::Vector, range::Vector{T}) where {T}
    Tout = Base.return_types(foos[1], (Float64,))[1]
    reduce!(funcomp{Float64,Tout}(ones(Float64, size(foos, 1)), Any[foos...], range))
end

#TODO:: All function must have the same output, it should be checked!
#TODO:: only the first elements of the Base.return_types, but is should not be a problem, because the inputstypes are well defined
#Ettől sokkal lassabb lett ha megadtam a kimenet tipusát: ::Base.return_types(fc.f[1],(PT,Float64))[1]
function (fc::funcomp{T,Tout})(t::Float64)::Tout where {T} where {Tout}
    @inbounds out = mapreduce(x -> x[1] .* x[2](t), +, zip(fc.scaler, fc.f))
    return out::Tout
end


#--------------- Functions --------------------
#keep the range of the first argument! (Mostly it is assumned that the ranges are the same)
function Base.:+(v::funcomp, w::funcomp)::funcomp
    return funcomp([v.scaler..., w.scaler...], [v.f..., w.f...], v.range)
end
function add!!(v::funcomp{T,Tout}, w::funcomp{T,Tout})::funcomp{T,Tout} where {T,Tout}
    return v+W
end
#keep the range of the first argument! (Mostly it is assumned that the ranges are the same)
function Base.:-(v::funcomp, w::funcomp)::funcomp
    return funcomp([v.scaler..., -w.scaler...], [v.f..., w.f...], v.range)
end

#Base.:*(α, v): multiply v with a scalar α, which can be of a different scalar type; in particular this method is used to create vectors similar to v but with a different type of underlying scalars.
function Base.:*(α::N, v::funcomp)::funcomp where {N<:Number}
    return funcomp(α .* v.scaler, v.f, v.range)
end
#foo = 10.0 * his
#foo([1,2], -1.3)


function Base.eltype(v::funcomp{T,Tout}) where {T,Tout}
    # #TODO: return an error if v is also empty - maybe it is not a problem
    return Tout
end
#Base.similar(v): a way to construct vectors which are exactly similar to v
function Base.similar(v::funcomp,element_type=eltype(v))::funcomp
    # #TODO: return an error if v is also empty - maybe it is not a problem
    wout = 0.0 * deepcopy(v)
    #empty!(v)
end
#TODO: this is  bug (or feature): foo([1,2], -1.3)

function zerovector!(v::funcomp{T,Tout}) where {T,Tout}
    return empty!(v)
end

function scalartype(v::funcomp{T,Tout}) where {T,Tout}
    return Tout
end

function scalartype(DT::Type{funcomp{T,Tout}}) where {T,Tout}
    return Tout
end
function Base.empty!(w::funcomp)::funcomp
    deleteat!(w.f, eachindex(w.f))#TODO:it is might be not necessary, but it is safe
    deleteat!(w.scaler, eachindex(w.scaler))#TODO:it is might be not necessary, but it is safe
    #w.range .= 0.0
    return w
end



#LinearAlgebra.mul!(w, v, α): out of place scalar multiplication; multiply vector v with scalar α and store the result in w
function LinearAlgebra.mul!(w::funcomp, v::funcomp, α::Number)::funcomp
    empty!(w)
    append!(w.f, v.f)
    append!(w.scaler, α .* v.scaler)
    w.range .= v.range
    return w
end


#LinearAlgebra.rmul!(v, α): in-place scalar multiplication of v with α; in particular with α = false, v is the corresponding zero vector
function LinearAlgebra.rmul!(v::funcomp, α::Number)::funcomp
    v.scaler .*= α
    return v::funcomp
end

#LinearAlgebra.axpy!(α, v, w): store in w the result of α*v + w
function LinearAlgebra.axpy!(α::Number, v::funcomp, w::funcomp)::funcomp

    # println("..............")
    # println(α)
    #
    # println("-----------")
    # println(size(w.f, 1))
    # println(size(v.f, 1))
    # println("-----------")
    # println(size(w.scaler, 1))
    # println(size(v.scaler, 1))
    # println(size(w.range, 1))
    # println(size(v.range, 1))


    append!(w.f, v.f)
    append!(w.scaler, α .* v.scaler)
    w.range[1] = maximum([v.range[1], w.range[1]])
    w.range[2] = minimum([v.range[2], w.range[2]])



    return reduce!(w)
end

#TODO: ezzel lenne jó convert(::Type{funcomp}, x::F) where F<:Function= funcomp(x)::funcomp
function LinearAlgebra.axpy!(α::Number, v::funcomp, fw::Function)::funcomp
    LinearAlgebra.axpy!(α, v, funcomp(fw, v.range))
end



#LinearAlgebra.axpby!(α, v, β, w): store in w the result of α*v + β*w
function LinearAlgebra.axpby!(α::Number, v::funcomp, β::Number, w::funcomp)::funcomp
    reduce!(axpy!(α, v, rmul!(w, β)))
end


## multi-level funcomp is flattended out
function collapse!(fc::funcomp)::funcomp
    for i in length(fc.f):-1:1
        if typeof(fc.f[i]) <: funcomp
            for kf in eachindex(fc.f[i].f)
                #push!(fc.f, pop!(fc.f[i].f))
                #push!(fc.scaler, fc.scaler[i] .* pop!(fc.f[i].scaler))
                push!(fc.f, fc.f[i].f[kf])
                push!(fc.scaler, fc.scaler[i] .* fc.f[i].scaler[kf])
            end

            fc.range[1] = maximum([fc.range[1], fc.f[i].range[1]])
            fc.range[2] = minimum([fc.range[2], fc.f[i].range[2]])

            deleteat!(fc.f, i)
            deleteat!(fc.scaler, i)
        end
    end
    return fc
end

# function compositions are reduced, usefull if the same function is apperas multiple times
function reduce!(fc::funcomp)::funcomp

    while any(typeof.(fc.f) .<: funcomp)
        collapse!(fc)
    end
    fsunique = unique(fc.f)

    newscales = [sum(fc.scaler[fc.f.==fu]) for fu in fsunique]

    empty!(fc)#range is not deleted!!!

    append!(fc.f, fsunique)
    append!(fc.scaler, newscales)

    return fc
end



#TODO: ehhez kell a struktúrába valami range (validity range-et definiálni)

#LinearAlgebra.dot(v,w): compute the inner product of two vectors
function LinearAlgebra.dot(v::funcomp, w::funcomp)::Float64
    #TODO: nem az egészet kell összeintegrálni, 
    #TODO: nem jók az idők, mert a függvény az csak a hsitory function. Csak a mappingben kell a sol és abból kivenni a cuccot...

    tstart = maximum([v.range[1], w.range[1]])
    tend = minimum([v.range[2], w.range[2]])

    #prob = IntegralProblem((t, p) -> v(t)' * w(t), tstart, tend)
    prob = IntegralProblem((t, p) -> dot(v(t), w(t)), tstart, tend)
    VW = solve(prob, HCubatureJL(), reltol=1e-5, abstol=1e-5).u# 
    #println(VW)
    return VW::Float64
end

#LinearAlgebra.norm(v): compute the 2-norm of a vector
function LinearAlgebra.norm(v::funcomp)::Float64
    return sqrt(dot(v, v))::Float64
end
