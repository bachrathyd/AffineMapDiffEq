5+5
includet("src\\fucntion_composition.jl")
#includet("fucntion_composition.jl")
#---------------- tests -------------------

foo1( t::Float64) = sin(t)
foo2( t::Float64) = 10.0
foo3( t::Float64) = t


his = funcomp(300.0, foo2, [-1.0, 10.0])
his = funcomp([300.0, 500], [foo1, foo2], [-1.0, 10.0])




foos = deepcopy(his)
empty!(foos)
foos
his


#TODO: this is  bug (or feature): foo([1,2], -1.3)
his = funcomp([300.0, 500], [foo1, foo2])
fooX = funcomp([5.0, 10.0, 20.0], [foo3, foo1, his], [-1.0, 0.5])
rmul!(fooX, 10.0)
mul!(fooX, his, 10.0)

his = funcomp([300.0, 500], [foo1, foo2])
for k in 1:5
    his = funcomp([5.0, 10.0, 20.0], [foo3, foo1, his], [-1.0, 0.5])#Already reducted during the construction
end
#@benchmark his(1.0)
#collapse!(his)
#@time his(1.0)
#reduce!(his)
#@time his(1.0)
#@benchmark his(1.0)

fc = fooX
foo = similar(his)

his = funcomp([300.0, 500], [foo1, foo2], [1.0, 2.0])
rmul!(his, 0.01)
his
his( -1.3)


his = funcomp([300.0, 500], [foo1, foo2], [-2.0, 0.5])
his2 = funcomp(foo3, [-1.0, 1.5])
axpy!(10.0, his, his2)
axpy!(10.0, his, foo2)



his = funcomp([300.0, 500], [foo1, foo2], [-2.0, 0.5])
his2 = funcomp(foo3, [-1.0, 1.5])
axpby!(10.0, his, 20.0, his2)
axpby!(10.0, his, 20.0, his2)


reduce!(his2)


his = funcomp([foo1, foo2], [-2.0, 0.5])
his1 = funcomp(foo3, [-0.0, 1.0])
his2 = funcomp(foo3, [-1.0, 1.5])

dot(his1, his2)

norm(his1)