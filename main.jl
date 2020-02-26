include("jsr.jl")

A = [1 0; 1 0]
B = [0 1; 0 1]

v = Matrix{Float64}[]
push!(v, A)
push!(v, B)

lb, ub = jsr(v)

println("[", lb, ", ", ub, "]")


