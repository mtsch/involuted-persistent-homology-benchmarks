using BenchmarkTools

#=
println("preparing Eirene")
include(joinpath(@__DIR__, "eirene-comparison.jl"))
e = run(EireneComparison.suite())
println("running Eirene")
BenchmarkTools.save("results-eirene.json", e)

println("preparing Matrix")
include(joinpath(@__DIR__, "matrix-reduction.jl"))
m = run(MatrixReduction.suite())
println("running Matrix")
BenchmarkTools.save("results-matrix.json", m)

println("preparing Scaling")
include(joinpath(@__DIR__, "scaling.jl"))
run_gcycle_and_random16_scaling()
run_subset_scaling()
=#
