module MatrixReduction

using BenchmarkTools
using Ripserer
using Ripserer: AbstractRipsFiltration, simplex_type, BoundaryMatrix, CoboundaryMatrix, compute_death_simplices!, compute_intervals!
using IterTools

include(joinpath(@__DIR__, "utils.jl"))

"""
    collect_simplices(flt::AbstractRipsFiltration, dim)

Collect all simplices in `dim`-th dimension of Rips-like filtration.
"""
function collect_simplices(flt::AbstractRipsFiltration, dim)
    res = simplex_type(flt, dim)[]
    for vs in IterTools.subsets(vertices(flt), Val(dim + 1))
        σ = simplex(flt, Val(dim), vs)
        !isnothing(σ) && push!(res, σ)
    end
    return res
end

function coboundary_matrix(flt, dim_max; field=Mod{2}, implicit=true)
    _, to_reduce, to_skip = Ripserer.zeroth_intervals(
        flt, 0, false, field, false
    )
    matrix = CoboundaryMatrix{implicit}(field, flt, to_reduce, to_skip)
    for dim in 1:(dim_max - 1)
        compute_intervals!(matrix, 0, false, false)
        matrix = Ripserer.next_matrix(matrix, false)
    end
    return matrix
end

function boundary_matrix(flt, dim; field=Mod{2}, implicit=false)
    columns = collect_simplices(flt, dim+1)
    sort!(columns)
    return BoundaryMatrix{implicit}(field, flt, columns)
end

function involuted_boundary_matrix(flt, dim; field=Mod{2}, implicit=false)
    d = coboundary_matrix(flt, dim; field)
    columns, _ = compute_death_simplices!(d, false, 0)
    return BoundaryMatrix{implicit}(field, flt, columns)
end

"""
    reduce!(matrix; verbose=false, reps=true)

Reduce a matrix. Note: this functions is mutating! Don't run it twice on the same matrix.
"""
function reduce!(matrix; verbose=false, reps=true)
    compute_intervals!(matrix, 0, verbose, reps)
end

"""
    add_suite!(suite, data, name, dim; seconds=100, samples=10)

Add a suite to the `BenchmarkGroup`.
"""
function add_suite!(suite, data, name, dim, threshold, hom; seconds=100, samples=10)
    msg = " preparing $name ($dim). "
    println("="^length(msg))
    println(msg)
    suite[name] = BenchmarkGroup()

    if isnothing(threshold)
        ℜ = Rips(data)
    else
        ℜ = Rips(data; threshold)
    end
    d = coboundary_matrix(ℜ, dim)
    println("  d has ", length(d.columns_to_reduce), " columns.")
    if hom
        ∂ = boundary_matrix(ℜ, dim)
        println("  ∂ has ", length(∂.columns_to_reduce), " columns.")
    end
    ∂′= involuted_boundary_matrix(ℜ, dim)
    println("  ∂′ has ", length(∂′.columns_to_reduce), " columns.")

    suite[name]["cohomology"] = @benchmarkable reduce!(M) setup=(
        print("  reducing d for ", $name, "                    "); M=deepcopy($d)
    ) samples=samples seconds=seconds gcsample=true

    if hom
        suite[name]["homology"] = @benchmarkable reduce!(M) setup=(
            print("  reducing ∂ for ", $name, "                    "); M=deepcopy($∂)
        ) samples=samples seconds=seconds gcsample=true
    end

    suite[name]["involuted"] = @benchmarkable reduce!(M) setup=(
        print("  reducing ∂ for ", $name, "                    "); M=deepcopy($∂′)
    ) samples=samples seconds=seconds gcsample=true

    return suite
end

"""
    MatrixReduction.suite()

The suite will benchmark the reduction of a matrix in the `dim`-th dimension for a `Rips`
filtration built on each dataset in `datasets`. Notes:

* All (co)boundary matrices are built on the fly.
* The coboundary matrix uses an implicit representation of reduced columns, while the
  representation of boundary matrices is more explicit (they're reduced faster that way).
* The coboundary matrix uses the twist algorithm to select which columns need to be
  reduced. Without it, coboundary matrix reduction would be orders of magnitude slower (and
  the timings would not reflect real-world scenarios)
* The boundary matrix is the full matrix. The assumption here is that `dim` is the highest
  dimension we are computing in, so the twist algorithm would not help.
* The involuted matrix sometimes has fewer columns than the coboundary matrix because it is
  only collected up to the time of the last non-zero interval.
* The time it takes to collect the simplices is not accounted for. Sorting the simplices
  is. Including this time would work in favor of cohomology.
* Due to how Ripserer is implemented, 1-dimensional columns in the coboundary matrix are not
  sorted.
"""
function suite(
    datasets=[
        # filename            name           dim threshold homology?
        ("gcycle.dist",       "gcycle (1)",   1,  nothing,  true),
        ("gcycle.dist",       "gcycle (3)",   3,  40,       true),
        ("2_sphere_100.pts",  "2 sphere",     2,  nothing,  true),
        ("rand_4_50.pts",     "rand_4",       3,  nothing,  true),
        ("random16.dist",     "random16 (3)", 3,  nothing,  true),
        ("celegans.dist",     "celegans",     2,  nothing,  true),
        ("dragon1000.dist",   "dragon1000",   1,  nothing,  true),
        ("o3_1024.pts",       "o3_1024",      3,  1.8,      false),
        ("o3_4096.pts",       "o3_4096",      3,  1.4,      false),
        ("random16.dist",     "random16 (7)", 7,  nothing,  false),
        ("hiv.dist",          "hiv",          1,  nothing,  true),
    ];
    seconds=1000, samples=10,
)
    suite = BenchmarkGroup()
    for (file, name, dim, threshold, hom) in datasets
        data = load_data(joinpath(@__DIR__, "..", "datasets", file))
        add_suite!(suite, data, name, dim, threshold, hom; seconds, samples)
    end
    return suite
end

end
