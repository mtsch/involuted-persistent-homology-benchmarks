module EireneComparison

using BenchmarkTools
using Eirene
using Ripserer

include(joinpath(@__DIR__, "utils.jl"))

function _eirene(data; dim_max=1, threshold=nothing, model)
    if isnothing(threshold)
        return eirene(data; maxdim=dim_max, model=model)
    else
        return eirene(data; maxdim=dim_max, maxrad=threshold, model=model)
    end
end

function add_suite!(suite, data, name, dim, threshold; seconds=100, samples=10)
    msg = " preparing $name ($dim). "
    println("="^length(msg))
    println(msg)
    suite[name] = BenchmarkGroup()

    rkw=(;dim_max=dim, reps=true)
    ekw=(;maxdim=dim)
    if !isnothing(threshold)
        rkw = (; rkw..., threshold)
        ekw = (; ekw..., maxrad=threshold)
    end
    if data isa Vector
        eirene_data = Ripserer.to_matrix(data)
        ekw = (;ekw..., model="pc")
    else
        eirene_data = data
        ekw = (;ekw..., model="vr")
    end
    suite[name]["cohomology"] = @benchmarkable ripserer($data; $rkw...) setup=(
        print("  cohomology for ", $name, "                ")
    ) samples=samples seconds=seconds gcsample=true

    suite[name]["involuted"] = @benchmarkable ripserer($data; $rkw..., alg=:involuted) setup=(
        print("  involuted for ", $name, "                ")
    ) samples=samples seconds=seconds gcsample=true

    suite[name]["eirene"] = @benchmarkable eirene($eirene_data; $ekw...) setup=(
        print("  eirene for ", $name, "                ")
    ) samples=samples seconds=seconds gcsample=true

    return suite
end

"""
    EireneComparison.suite()

The suite will compare the computation time of persistent homology in the `dim`-th dimension
for a `Rips` filtration built on each dataset in `datasets`. Notes:

* Ripserer is benchmarked for both persistent cohomology and involuted persistent homology.
* When the data set is a point cloud, calculating the distance matrix is also included in
  the benchmarks. This time is expected to be negligible, but it is done in all cases in the
  interest of fairness.
* Ripserer is set to compute representatives in both cases. The effect this has on
  performance is minor, but is kept in the interest of fairness.
* Eirene will outperform ripserer on `rand` benchmarks because it can do collapses on the
  complex. Ripserer can also perform edge collapses, but this was left out to keep the
  number of benchmarks lower.
* Some Ripserer times could be improved by using a sparse filtration. This was left out.
"""
function suite(
    datasets=[
        # filename            name           dim threshold
        ("gcycle.dist",       "gcycle (1)",   1,  nothing),
        ("gcycle.dist",       "gcycle (3)",   3,  40),
        ("2_sphere_100.pts",  "2 sphere",     2,  nothing),
        ("rand_4_50.pts",     "rand_4",       3,  nothing),
        ("random16.dist",     "random16 (3)", 3,  nothing),
        ("celegans.dist",     "celegans",     2,  nothing),
        ("dragon1000.dist",   "dragon1000",   1,  nothing),
        ("o3_1024.pts",       "o3_1024",      3,  1.8),
        ("o3_4096.pts",       "o3_4096",      3,  1.4),
        ("random16.dist",     "random16 (7)", 7,  nothing),
        ("hiv.dist",          "hiv",          1,  nothing),
    ];
    seconds=1000, samples=10,
)
    suite = BenchmarkGroup()
    for (file, name, dim, threshold) in datasets
        data = load_data(joinpath(@__DIR__, "..", "datasets", file))
        add_suite!(suite, data, name, dim, threshold; seconds, samples)
    end
    return suite
end

end
