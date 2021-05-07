using BenchmarkTools
using LightGraphs
using Ripserer
using Eirene
using Random

include(joinpath(@__DIR__, "utils.jl"))

function gcycle_scaling(Ns=[[10]; 100:100:1000]; seconds=1000, samples=10)
    map(Ns) do N
        dist = floyd_warshall_shortest_paths(cycle_graph(N)).dists

        println("N=$N")
        println("running cohomology")
        bc = @benchmark ripserer(
            $dist
        ) seconds=seconds samples=samples evals=1 gcsample=true
        tc = minimum(bc).time

        println("  running involuted ")
        bi = @benchmark ripserer(
            $dist; alg=:involuted
        ) seconds=seconds samples=samples evals=1 gcsample=true
        ti = minimum(bi).time

        println("  running eirene    ")
        be = @benchmark eirene(
            $dist; model="vr"
        ) seconds=seconds samples=samples evals=1 gcsample=true
        te = minimum(be).time
        println("  done.             ")

        println((;N, tc, ti, te))
        return (;N, tc, ti, te)
    end
end

function random16_scaling(dims=1:8)
    dist = load_data(joinpath(@__DIR__, "..", "datasets", "random16.dist"))

    map(dims) do dim
        println("dim=$dim")

        println("running cohomology")
        bc = @benchmark ripserer(
            $dist; dim_max=$dim
        ) seconds=1000 samples=10 evals=1 gcsample=true
        tc = minimum(bc).time

        println("  running involuted ")
        bi = @benchmark ripserer(
            $dist; dim_max=$dim, alg=:involuted
        ) seconds=1000 samples=10 evals=1 gcsample=true
        ti = minimum(bi).time

        println("  running eirene    ")
        be = @benchmark eirene(
            $dist; model=$("vr"), maxdim=$dim
        ) seconds=1000 samples=10 evals=1 gcsample=true
        te = minimum(be).time
        println("  done.             ")

        (;dim, tc, ti, te)
    end
end

function run_gcycle_and_random16_scaling()
    CSV.write("random16.csv", random16_scaling())
    CSV.write("gcycle.csv", gcycle_scaling())
end

"""
    subset_scaling(data, subsets; eirene_max=subsets[end], samples=10, dim_max=1)

Benchmark Ripserer and Eirene at various subsets of the `data`. `data` is shuffled first,
then increasingly large subsets are taken.
"""
function subset_scaling(
    data, subsets; eirene_max=subsets[end], samples=10, seconds=1000, dim_max=1, threshold=nothing
)
    data = perturb(data)
    if !isnothing(threshold)
        thr_r = (;threshold)
        thr_e = (;maxrad=threshold)
    else
        thr_r = NamedTuple()
        thr_e = NamedTuple()
    end
    map(subsets) do N
        ripserer_data = subset(data, N)
        if ripserer_data isa Vector
            eirene_data = Ripserer.to_matrix(ripserer_data)
            model="pc"
        else
            eirene_data = ripserer_data
            model="vr"
        end

        println("  running cohomology on $N points")
        bc = @benchmark ripserer(
            $ripserer_data, dim_max=$dim_max, $thr_r...
        ) seconds=seconds samples=samples evals=1 gcsample=true
        tc = minimum(bc).time

        println("  running involuted on $N points")
        bi = @benchmark ripserer(
            $ripserer_data; alg=:involuted, dim_max=$dim_max, $thr_r...
        ) seconds=seconds samples=samples evals=1 gcsample=true
        ti = minimum(bi).time

        if N ≤ eirene_max
            println("  running eirene on $N points")
            be = @benchmark eirene(
                $eirene_data; model=$model, maxdim=$dim_max, $thr_e...
            ) seconds=seconds samples=samples evals=1 gcsample=true
            te = minimum(be).time
        else
            te = nothing
        end

        (;N, tc, ti, te)
     end
end

"""
    dim_scaling(data, dims; eirene_max, samples=10)

Run Ripserer and Eirene on the same data in various dimensions.
"""
function dim_scaling(data, dims; eirene_max, samples=10, seconds=1000, threshold=nothing)
    if !isnothing(threshold)
        thr_r = (;threshold)
        thr_e = (;maxrad=threshold)
    else
        thr_r = NamedTuple()
        thr_e = NamedTuple()
    end
    map(dims) do dim_max
        model = data isa AbstractMatrix ? "vr" : "pc"
        bc = @benchmark ripserer(
            $data; dim_max=$dim_max, $thr_r...
        ) seconds=seconds samples=samples evals=1 gcsample=true
        tc = minimum(bc).time

        println("  running involuted ")
        bi = @benchmark ripserer(
            $data; alg=:involuted, dim_max=$dim_max, $thr_r...
        ) seconds=seconds samples=samples evals=1 gcsample=true
        ti = minimum(bi).time

        println("  running eirene    ")
        be = @benchmark eirene(
            $data; model=$model, maxdim=$dim_max, $thr_e...
        ) seconds=seconds samples=samples evals=1 gcsample=true
        te = minimum(be).time

        (;N, tc, ti, te)
     end
end

perturb(v::AbstractVector) = shuffle(v)
function perturb(m::AbstractMatrix)
    ord = shuffle(1:size(m, 1))
    return m[ord, ord]
end
subset(v::AbstractVector, N) = v[1:N]
subset(m::AbstractMatrix, N) = m[1:N, 1:N]

function run_subset_scaling()
    for (file,              name,         dim_max, threshold, subsets) in [
        ("dragon1000.dist", "dragon1000", 1,       nothing,   100:100:1000),
        ("o3_1024.pts",     "o3_1024",    3,       1.8,       124:100:1024),
        ("o3_4096.pts",     "o3_4096",    3,       1.4,       196:100:4096),
        ("hiv.dist",        "hiv",        1,       nothing,   188:100:1088),
    ]
        println(" ++ $name ++")
        data = load_data(joinpath(@__DIR__, "..", "datasets", file))
        res = subset_scaling(
            data, subsets; eirene_max=subsets[end], samples=10, dim_max=dim_max
        )
        CSV.write("$name.csv", res)
    end
end
