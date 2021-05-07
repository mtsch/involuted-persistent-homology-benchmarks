using LinearAlgebra
using Plots
using Random
using Ripserer
pgfplotsx()

function circ_rand(n, r, center)
    map(rand(n)) do φ
        (r * sin(2π * φ), r * cos(2π * φ)) .+ center
    end
end

function circ_equi(n, r, center)
    map(range(0, 1, length=n+1)[1:end-1]) do φ
        (r * sin(2π * φ), r * cos(2π * φ)) .+ center
    end
end

function annulus(n, r1, r2, center)
    points = Tuple{Float64,Float64}[]
    while length(points) < n
        v = 2r1 * (rand(2) .- 0.5)
        if r1 ≥ norm(v) > r2
            push!(points, (v[1], v[2]) .+ center)
        end
    end
    return points
end

Random.seed!(1337)
data = [
    circ_rand(20, 1, (0, 2));
    circ_equi(12, 1.5, (0, -2));
    annulus(40, 2, 1, (2 * √3, 0));
]

res_co = ripserer(data; reps=true)[2][end-2:end]
res_ho = ripserer(data; alg=:involuted)[2][end-2:end]

function plot_representatives(data, ints; kwargs...)
    plt = scatter(data; legend=false, aspect_ratio=1, markersize=2, kwargs...)
    foreach(ints) do int
        plot!(plt, int, data; threshold_strict=death(int))
    end
    return plt
end

plot(
    plot_representatives(data, res_co; title="Representative Cocycles"),
    plot_representatives(data, res_ho; title="Representative Cycles"),
    size=(600, 300),
)
savefig(joinpath(@__DIR__, "figures/cycle-cocycle-comparison.pdf"))
