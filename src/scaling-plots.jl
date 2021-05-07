using CSV
using Plots
pgfplotsx()

function plot_random16(
    results=CSV.File(joinpath(@__DIR__, "..", "results", "random16.csv"));
    kwargs...
)
    dim = getproperty.(results, :dim)
    te = getproperty.(results, :te) ./ 1e9
    ti = getproperty.(results, :ti) ./ 1e9
    tc = getproperty.(results, :tc) ./ 1e9

    plt = plot(
        dim, tc;
        label="cohomology", markershape=:rect,
        legend=:topleft, xlabel="dim", ylabel="t[s]",
        title="",
        kwargs...,
    )
    plot!(plt, dim, te, label="Eirene", markershape=:circle, color=3)
    plot!(plt, dim, ti, label="involuted", markershape=:diamond, color=2)

    plt_log = plot(
        dim, tc;
        label="cohomology", markershape=:rect,
        legend=false, xlabel="dim", ylabel="t[s]",
        title="",
        yscale=:log10,
        kwargs...,
    )
    plot!(plt_log, dim, te, label="Eirene", markershape=:circle, color=3)
    plot!(plt_log, dim, ti, label="involuted", markershape=:diamond, color=2)
    return plot(plt, plt_log, size=(600, 300))
end

function plot_scaling(
    results=CSV.File(joinpath(@__DIR__, "..", "results", "gcycle.csv"));
    left_ylim=nothing,
    right_ylim=nothing,
    kwargs...,
)
    N = getproperty.(results, :N)
    te = getproperty.(results, :te) ./ 1e9
    ti = getproperty.(results, :ti) ./ 1e9
    tc = getproperty.(results, :tc) ./ 1e9
    if !isnothing(left_ylim)
        left_kwargs = (; kwargs..., ylim=left_ylim)
    else
        left_kwargs = kwargs
    end
    if !isnothing(right_ylim)
        right_kwargs = (; kwargs..., ylim=right_ylim)
    else
        right_kwargs = kwargs
    end

    plt = plot(
        N, tc;
        label="cohomology", markershape=:rect,
        legend=:topleft, xlabel="N", ylabel="t[s]",
        title="",
        left_kwargs...,
    )
    plot!(plt, N, te, label="Eirene", markershape=:circle, color=3)
    plot!(plt, N, ti, label="involuted", markershape=:diamond, color=2)

    plt_log = plot(
        N, tc;
        label="cohomology", markershape=:rect,
        legend=false, xlabel="N", ylabel="t[s]",
        title="",
        yscale=:log10,
        right_kwargs...,
    )
    plot!(plt_log, N, te, label="Eirene", markershape=:circle, color=3)
    plot!(plt_log, N, ti, label="involuted", markershape=:diamond, color=2)
    return plot(plt, plt_log, size=(600,300))
end
