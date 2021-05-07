using CSV

# Utility file I/O functions. Not fit for general purpose use.
function load_dipha(filename)
    open(filename, "r") do f
        magic = read(f, Int)
        @assert magic == 8067171840
        type = read(f, Int)
        @assert type == 1
        len = read(f, Int)
        dim = read(f, Int)
        if dim == 2
            m = read(f, Int)
            n = read(f, Int)
            result = Array{Float64}(undef, (m, n))
            read!(f, result)

            return result
        elseif dim == 3
            m = read(f, Int)
            n = read(f, Int)
            o = read(f, Int)
            result = Array{Float64}(undef, (m, n, o))
            read!(f, result)

            return result
        else
            error("not implemented for dim=$dim")
        end
    end
end

function save_dipha(filename, data)
    open(filename, "w") do f
        write(f, 8067171840)
        write(f, 1)
        write(f, length(data))
        write(f, length(size(data)))
        for s in size(data)
            write(f, s)
        end
        for i in data
            write(f, Float64(i))
        end
    end
end

function load_points(filename)
    # ripser uses Float32
    file = CSV.File(filename, header=0, type=Float32, delim=' ')
    nrow = length(file)
    ncol = length(file[1])
    result = Vector{NTuple{ncol, Float32}}(undef, nrow)
    for i in 1:nrow
        result[i] = Tuple(file[i, :]...)
    end
    return result
end

function load_dist(filename)
    file = CSV.File(filename, header=0, type=Float32, delim=' ')
    result = zeros(Float32, (length(file), length(file)))
    for (i, row) in enumerate(file)
        result[i, :] .= collect(row)
    end
    return result
end

function save_sparse(filename, data)
    I, J, V = findnz(sparse(LowerTriangular(data)))
    return CSV.write(joinpath(filename), (I=I.-1, J=J.-1, V=V), writeheader=false, delim=' ')
end

function load_sparse(filename)
    table = CSV.File(filename, header=0, type=Float32, delim=' ')
    I = table[:, 1] .+ 1
    J = table[:, 2] .+ 1
    V = table[:, 3] .+ 1
    sparse([I; J], [J; I], [V; V])
end

function load_data(filename)
    ext = splitext(filename)[2]
    if ext == ".pts" || ext == ".alpha"
        return load_points(filename)
    elseif ext == ".dist"
        return load_dist(filename)
    elseif ext == ".spdist"
        return load_sparse(filename)
    elseif ext == ".dipha"
        return load_dipha(filename)
    else
        error("unsupported filetype \"$(ext)\"")
    end
end

using BenchmarkTools: prettytime

function print_table_eirene(res)
    for k in keys(res)
        r = res[k]
        println(
            k, " & & & & ", prettytime(minimum(r["cohomology"]).time),
            " & ", prettytime(minimum(r["involuted"]).time),
            " & ", prettytime(minimum(r["eirene"]).time),
            " & ", minimum(r["involuted"]).time / minimum(r["eirene"]).time
        )
    end
end

function print_table_matrix(res)
    for k in keys(res)
        r = res[k]
        if haskey(r, "homology")
            println(
                k, " & & & & ", prettytime(minimum(r["cohomology"]).time),
                " & ", prettytime(minimum(r["homology"]).time),
                " & ", prettytime(minimum(r["involuted"]).time),
            )
        else
            println(
                k, " & & & & ", prettytime(minimum(r["cohomology"]).time),
                " & & ", prettytime(minimum(r["involuted"]).time),
            )
        end
    end
end
