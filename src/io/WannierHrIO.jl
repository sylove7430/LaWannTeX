module WannierHrIO

using ..WannierTypes: HrBlocks, RKey

export read_hr

function _get_hopping!(dict::Dict{RKey, Matrix{ComplexF64}}, R::RKey, nw::Int)
    return get!(dict, R) do
        zeros(ComplexF64, nw, nw)
    end
end

function _skip_ndegen(io::IO, nrpts::Int)
    seen = 0
    while seen < nrpts
        eof(io) && error("Unexpected EOF while skipping ndegen")
        line = strip(readline(io))
        isempty(line) && continue
        seen += length(split(line))
    end
    return nothing
end

function _parse_hr_line(line::AbstractString)
    fields = split(strip(line))
    length(fields) == 7 || error("Malformed hr line: $line")
    R = (parse(Int, fields[1]), parse(Int, fields[2]), parse(Int, fields[3]))
    m = parse(Int, fields[4])
    n = parse(Int, fields[5])
    value = ComplexF64(parse(Float64, fields[6]), parse(Float64, fields[7]))
    return R, m, n, value
end

function read_hr(path::AbstractString)::HrBlocks
    open(path, "r") do io
        eof(io) && error("Empty hr file: $(abspath(path))")
        header = strip(readline(io))
        num_wann = parse(Int, strip(readline(io)))
        nrpts = parse(Int, strip(readline(io)))
        _skip_ndegen(io, nrpts)

        hoppings = Dict{RKey, Matrix{ComplexF64}}()

        for ir in 1:nrpts
            block_R = nothing
            for _ in 1:(num_wann * num_wann)
                eof(io) && error("Unexpected EOF while reading hopping block $ir from $(abspath(path))")
                R, m, n, value = _parse_hr_line(readline(io))
                if isnothing(block_R)
                    block_R = R
                    haskey(hoppings, R) && error("Duplicate hopping block for R=$R in $(abspath(path))")
                elseif R != block_R
                    error("Encountered mixed R block in $(abspath(path)); expected $block_R, got $R")
                end
                H = _get_hopping!(hoppings, R, num_wann)
                H[m, n] = value
            end
        end

        return HrBlocks(header, num_wann, nrpts, hoppings)
    end
end
end
