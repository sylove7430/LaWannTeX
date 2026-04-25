module LatticeIO

using TOML

using ..Model: LatexConfig
using ..PoscarIO: read_poscar
using ..WannierWinIO: read_win_lattice

export resolve_lattice

function _parse_lattice(raw, key::AbstractString)::Matrix{Float64}
    raw isa AbstractVector || error("$key must be a 3x3 nested array")
    length(raw) == 3 || error("$key must have three basis vectors")

    lattice = Matrix{Float64}(undef, 3, 3)
    for col in 1:3
        vec = raw[col]
        vec isa AbstractVector || error("Each lattice vector in $key must be an array")
        length(vec) == 3 || error("Each lattice vector in $key must have length 3")
        lattice[:, col] = Float64.(vec)
    end
    return lattice
end

function _extract_structure_table(cfg)
    structure_tbl = get(cfg, "structure", nothing)
    if structure_tbl isa AbstractDict && haskey(structure_tbl, "lattice")
        return structure_tbl
    end
    return nothing
end

function _read_lattice_toml(path::AbstractString)
    cfg = TOML.parsefile(path)
    structure_tbl = _extract_structure_table(cfg)
    isnothing(structure_tbl) && error("Missing [structure].lattice in $(abspath(path))")
    return _parse_lattice(structure_tbl["lattice"], "structure.lattice")
end

function _read_lattice(path::AbstractString)
    if lowercase(splitext(path)[2]) == ".toml"
        return _read_lattice_toml(path)
    end
    return read_poscar(path).lattice
end

function _read_inline_lattice(path::AbstractString)
    lowercase(splitext(path)[2]) == ".toml" || return nothing

    cfg = TOML.parsefile(path)
    isnothing(_extract_structure_table(cfg)) && return nothing
    return _read_lattice_toml(path)
end

function resolve_lattice(cfg::LatexConfig)
    try
        lattice = _read_inline_lattice(cfg.input_path)
        !isnothing(lattice) && return lattice
    catch
    end

    if !isnothing(cfg.win_path) && isfile(cfg.win_path)
        try
            return read_win_lattice(cfg.win_path)
        catch
        end
    end

    if !isnothing(cfg.structure_path) && isfile(cfg.structure_path)
        try
            return _read_lattice(cfg.structure_path)
        catch
        end
    end

    return nothing
end

end
