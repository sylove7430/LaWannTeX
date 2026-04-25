include(joinpath(@__DIR__, "InputParsing.jl"))

module InputIO

using TOML

using ..InputParsing: optional_int, optional_string, parse_vec3_float, required_bool
using ..InputParsing: required_float, required_string, required_string_vector, required_table, resolve_path
using ..LatexConstants: DEFAULT_ATOL, DEFAULT_INCLUDE_SECOND_QUANT, expand_orbital_symbol
using ..Model: LatexConfig, LegacySelection, OrbitalSelection, OrbitalSpec, SelectionMode

export read_input

const ALLOWED_TOP_LEVEL_KEYS = ("files", "render", "selection", "structure", "orbitals")
const ALLOWED_FILE_KEYS = ("hr", "win", "structure", "out")
const ALLOWED_RENDER_KEYS = ("atol", "max_rblocks", "include_second_quant")
const ALLOWED_SELECTION_KEYS = ("mode", "wannier_indices", "max_nn", "distance_tol")
const ALLOWED_STRUCTURE_KEYS = ("lattice",)
const ALLOWED_ORBITAL_KEYS = ("name", "position_frac", "orbitals", "spins")

function _reject_unknown_keys(tbl, allowed_keys, context::AbstractString)
    for raw_key in keys(tbl)
        key = String(raw_key)
        key in allowed_keys && continue
        error("$context has an unknown key \"$key\"")
    end
    return nothing
end

function _build_orbital_label(name::AbstractString, orbital::AbstractString, spin::AbstractString)::String
    spin_text = lowercase(strip(String(spin)))
    if spin_text in ("", "none", "spinless")
        return "$(name):$(orbital)"
    end
    return "$(name):$(orbital):$(spin_text)"
end

function _parse_input_orbital_specs(cfg)::Union{Nothing, Vector{OrbitalSpec}}
    raw = get(cfg, "orbitals", nothing)
    isnothing(raw) && return nothing
    raw isa AbstractVector || error("Use repeated [[orbitals]] tables.")
    isempty(raw) && error("[[orbitals]] cannot be empty.")

    specs = OrbitalSpec[]
    next_index = 1
    for (idx, item) in enumerate(raw)
        context = "[[orbitals]] #$idx"
        item isa AbstractDict || error("$context must be a table")
        _reject_unknown_keys(item, ALLOWED_ORBITAL_KEYS, context)

        name = required_string(item, "name"; context=context)
        position_frac = haskey(item, "position_frac") ? parse_vec3_float(item["position_frac"], "$context.position_frac") : nothing
        orbitals = required_string_vector(item, "orbitals"; context=context)
        spins = haskey(item, "spins") ? required_string_vector(item, "spins"; context=context) : ["none"]

        for raw_orb in orbitals
            for orbital in expand_orbital_symbol(raw_orb)
                for spin in spins
                    push!(
                        specs,
                        OrbitalSpec(
                            next_index,
                            _build_orbital_label(name, orbital, spin),
                            position_frac,
                        ),
                    )
                    next_index += 1
                end
            end
        end
    end

    return specs
end

function _normalized_unique_indices(values::Vector{String})::Vector{Int}
    indices = Int[]
    seen = Set{Int}()
    for (idx, value) in enumerate(values)
        parsed = try
            parse(Int, value)
        catch
            error("[selection]: wannier_indices must contain integers (bad entry at index $idx)")
        end
        parsed > 0 || error("[selection]: wannier_indices must be positive")
        parsed in seen && error("[selection]: wannier_indices must not contain duplicates")
        push!(indices, parsed)
        push!(seen, parsed)
    end
    sort!(indices)
    return indices
end

function _parse_selection(cfg, files_tbl, render_tbl, base_dir)::SelectionMode
    selection_tbl = get(cfg, "selection", Dict{String, Any}())
    selection_tbl isa AbstractDict || error("[selection] must be a table")
    _reject_unknown_keys(selection_tbl, ALLOWED_SELECTION_KEYS, "[selection]")

    mode = lowercase(optional_string(selection_tbl, "mode"; default="legacy", context="[selection]"))
    if mode == "legacy"
        max_rblocks = optional_int(render_tbl, "max_rblocks")
        return LegacySelection(max_rblocks)
    elseif mode == "orbital"
        haskey(selection_tbl, "wannier_indices") || error("[selection]: Missing required key \"wannier_indices\"")
        haskey(selection_tbl, "max_nn") || error("[selection]: Missing required key \"max_nn\"")
        haskey(render_tbl, "max_rblocks") && error("[render]: max_rblocks is legacy-only and cannot be used with selection.mode = \"orbital\"")

        raw_indices = selection_tbl["wannier_indices"]
        raw_indices isa AbstractVector || error("[selection]: wannier_indices must be an array")
        indices_as_text = String[string(item) for item in raw_indices]
        wannier_indices = _normalized_unique_indices(indices_as_text)

        max_nn = optional_int(selection_tbl, "max_nn"; context="[selection]")
        isnothing(max_nn) && error("[selection]: max_nn must be an integer")
        something(max_nn) >= 0 || error("[selection]: max_nn must be >= 0")

        distance_tol = haskey(selection_tbl, "distance_tol") ? required_float(selection_tbl, "distance_tol"; context="[selection]") : 1e-5

        return OrbitalSelection(
            wannier_indices,
            something(max_nn),
            distance_tol,
        )
    end

    error("[selection]: mode must be \"legacy\" or \"orbital\"")
end

function read_input(path::AbstractString)::LatexConfig
    cfg = TOML.parsefile(path)
    input_path = abspath(path)
    base_dir = dirname(input_path)
    _reject_unknown_keys(cfg, ALLOWED_TOP_LEVEL_KEYS, "top-level TOML")

    files_tbl = required_table(cfg, "files")
    render_tbl = get(cfg, "render", Dict{String, Any}())
    render_tbl isa AbstractDict || error("[render] must be a table")

    _reject_unknown_keys(files_tbl, ALLOWED_FILE_KEYS, "[files]")
    _reject_unknown_keys(render_tbl, ALLOWED_RENDER_KEYS, "[render]")

    if haskey(cfg, "structure")
        structure_tbl = cfg["structure"]
        structure_tbl isa AbstractDict || error("[structure] must be a table")
        _reject_unknown_keys(structure_tbl, ALLOWED_STRUCTURE_KEYS, "[structure]")
    end

    hr_path = resolve_path(base_dir, required_string(files_tbl, "hr"))
    win_path = resolve_path(base_dir, optional_string(files_tbl, "win"))
    structure_path = resolve_path(base_dir, optional_string(files_tbl, "structure"))
    out_path = resolve_path(base_dir, required_string(files_tbl, "out"))
    orbital_specs = _parse_input_orbital_specs(cfg)

    atol = haskey(render_tbl, "atol") ? required_float(render_tbl, "atol") : DEFAULT_ATOL
    include_second_quant = haskey(render_tbl, "include_second_quant") ? required_bool(render_tbl, "include_second_quant") : DEFAULT_INCLUDE_SECOND_QUANT
    selection = _parse_selection(cfg, files_tbl, render_tbl, base_dir)

    return LatexConfig(
        input_path,
        hr_path,
        win_path,
        structure_path,
        out_path,
        orbital_specs,
        atol,
        include_second_quant,
        selection,
    )
end

end
