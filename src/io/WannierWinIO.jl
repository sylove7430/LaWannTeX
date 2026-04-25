module WannierWinIO

using ..LatexConstants: BOHR_TO_ANG, expand_orbital_symbol
using ..Model: OrbitalSpec

export read_win_lattice, read_win_orbital_specs

function _clean_line(line::AbstractString)::String
    return strip(replace(String(line), r"[!#].*$" => ""))
end

function _read_block(path::AbstractString, block_name::AbstractString)::Vector{String}
    lines = String[]
    begin_token = "begin $(lowercase(block_name))"
    end_token = "end $(lowercase(block_name))"
    in_block = false

    for raw in readlines(path)
        stripped = _clean_line(raw)
        isempty(stripped) && continue

        lower = lowercase(stripped)
        if lower == begin_token
            in_block = true
            continue
        elseif lower == end_token
            break
        end

        in_block || continue
        push!(lines, stripped)
    end

    isempty(lines) && error("Missing $block_name block in $(abspath(path))")
    return lines
end

function _optional_block(path::AbstractString, block_name::AbstractString)::Union{Nothing, Vector{String}}
    try
        return _read_block(path, block_name)
    catch
        return nothing
    end
end

function _parse_vec3_line(line::AbstractString, key::AbstractString)::Vector{Float64}
    parts = split(replace(line, "," => " "))
    length(parts) == 3 || error("$key must contain exactly three numbers")
    return [parse(Float64, part) for part in parts]
end

function _read_block_with_optional_unit(path::AbstractString, block_name::AbstractString)
    block = _read_block(path, block_name)
    start_idx = 1
    scale = 1.0

    if !isempty(block)
        unit = lowercase(replace(block[1], r"\s+" => ""))
        if unit in ("ang", "angstrom")
            start_idx = 2
        elseif unit in ("bohr", "au", "a.u.")
            start_idx = 2
            scale = BOHR_TO_ANG
        end
    end
    return block, start_idx, scale
end

function read_win_lattice(path::AbstractString)::Matrix{Float64}
    block, start_idx, scale = _read_block_with_optional_unit(path, "unit_cell_cart")
    length(block) - start_idx + 1 >= 3 || error("unit_cell_cart in $(abspath(path)) must contain three lattice vectors")

    lattice = Matrix{Float64}(undef, 3, 3)
    for col in 1:3
        lattice[:, col] = scale .* _parse_vec3_line(block[start_idx + col - 1], "unit_cell_cart line $col")
    end
    return lattice
end

function _read_atoms_block(path::AbstractString, block_name::AbstractString)::Union{Nothing, Vector{Tuple{String, NTuple{3, Float64}}}}
    block = _optional_block(path, block_name)
    isnothing(block) && return nothing

    unit = block_name == "atoms_cart" ? 1.0 : 1.0
    start_idx = 1
    if !isempty(block)
        first_line = lowercase(replace(block[1], r"\s+" => ""))
        if block_name == "atoms_cart"
            if first_line in ("ang", "angstrom")
                start_idx = 2
            elseif first_line in ("bohr", "au", "a.u.")
                start_idx = 2
                unit = BOHR_TO_ANG
            end
        end
    end

    atoms = Tuple{String, NTuple{3, Float64}}[]
    for (offset, line) in enumerate(block[start_idx:end])
        fields = split(replace(line, "," => " "))
        length(fields) == 4 || error("$block_name line $(offset) in $(abspath(path)) must contain species and three coordinates")
        species = fields[1]
        coords = (unit * parse(Float64, fields[2]), unit * parse(Float64, fields[3]), unit * parse(Float64, fields[4]))
        push!(atoms, (species, coords))
    end
    return atoms
end

function _cart_to_frac(lattice::Matrix{Float64}, cart::NTuple{3, Float64})::NTuple{3, Float64}
    frac = lattice \ collect(cart)
    return Tuple(frac)
end

function _build_site_position_lookup(path::AbstractString)::Dict{String, Vector{NTuple{3, Float64}}}
    lattice = read_win_lattice(path)
    atoms_frac = _read_atoms_block(path, "atoms_frac")
    atoms_cart = _read_atoms_block(path, "atoms_cart")
    !isnothing(atoms_frac) || !isnothing(atoms_cart) || error("No atoms_frac/atoms_cart block found in $(abspath(path))")

    lookup = Dict{String, Vector{NTuple{3, Float64}}}()

    if !isnothing(atoms_frac)
        for (species, frac) in something(atoms_frac)
            push!(get!(lookup, species, NTuple{3, Float64}[]), frac)
        end
        return lookup
    end

    for (species, cart) in something(atoms_cart)
        push!(get!(lookup, species, NTuple{3, Float64}[]), _cart_to_frac(lattice, cart))
    end
    return lookup
end

function _parse_projection_line(line::AbstractString)
    # f=x_coord,y_coord,z_coord:orbital:local axes
    # Atom:orbital
    parts = split(line, ":"; limit=3)
    length(parts) >= 2 || return nothing
    coord_part = strip(parts[1])
    orb_part = strip(parts[2])
    return coord_part, orb_part
end

function _projection_orbitals(orb_part::AbstractString)::Vector{String}
    normalized = strip(orb_part)
    if occursin(r"(?i)\bl\s*=", normalized)
        match_obj = match(r"(?i)l\s*=\s*([a-zA-Z0-9_,\-\(\)]+)", normalized)
        isnothing(match_obj) && return String[]
        normalized = something(match_obj.captures[1])
    end

    orbitals = String[]
    for raw_orb in split(normalized, ",")
        isempty(strip(raw_orb)) && continue
        append!(orbitals, expand_orbital_symbol(raw_orb))
    end
    return orbitals
end

function _projection_site_counts(lines::Vector{String})::Dict{String, Int}
    counts = Dict{String, Int}()
    for stripped in lines
        parsed = _parse_projection_line(stripped)
        isnothing(parsed) && continue
        coord_part, orb_part = parsed
        isempty(_projection_orbitals(orb_part)) && continue
        startswith(lowercase(coord_part), "f=") && continue
        site_name = strip(coord_part)
        counts[site_name] = get(counts, site_name, 0) + 1
    end
    return counts
end

function _site_instance_name(site_name::AbstractString, instance_idx::Int, total_instances::Int)::String
    total_instances <= 1 && return String(site_name)
    return String(site_name) * string(instance_idx)
end

function read_win_orbital_specs(path::String, nw::Int)::Vector{OrbitalSpec}
    lines = _read_block(path, "projections")
    projection_site_counts = _projection_site_counts(lines)
    site_positions = try
        _build_site_position_lookup(path)
    catch
        Dict{String, Vector{NTuple{3, Float64}}}()
    end
    site_counters = Dict{String, Int}()

    specs = OrbitalSpec[]
    next_index = 1
    site_idx = 1

    for stripped in lines
        parsed = _parse_projection_line(stripped)
        isnothing(parsed) && continue
        coord_part, orb_part = parsed
        orbitals = _projection_orbitals(orb_part)
        isempty(orbitals) && continue

        if startswith(lowercase(coord_part), "f=")
            coords = _parse_vec3_line(coord_part[3:end], "projection fractional coordinate")
            positions = [(coords[1], coords[2], coords[3])]
            site_names = ["P$(site_idx)"]
        else
            site_name = strip(coord_part)
            positions = get(site_positions, site_name, nothing)
            isnothing(positions) && error("Unable to resolve projection site \"$site_name\" from atoms_frac/atoms_cart in $(abspath(path))")
            total_instances = length(positions)
            projection_count = get(projection_site_counts, site_name, 0)
            projection_count > 0 || error("Projection site \"$site_name\" was not counted in $(abspath(path))")

            if projection_count == 1
                site_names = [_site_instance_name(site_name, idx, total_instances) for idx in 1:total_instances]
            else
                counter = get(site_counters, site_name, 0) + 1
                counter <= length(positions) || error("Projection site \"$site_name\" appears more times than atoms listed in $(abspath(path))")
                site_counters[site_name] = counter
                positions = [positions[counter]]
                site_names = [_site_instance_name(site_name, counter, total_instances)]
            end
        end

        for (site_name, position) in zip(site_names, positions)
            for orbital in orbitals
                push!(specs, OrbitalSpec(next_index, "$(site_name):$(orbital)", position))
                next_index += 1
            end
        end
        site_idx += 1
    end

    isempty(specs) && error("No projections found in $(abspath(path))")
    length(specs) == nw || error("Win file has $(length(specs)) projections but hr.dat has num_wann=$nw")
    return specs
end

end
