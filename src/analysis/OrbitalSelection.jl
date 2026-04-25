module OrbitalSelection

using LinearAlgebra

using ..Model
using ..Model: HoppingEntry, OrbitalSpec
using ..Table: build_legacy_entries
using ..WannierTypes: HrBlocks, RKey

export build_legacy_entries, build_orbital_entries

function _validate_selected_indices(selection::Model.OrbitalSelection, nw::Int)
    for idx in selection.wannier_indices
        1 <= idx <= nw || error("selection.wannier_indices contains out-of-range index $idx for num_wann=$nw")
    end
    return nothing
end

function _spec_position(spec::OrbitalSpec)::NTuple{3, Float64}
    isnothing(spec.position_frac) && error("Missing fractional position for orbital $(spec.label)")
    return something(spec.position_frac)
end

function _raw_cartesian_delta(
    R::RKey,
    m::Int,
    n::Int,
    specs::Vector{OrbitalSpec},
    lattice::Matrix{Float64},
)::Vector{Float64}
    tau_m = collect(_spec_position(specs[m]))
    tau_n = collect(_spec_position(specs[n]))
    frac_delta = collect(R) .+ tau_n .- tau_m
    return lattice * frac_delta
end

function _raw_distance(
    R::RKey,
    m::Int,
    n::Int,
    specs::Vector{OrbitalSpec},
    lattice::Matrix{Float64},
)::Float64
    return norm(_raw_cartesian_delta(R, m, n, specs, lattice))
end

function _is_strict_onsite(entry::HoppingEntry, tol::Float64)
    return entry.R == (0, 0, 0) && !isnothing(entry.distance) && something(entry.distance) <= tol
end

function _assign_shells!(entries::Vector{HoppingEntry}; tol::Float64=1e-5)
    radii = Float64[]
    for (idx, entry) in enumerate(entries)
        if _is_strict_onsite(entry, tol)
            entries[idx] = HoppingEntry(entry.R, entry.m, entry.n, entry.t, entry.distance, 0)
            continue
        end
        dist = something(entry.distance)
        dist > tol || error("Encountered non-home-cell zero-distance bond at (R=$(entry.R), m=$(entry.m), n=$(entry.n))")
        if !any(abs(dist - radius) <= tol for radius in radii)
            push!(radii, dist)
        end
    end

    sort!(radii)
    for (idx, entry) in enumerate(entries)
        isnothing(entry.nn_shell) || continue
        dist = something(entry.distance)
        shell = findfirst(radius -> abs(dist - radius) <= tol, radii)
        isnothing(shell) && error("Failed to assign shell for distance $dist")
        entries[idx] = HoppingEntry(entry.R, entry.m, entry.n, entry.t, entry.distance, shell)
    end

    return entries
end

function _conjugate_R(R::RKey)::RKey
    return (-R[1], -R[2], -R[3])
end

function _ensure_hermitian_complete(entries::Vector{HoppingEntry}, hr::HrBlocks)::Vector{HoppingEntry}
    lookup = Dict((entry.R, entry.m, entry.n) => entry for entry in entries)
    completed = copy(entries)
    for entry in entries
        conj_key = (_conjugate_R(entry.R), entry.n, entry.m)
        haskey(lookup, conj_key) && continue

        haskey(hr.hoppings, conj_key[1]) || error("Missing Hermitian conjugate R block for entry ($(entry.R), $(entry.m), $(entry.n))")
        conj_matrix = hr.hoppings[conj_key[1]]
        conj_value = conj_matrix[conj_key[2], conj_key[3]]
        if abs(conj_value - conj(entry.t)) > 1e-8
            error("Hermitian conjugate mismatch for entry ($(entry.R), $(entry.m), $(entry.n))")
        end

        conj_entry = HoppingEntry(conj_key[1], conj_key[2], conj_key[3], conj_value, entry.distance, entry.nn_shell)
        push!(completed, conj_entry)
        lookup[conj_key] = conj_entry
    end
    return completed
end

function _sort_entries!(entries::Vector{HoppingEntry})
    sort!(
        entries;
        by = e -> (
            something(e.nn_shell, typemax(Int)),
            -abs(e.t),
            e.m,
            e.n,
            e.R[1], e.R[2], e.R[3],
        ),
    )
    return entries
end

function build_orbital_entries(
    hr::HrBlocks,
    specs::Vector{OrbitalSpec},
    lattice::Matrix{Float64},
    selection::Model.OrbitalSelection;
    atol::Float64=0.0,
)::Vector{HoppingEntry}
    _validate_selected_indices(selection, hr.num_wann)
    selected = Set(selection.wannier_indices)

    entries = HoppingEntry[]
    for (R, H) in hr.hoppings
        for m in selection.wannier_indices, n in selection.wannier_indices
            m in selected || continue
            n in selected || continue
            value = H[m, n]
            abs(value) <= atol && continue
            distance = _raw_distance(R, m, n, specs, lattice)
            push!(entries, HoppingEntry(R, m, n, value, distance, nothing))
        end
    end

    entries = filter(entry -> abs(entry.t) > atol, entries)
    _assign_shells!(entries; tol=selection.distance_tol)
    entries = filter(entry -> begin
        shell = something(entry.nn_shell)
        shell == 0 || shell <= selection.max_nn
    end, entries)
    entries = _ensure_hermitian_complete(entries, hr)
    _sort_entries!(entries)
    return entries
end

end
