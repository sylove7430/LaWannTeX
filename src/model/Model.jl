module Model

using ..WannierTypes: HrBlocks, RKey

export SelectionMode, LegacySelection, OrbitalSelection
export LatexConfig, OrbitalSpec, HoppingEntry, LatexRunResult
export orbital_labels, active_indices

abstract type SelectionMode end

struct LegacySelection <: SelectionMode
    max_rblocks::Union{Nothing, Int}
end

struct OrbitalSelection <: SelectionMode
    wannier_indices::Vector{Int}
    max_nn::Int
    distance_tol::Float64
end

struct OrbitalSpec
    index::Int
    label::String
    position_frac::Union{Nothing, NTuple{3, Float64}}
end

struct LatexConfig
    input_path::String
    hr_path::String
    win_path::Union{Nothing, String}
    structure_path::Union{Nothing, String}
    out_path::String
    orbital_specs::Union{Nothing, Vector{OrbitalSpec}}
    atol::Float64
    include_second_quant::Bool
    selection::SelectionMode
end

struct HoppingEntry
    R::RKey
    m::Int
    n::Int
    t::ComplexF64
    distance::Union{Nothing, Float64}
    nn_shell::Union{Nothing, Int}
end

struct LatexRunResult
    config::LatexConfig
    hr::HrBlocks
    orbital_specs::Vector{OrbitalSpec}
    orbital_labels::Vector{String}
    label_source::String
    entries::Vector{HoppingEntry}
    active_indices::Vector{Int}
    output_path::String
end

orbital_labels(specs::Vector{OrbitalSpec}) = [spec.label for spec in specs]

function active_indices(cfg::LatexConfig, nw::Int)::Vector{Int}
    if cfg.selection isa LegacySelection
        return collect(1:nw)
    end
    return (cfg.selection::OrbitalSelection).wannier_indices
end

function Base.getproperty(cfg::LatexConfig, name::Symbol)
    if name === :orbital_labels
        specs = getfield(cfg, :orbital_specs)
        return isnothing(specs) ? nothing : orbital_labels(specs)
    elseif name === :max_rblocks
        selection = getfield(cfg, :selection)
        return selection isa LegacySelection ? selection.max_rblocks : nothing
    end
    return getfield(cfg, name)
end

end
