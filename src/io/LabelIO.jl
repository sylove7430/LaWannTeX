module LabelIO

using ..Model: LatexConfig, LegacySelection, OrbitalSelection, OrbitalSpec
using ..WannierWinIO: read_win_orbital_specs

export resolve_orbital_specs

function _default_orbital_specs(nw::Int)::Vector{OrbitalSpec}
    specs = OrbitalSpec[]
    for i in 1:nw
        idx = i - 1
        label = ""
        while true
            label = string(Char('A' + (idx % 26)), label)
            idx = fld(idx, 26) - 1
            idx < 0 && break
        end
        push!(specs, OrbitalSpec(i, label, nothing))
    end
    return specs
end

function _validate_positions!(cfg::LatexConfig, specs::Vector{OrbitalSpec})
    cfg.selection isa OrbitalSelection || return specs
    for spec in specs
        isnothing(spec.position_frac) && error("orbital mode requires fractional positions for every orbital; missing position for $(spec.label)")
    end
    return specs
end

function resolve_orbital_specs(cfg::LatexConfig, nw::Int)::Tuple{Vector{OrbitalSpec}, String}
    if !isnothing(cfg.orbital_specs)
        length(cfg.orbital_specs) == nw || error("[[orbitals]] expands to $(length(cfg.orbital_specs)) orbitals but hr.dat has num_wann=$nw")
        return _validate_positions!(cfg, cfg.orbital_specs), "from [[orbitals]]"
    end
    if !isnothing(cfg.win_path)
        specs = read_win_orbital_specs(cfg.win_path, nw)
        return _validate_positions!(cfg, specs), "from Wannier90 projections $(cfg.win_path)"
    end
    specs = _default_orbital_specs(nw)
    cfg.selection isa LegacySelection || error("orbital mode requires either [[orbitals]] with position_frac or files.win with projections/atoms metadata")
    return specs, "default labels"
end

end
