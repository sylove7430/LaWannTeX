module Service

using ..LatexConstants: DEFAULT_ENTRY_LIMIT
using ..LabelIO
using ..LatticeIO
using ..Model
using ..Model: LatexConfig, LatexRunResult, LegacySelection, active_indices, orbital_labels
using ..OrbitalSelection: build_orbital_entries
using ..Render
using ..Table
using ..WannierHrIO

export run, entry_limit

function entry_limit(cfg::LatexConfig)::Int
    if cfg.selection isa LegacySelection
        max_rblocks = (cfg.selection::LegacySelection).max_rblocks
        isnothing(max_rblocks) && return DEFAULT_ENTRY_LIMIT
        return max(0, something(max_rblocks))
    end
    return typemax(Int)
end

function run(cfg::LatexConfig)::LatexRunResult
    hr = WannierHrIO.read_hr(cfg.hr_path)
    specs, label_source = LabelIO.resolve_orbital_specs(cfg, hr.num_wann)
    lattice = LatticeIO.resolve_lattice(cfg)

    entries, active = if cfg.selection isa LegacySelection
        (Table.build_legacy_entries(hr; atol=cfg.atol), collect(1:hr.num_wann))
    else
        lattice === nothing && error("orbital mode requires lattice information")
        (
            build_orbital_entries(hr, specs, lattice, cfg.selection::Model.OrbitalSelection; atol=cfg.atol),
            active_indices(cfg, hr.num_wann),
        )
    end

    out_path = cfg.out_path
    source_label = Render.display_source_label(hr, cfg.hr_path)
    latex_src = Render.build_latex(
        hr,
        entries;
        top=entry_limit(cfg),
        orbital_labels=orbital_labels(specs),
        active_indices=active,
        include_second_quant=cfg.include_second_quant,
        selection=cfg.selection,
        source_label=source_label,
        lattice=lattice,
    )
    mkpath(dirname(abspath(out_path)))
    write(out_path, latex_src)
    return LatexRunResult(cfg, hr, specs, orbital_labels(specs), label_source, entries, active, out_path)
end

end
