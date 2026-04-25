module LatexConstants

export DEFAULT_ATOL, DEFAULT_INCLUDE_SECOND_QUANT, DEFAULT_ENTRY_LIMIT, BOHR_TO_ANG
export ORBITAL_EXPANSIONS, ALIGN_ROW_END, ORBITAL_LATEX
export expand_orbital_symbol

const DEFAULT_ATOL = 1e-6
const DEFAULT_INCLUDE_SECOND_QUANT = true
const DEFAULT_ENTRY_LIMIT = 40
const BOHR_TO_ANG = 0.529177210903

const ORBITAL_EXPANSIONS = Dict(
    "s" => ("s",),
    "p" => ("pz", "px", "py"),
    "d" => ("dz2", "dxz", "dyz", "dx2-y2", "dxy"),
    # Real f harmonics following the Wannier90 projection naming convention.
    "f" => ("fz3", "fxz2", "fyz2", "fz(x2-y2)", "fxyz", "fx(x2-3y2)", "fy(3x2-y2)"),
)

const ALIGN_ROW_END = " \\\\"

const ORBITAL_LATEX = Dict(
    "dz2"    => "d_{z^2}",
    "dxz"    => "d_{xz}",
    "dyz"    => "d_{yz}",
    "dx2-y2" => "d_{x^2-y^2}",
    "dxy"    => "d_{xy}",
    "pz"     => "p_z",
    "px"     => "p_x",
    "py"     => "p_y",
    "s"      => "s",
    "fz3" => "f_{z^3}",
    "fxz2" => "f_{xz^2}",
    "fyz2" => "f_{yz^2}",
    "fz(x2-y2)" => "f_{z(x^2-y^2)}",
    "fxyz" => "f_{xyz}",
    "fx(x2-3y2)" => "f_{x(x^2-3y^2)}",
    "fy(3x2-y2)" => "f_{y(3x^2-y^2)}",
)

function expand_orbital_symbol(orb::AbstractString)::Vector{String}
    normalized = lowercase(strip(String(orb)))
    return collect(get(ORBITAL_EXPANSIONS, normalized, (normalized,)))
end

end
