module Render

using Printf

using ..LatexConstants: ALIGN_ROW_END, ORBITAL_LATEX
using ..Model: HoppingEntry, LegacySelection, OrbitalSelection, SelectionMode
using ..WannierTypes: HrBlocks, RKey

export build_latex, display_source_label

function _fmt_complex(z::ComplexF64; digits::Int=6)::String
    re = round(real(z); digits=digits)
    im = round(imag(z); digits=digits)
    tol = 10.0^(-digits)
    if abs(im) < tol
        return string(re)
    elseif abs(re) < tol
        return "$(im)i"
    end
    sign = im >= 0 ? "+" : "-"
    return "$(re) $(sign) $(abs(im))i"
end

function _R_to_latex(R::RKey)::String
    parts = String[]
    labels = ("\\mathbf{a}_1", "\\mathbf{a}_2", "\\mathbf{a}_3")
    for (idx, ri) in enumerate(R)
        ri == 0 && continue
        abs_coeff = abs(ri) == 1 ? "" : string(abs(ri))
        term = abs_coeff * labels[idx]
        if isempty(parts)
            push!(parts, ri < 0 ? "-" * term : term)
        else
            push!(parts, ri < 0 ? " - " * term : " + " * term)
        end
    end
    return isempty(parts) ? "\\mathbf{0}" : join(parts)
end

function _looks_like_timestamp_header(text::AbstractString)::Bool
    normalized = lowercase(strip(text))
    startswith(normalized, "written on ") && return true
    startswith(normalized, "supercell(") && occursin("written on ", normalized) && return true
    return false
end

function display_source_label(hr::HrBlocks, hr_path::AbstractString)::String
    header = strip(hr.header)
    if isempty(header) || _looks_like_timestamp_header(header)
        return basename(hr_path)
    end
    return header
end

function _latex_text_escape(text::AbstractString)::String
    buf = IOBuffer()
    for ch in text
        if ch == '_'
            print(buf, raw"\_")
        elseif ch == '&'
            print(buf, raw"\&")
        elseif ch == '%'
            print(buf, raw"\%")
        elseif ch == '#'
            print(buf, raw"\#")
        elseif ch == '$'
            print(buf, raw"\$")
        elseif ch == '{'
            print(buf, raw"\{")
        elseif ch == '}'
            print(buf, raw"\}")
        else
            print(buf, ch)
        end
    end
    return String(take!(buf))
end

function _format_orbital(orb::String)::String
    return get(ORBITAL_LATEX, lowercase(orb), "\\text{" * _latex_text_escape(orb) * "}")
end

function _format_spin(spin::AbstractString)::String
    normalized = lowercase(strip(String(spin)))
    normalized in ("", "none", "spinless") && return ""
    normalized in ("up", "u", "spin-up", "spinup", "↑") && return raw"\uparrow"
    normalized in ("down", "dn", "d", "spin-down", "spindown", "↓") && return raw"\downarrow"
    return "\\text{:" * _latex_text_escape(String(spin)) * "}"
end

function _latex_subscript(label::String)::String
    parts = split(label, ":")
    if length(parts) == 2
        site, orb = parts
        return "\\text{" * _latex_text_escape(String(site)) * ":}" * _format_orbital(String(orb))
    elseif length(parts) == 3
        site, orb, spin = parts
        return "\\text{" * _latex_text_escape(String(site)) * ":}" * _format_orbital(String(orb)) * _format_spin(String(spin))
    end
    all(ch -> isletter(ch) || isdigit(ch), label) && return label
    return "\\text{" * _latex_text_escape(label) * "}"
end

function _latex_label_list(labels::Vector{String}, indices::Vector{Int})::String
    return join(["\$" * _latex_subscript(labels[idx]) * "\$" for idx in indices], ", ")
end

function _fmt_real(x::Real; digits::Int=6)::String
    rounded = round(Float64(x); digits=digits)
    abs(rounded) < 10.0^(-digits) && return "0"
    return @sprintf("%.6g", rounded)
end

function _vector_to_latex(vec::AbstractVector{<:Real})::String
    return "(" * join([_fmt_real(item) for item in vec], ", ") * ")"
end

function _build_lattice_latex!(buf::IOBuffer, lattice::Matrix{Float64})
    println(
        buf,
        raw"The direct-lattice basis vectors used in $\mathbf{R} = R_1\mathbf{a}_1 + R_2\mathbf{a}_2 + R_3\mathbf{a}_3$ are",
    )
    println(buf)
    println(buf, raw"\begin{align}")
    for idx in 1:3
        if idx == 3
            println(buf, "  \\mathbf{a}_$(idx) &= " * _vector_to_latex(lattice[:, idx]))
        else
            println(buf, "  \\mathbf{a}_$(idx) &= " * _vector_to_latex(lattice[:, idx]) * ALIGN_ROW_END)
        end
    end
    println(buf, raw"\end{align}")
    println(buf)
    return nothing
end

function _entry_dict(entries::Vector{HoppingEntry})
    mn_dict = Dict{Tuple{Int, Int}, Vector{Tuple{RKey, ComplexF64}}}()
    for entry in entries
        push!(get!(mn_dict, (entry.m, entry.n), Tuple{RKey, ComplexF64}[]), (entry.R, entry.t))
    end
    return mn_dict
end

function _wrap_signed_terms(lhs::AbstractString, signed_terms::Vector{String}; terms_per_line::Int=3, final_suffix::AbstractString="")
    isempty(signed_terms) && return [String(lhs) * " &= 0" * String(final_suffix)]

    lines = String[]
    total_chunks = cld(length(signed_terms), terms_per_line)

    for chunk_idx in 1:total_chunks
        start_idx = (chunk_idx - 1) * terms_per_line + 1
        end_idx = min(chunk_idx * terms_per_line, length(signed_terms))
        chunk = signed_terms[start_idx:end_idx]

        prefix = chunk_idx == 1 ? String(lhs) * " &= " : raw"  &\quad "
        body = if chunk_idx == 1
            first_term = chunk[1]
            trailing_terms = chunk[2:end]
            isempty(trailing_terms) ? first_term : first_term * "\n    " * join(trailing_terms, "\n    ")
        else
            join(chunk, "\n    ")
        end

        suffix = chunk_idx == total_chunks ? String(final_suffix) : ALIGN_ROW_END
        push!(lines, prefix * body * suffix)
    end

    return lines
end

function _is_conjugate_pair(mn_dict, m::Int, n::Int)
    haskey(mn_dict, (m, n)) || return false
    haskey(mn_dict, (n, m)) || return false
    terms_mn = mn_dict[(m, n)]
    terms_nm = mn_dict[(n, m)]
    length(terms_mn) == length(terms_nm) || return false
    lookup = Dict(R => t for (R, t) in terms_mn)
    for (R, t_nm) in terms_nm
        t_mn = get(lookup, (-R[1], -R[2], -R[3]), nothing)
        isnothing(t_mn) && return false
        abs(t_nm - conj(t_mn)) > 1e-8 && return false
    end
    return true
end

function _build_second_quant_latex!(
    buf::IOBuffer,
    hr::HrBlocks,
    entries::Vector{HoppingEntry};
    orbital_labels::Vector{String},
    active_indices::Vector{Int},
)
    println(buf, "The Hamiltonian in second-quantization reads")
    println(buf)
    println(buf, raw"\begin{equation}")
    println(buf, raw"  H = \sum_{\mathbf{k}} \sum_{m,n} H_{mn}(\mathbf{k})\, c^\dagger_{m\mathbf{k}} c_{n\mathbf{k}},")
    println(buf, raw"\end{equation}")
    println(buf)
    println(buf, raw"where $H_{mn}(\mathbf{k}) = \sum_{\mathbf{R}} H_{mn}(\mathbf{R})\, e^{i\mathbf{k}\cdot\mathbf{R}}$.")
    println(buf)

    mn_dict = _entry_dict(entries)
    H0 = get(hr.hoppings, (0, 0, 0), zeros(ComplexF64, hr.num_wann, hr.num_wann))
    shown = Set{Tuple{Int, Int}}()

    println(buf, "Expanded term by term (orbital labels: " * _latex_label_list(orbital_labels, active_indices) * "):")
    println(buf)
    println(buf, raw"\begin{align}")
    print(buf, raw"H &= \sum_{\mathbf{k}} \Bigl[")
    println(buf)

    term_strs = String[]
    for m in active_indices
        label = _latex_subscript(orbital_labels[m])
        push!(term_strs, "  \\varepsilon_{$(label)}\\, c^\\dagger_{$(label),\\mathbf{k}} c_{$(label),\\mathbf{k}}")
        push!(shown, (m, m))
    end

    for m in active_indices, n in active_indices
        m == n && continue
        (m, n) in shown && continue
        (n, m) in shown && continue

        label_m = _latex_subscript(orbital_labels[m])
        label_n = _latex_subscript(orbital_labels[n])

        if m < n && _is_conjugate_pair(mn_dict, m, n)
            push!(term_strs, "  f_{$(label_m) \\to $(label_n)}(\\mathbf{k})\\, c^\\dagger_{$(label_m),\\mathbf{k}} c_{$(label_n),\\mathbf{k}}")
            push!(term_strs, "  f_{$(label_m) \\to $(label_n)}^*(\\mathbf{k})\\, c^\\dagger_{$(label_n),\\mathbf{k}} c_{$(label_m),\\mathbf{k}}")
            push!(shown, (m, n))
            push!(shown, (n, m))
        elseif haskey(mn_dict, (m, n))
            push!(term_strs, "  f_{$(label_m) \\to $(label_n)}(\\mathbf{k})\\, c^\\dagger_{$(label_m),\\mathbf{k}} c_{$(label_n),\\mathbf{k}}")
            push!(shown, (m, n))
        end
    end

    for (idx, term) in enumerate(term_strs)
        is_last = idx == length(term_strs)
        line_break = (!is_last && idx % 2 == 0) ? " \\\\" : ""
        prefix = if idx == 1
            "    "
        elseif (idx - 1) % 2 == 0
            raw"  &\quad + "
        else
            "    + "
        end
        println(buf, prefix * term * line_break)
    end
    println(buf, raw"  \Bigr],")
    println(buf, raw"\end{align}")
    println(buf)

    println(buf, raw"where the on-site energies are")
    println(buf)
    println(buf, raw"\begin{align}")
    for (idx, m) in enumerate(active_indices)
        label = _latex_subscript(orbital_labels[m])
        eps = real(H0[m, m])
        eps_str = @sprintf("%.6g", eps)
        suffix = idx == length(active_indices) ? "" : ALIGN_ROW_END
        println(buf, "  \\varepsilon_{$(label)} &= \\operatorname{Re}[H_{$(label),$(label)}(\\mathbf{0})] = $(eps_str)\\text{ eV}" * suffix)
    end
    println(buf, raw"\end{align}")
    println(buf)

    def_pairs = Tuple{Int, Int}[]
    seen = Set{Tuple{Int, Int}}()
    for m in active_indices, n in active_indices
        m == n && continue
        canon = m < n ? (m, n) : (n, m)
        canon in seen && continue
        if (m, n) in shown || (n, m) in shown
            push!(def_pairs, canon)
            push!(seen, canon)
        end
    end

    isempty(def_pairs) && return nothing

    println(buf, raw"The structure functions are defined as")
    println(buf)
    println(buf, raw"\begin{align}")
    for (pair_idx, (m, n)) in enumerate(def_pairs)
        label_m = _latex_subscript(orbital_labels[m])
        label_n = _latex_subscript(orbital_labels[n])
        terms = get(mn_dict, (m, n), Tuple{RKey, ComplexF64}[])
        isempty(terms) && continue

        lhs = "  f_{$(label_m) \\to $(label_n)}(\\mathbf{k})"
        signed_terms = String[]
        for (idx, (R, t)) in enumerate(terms)
            re_t = round(real(t); digits=6)
            im_t = round(imag(t); digits=6)
            tol = 1e-9

            if abs(im_t) < tol
                coeff_str = abs(abs(re_t) - 1.0) < tol ? (re_t < 0 ? "-" : "") : @sprintf("%.6g", re_t) * "\\,"
            elseif abs(re_t) < tol
                coeff_str = @sprintf("%.6g", im_t) * "i\\,"
            else
                coeff_str = "(" * @sprintf("%.6g", re_t) * (im_t >= 0 ? "+" : "") * @sprintf("%.6g", im_t) * "i)\\,"
            end

            phase = R == (0, 0, 0) ? "" : "e^{i\\mathbf{k}\\cdot($(_R_to_latex(R)))}"
            term = coeff_str * phase
            if idx == 1
                push!(signed_terms, term)
            elseif startswith(coeff_str, "-")
                push!(signed_terms, "- " * term[2:end])
            else
                push!(signed_terms, "+ " * term)
            end
        end
        suffix = pair_idx == length(def_pairs) ? "" : ALIGN_ROW_END
        for line in _wrap_signed_terms(lhs, signed_terms; terms_per_line=3, final_suffix=suffix)
            println(buf, line)
        end
    end
    println(buf, raw"\end{align}")
    println(buf)
    return nothing
end

function _build_table!(
    buf::IOBuffer,
    entries::Vector{HoppingEntry},
    orbital_labels::Vector{String},
    display_label::AbstractString,
    selection::SelectionMode,
)
    is_orbital = selection isa OrbitalSelection
    if is_orbital
        println(buf, raw"  \begin{tabular}{llrrrrr S[table-format=+1.6] S[table-format=+1.6] S[table-format=1.6]}")
        println(buf, raw"    \toprule")
        println(buf, "    from & to & NN & distance (\$\\AA\$) & \$R_1\$ & \$R_2\$ & \$R_3\$ & {\$\\operatorname{Re}[t]\$ (eV)} & {\$\\operatorname{Im}[t]\$ (eV)} & {\$|t|\$ (eV)} \\\\")
    else
        println(buf, raw"  \begin{tabular}{llrrr S[table-format=+1.6] S[table-format=+1.6] S[table-format=1.6]}")
        println(buf, raw"    \toprule")
        println(buf, "    from & to & \$R_1\$ & \$R_2\$ & \$R_3\$ & {\$\\operatorname{Re}[t]\$ (eV)} & {\$\\operatorname{Im}[t]\$ (eV)} & {\$|t|\$ (eV)} \\\\")
    end
    println(buf, raw"    \midrule")
    for entry in entries
        common =
            "    \$" * _latex_subscript(orbital_labels[entry.m]) * "\$" *
            " & \$" * _latex_subscript(orbital_labels[entry.n]) * "\$"
        if is_orbital
            println(
                buf,
                common *
                " & $(something(entry.nn_shell)) & " * @sprintf("%.6f", something(entry.distance)) *
                " & $(entry.R[1]) & $(entry.R[2]) & $(entry.R[3]) & " *
                @sprintf("% .6f", real(entry.t)) * " & " *
                @sprintf("% .6f", imag(entry.t)) * " & " *
                @sprintf("%.6f", abs(entry.t)) * ALIGN_ROW_END,
            )
        else
            println(
                buf,
                common *
                " & $(entry.R[1]) & $(entry.R[2]) & $(entry.R[3]) & " *
                @sprintf("% .6f", real(entry.t)) * " & " *
                @sprintf("% .6f", imag(entry.t)) * " & " *
                @sprintf("%.6f", abs(entry.t)) * ALIGN_ROW_END,
            )
        end
    end
    println(buf, raw"    \bottomrule")
    println(buf, raw"  \end{tabular}")
    println(buf, raw"\end{table}")
    println(buf)
    return nothing
end

function build_latex(
    hr::HrBlocks,
    entries::Vector{HoppingEntry};
    top::Int=typemax(Int),
    orbital_labels::Vector{String},
    active_indices::Vector{Int},
    include_second_quant::Bool=true,
    selection::SelectionMode=LegacySelection(nothing),
    source_label::AbstractString=strip(hr.header),
    lattice::Union{Nothing, Matrix{Float64}}=nothing,
)
    nw = hr.num_wann
    render_explicit_matrix_form = false
    length(orbital_labels) == nw || error("orbital_labels length $(length(orbital_labels)) != num_wann $nw")
    all(1 <= idx <= nw for idx in active_indices) || error("active_indices must lie within 1:num_wann")
    display_label = strip(String(source_label))

    buf = IOBuffer()
    println(buf, raw"% Auto-generated by latex from " * display_label)
    println(buf, raw"\documentclass[12pt]{article}")
    println(buf, raw"\usepackage{amsmath, booktabs, siunitx, geometry}")
    println(buf, raw"\geometry{margin=2.5cm}")
    println(buf, raw"\begin{document}")
    println(buf)

    println(buf, raw"\section*{Tight-Binding Hamiltonian}")
    println(buf)
    println(buf, raw"The real-space matrix elements define")
    println(buf)
    println(buf, raw"\begin{equation}")
    println(buf, raw"  H(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}\, H(\mathbf{R}),")
    println(buf, raw"\end{equation}")
    println(buf)
    println(buf, raw"where $\mathbf{k}$ is given in reduced coordinates and $H_{mn}(\mathbf{R}) = \langle \mathbf{0}m | \hat{H} | \mathbf{R}n \rangle$.")
    println(buf)
    !isnothing(lattice) && _build_lattice_latex!(buf, lattice)

    shown_entries = entries[1:min(top, length(entries))]
    if !isempty(shown_entries)
        println(buf, raw"\medskip")
        if selection isa OrbitalSelection
            println(buf, raw"The selected contributions up to the requested neighbor shell are")
        else
            println(buf, raw"The dominant contributions (" * string(length(shown_entries)) * raw" largest $|t_{mn}|$) are")
        end
        println(buf)
        println(buf, raw"\begin{equation}")
        println(buf, raw"  H(\mathbf{k}) \approx \sum_{\text{listed}} t_{mn}(\mathbf{R})\, e^{i\mathbf{k}\cdot\mathbf{R}}\, |m\rangle\langle n| + h.c.,")
        println(buf, raw"\end{equation}")
        println(buf)
        println(buf, raw"with $t_{mn}(\mathbf{R}) = H_{mn}(\mathbf{R})$ in eV:")
        println(buf)
        println(buf, raw"\begin{table}[h]")
        println(buf, raw"  \centering")
        println(buf, raw"  \caption{Hopping parameters extracted from \texttt{" * _latex_text_escape(display_label) * raw"}.}")
        println(buf, raw"  \label{tab:hoppings}")
        _build_table!(buf, shown_entries, orbital_labels, display_label, selection)

        if render_explicit_matrix_form
            println(buf, raw"\subsection*{Explicit matrix form}")
            println(buf)
            println(buf, raw"Collecting the $R$-resolved contributions and writing $\phi_{\mathbf{R}} \equiv e^{i\mathbf{k}\cdot\mathbf{R}}$:")
            println(buf)

            mn_dict = _entry_dict(shown_entries)
            println(buf, raw"\begin{align}")
            active_set = Set(active_indices)
            eq_count = 0
            total_eq = count(((m, n),) -> (m in active_set && n in active_set), collect(keys(mn_dict)))
            for m in active_indices, n in active_indices
                terms = get(mn_dict, (m, n), nothing)
                isnothing(terms) && continue
                eq_count += 1
                lhs = "  H_{$(_latex_subscript(orbital_labels[m])),$(_latex_subscript(orbital_labels[n]))}(\\mathbf{k})"
                signed_terms = String[]
                for (idx, (R, t)) in enumerate(terms)
                    phase = R == (0, 0, 0) ? "" : "e^{i\\mathbf{k}\\cdot($(_R_to_latex(R)))}"
                    term = isempty(phase) ? _fmt_complex(t) : _fmt_complex(t) * "\\," * phase
                    if idx == 1
                        push!(signed_terms, term)
                    elseif startswith(term, "-")
                        push!(signed_terms, "- " * term[2:end])
                    else
                        push!(signed_terms, "+ " * term)
                    end
                end
                suffix = eq_count == total_eq ? "" : ALIGN_ROW_END
                for line in _wrap_signed_terms(lhs, signed_terms; terms_per_line=3, final_suffix=suffix)
                    println(buf, line)
                end
            end
            println(buf, raw"\end{align}")
            println(buf)
        end

        if include_second_quant
            println(buf, raw"\subsection*{Second-quantization form}")
            println(buf)
            _build_second_quant_latex!(buf, hr, shown_entries; orbital_labels=orbital_labels, active_indices=active_indices)
        end
    end

    println(buf, raw"\end{document}")
    return String(take!(buf))
end

end
