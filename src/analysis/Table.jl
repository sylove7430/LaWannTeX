module Table

using Printf

using ..Model: HoppingEntry
using ..WannierTypes: HrBlocks

export build_legacy_entries, print_table

function build_legacy_entries(hr::HrBlocks; atol::Float64=0.0)::Vector{HoppingEntry}
    entries = HoppingEntry[]
    for (R, H) in hr.hoppings
        for m in 1:hr.num_wann, n in 1:hr.num_wann
            value = H[m, n]
            abs(value) <= atol && continue
            push!(entries, HoppingEntry(R, m, n, value, nothing, nothing))
        end
    end
    sort!(entries; by = entry -> (-abs(entry.t), entry.R[1], entry.R[2], entry.R[3], entry.m, entry.n))
    return entries
end

function print_table(entries::Vector{HoppingEntry}; top::Int=typemax(Int))
    n = min(top, length(entries))
    println()
    if all(isnothing(entry.nn_shell) for entry in entries)
        println("  Dominant hopping parameters (|t| descending)")
    else
        println("  Selected hopping parameters (nn shell, then |t|)")
    end
    println()
    has_distance = any(!isnothing(entry.distance) for entry in entries)
    header =
        " " * lpad("Ind", 4) * lpad("m", 4) * lpad("n", 4) *
        (has_distance ? lpad("NN", 5) * lpad("dist", 12) : "") *
        lpad("R1", 5) * lpad("R2", 5) * lpad("R3", 5) *
        rpad("   Re[t] (eV)", 18) *
        rpad("Im[t] (eV)", 16) *
        rpad("|t| (eV)", 12)
    println(header)
    println("  " * "-"^max(80, length(header)))
    for (idx, entry) in enumerate(entries[1:n])
        row =
            " " * lpad(idx, 4) *
            lpad(entry.m, 4) * lpad(entry.n, 4) *
            (has_distance ? lpad(something(entry.nn_shell, -1), 5) * lpad(isnothing(entry.distance) ? "-" : @sprintf("%.6f", something(entry.distance)), 12) : "") *
            lpad(entry.R[1], 5) * lpad(entry.R[2], 5) * lpad(entry.R[3], 5) *
            lpad(@sprintf("% .6f", real(entry.t)), 18) *
            lpad(@sprintf("% .6f", imag(entry.t)), 16) *
            lpad(@sprintf("%.6f", abs(entry.t)), 12)
        println(row)
    end
    println()
    return nothing
end

end
