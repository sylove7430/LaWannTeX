module CrystalCellModel

export CrystalCell
export wrap_fractional

struct CrystalCell
    comment::String
    lattice::Matrix{Float64}
    frac_positions::Matrix{Float64}
    species_names::Vector{String}
    species_ids::Vector{Int}
    source::String
end

function wrap_fractional(x::AbstractVector{<:Real}; tol::Float64=1e-10)
    y = Vector{Float64}(undef, length(x))
    y .= mod.(Float64.(x), 1.0)
    y[isapprox.(y, 1.0; atol=tol)] .= 0.0
    y[isapprox.(y, 0.0; atol=tol)] .= 0.0
    return y
end

function wrap_fractional(xs::AbstractMatrix{<:Real}; tol::Float64=1e-10)
    ys = Matrix{Float64}(undef, size(xs)...)
    ys .= mod.(Float64.(xs), 1.0)
    ys[isapprox.(ys, 1.0; atol=tol)] .= 0.0
    ys[isapprox.(ys, 0.0; atol=tol)] .= 0.0
    return ys
end

end
