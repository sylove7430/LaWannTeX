module WannierTypes

export RKey, HrBlocks

const RKey = NTuple{3, Int}

struct HrBlocks
    header::String
    num_wann::Int
    nrpts::Int
    hoppings::Dict{RKey, Matrix{ComplexF64}}

    function HrBlocks(
        header::AbstractString,
        num_wann::Integer,
        nrpts::Integer,
        hoppings::Dict{RKey, Matrix{ComplexF64}},
    )
        return new(String(header), Int(num_wann), Int(nrpts), hoppings)
    end
end

end
