module InputParsing

export required_table
export required_string, optional_string, required_string_vector
export required_float, optional_int, required_bool
export resolve_path
export parse_vec3_float

_prefix(context::AbstractString) = isempty(context) ? "" : string(context, ": ")

function required_table(tbl, key::AbstractString; context::AbstractString="")
    value = get(tbl, key, nothing)
    value isa AbstractDict || error(_prefix(context) * "Missing [$key] table")
    return value
end

function required_string(tbl, key::AbstractString; context::AbstractString="")
    haskey(tbl, key) || error(_prefix(context) * "Missing required key \"$key\"")
    value = tbl[key]
    value isa AbstractString || error(_prefix(context) * "Key \"$key\" must be a string")
    text = strip(String(value))
    isempty(text) && error(_prefix(context) * "Key \"$key\" cannot be empty")
    return text
end

function optional_string(tbl, key::AbstractString; default=nothing, context::AbstractString="")
    haskey(tbl, key) || return default
    value = tbl[key]
    value isa AbstractString || error(_prefix(context) * "Key \"$key\" must be a string")
    text = strip(String(value))
    return isempty(text) ? default : text
end

function required_string_vector(tbl, key::AbstractString; context::AbstractString="")
    haskey(tbl, key) || error(_prefix(context) * "Missing required key \"$key\"")
    raw = tbl[key]
    raw isa AbstractVector || error(_prefix(context) * "Key \"$key\" must be an array of strings")
    values = String[]
    for (i, item) in enumerate(raw)
        item isa AbstractString || error(_prefix(context) * "Key \"$key\" must contain only strings (bad entry at index $i)")
        text = strip(String(item))
        isempty(text) && error(_prefix(context) * "Key \"$key\" cannot contain empty strings")
        push!(values, text)
    end
    isempty(values) && error(_prefix(context) * "Key \"$key\" cannot be an empty array")
    return values
end

function required_float(tbl, key::AbstractString; context::AbstractString="")
    haskey(tbl, key) || error(_prefix(context) * "Missing required key \"$key\"")
    value = tbl[key]
    value isa Real || error(_prefix(context) * "Key \"$key\" must be numeric")
    return Float64(value)
end

function optional_int(tbl, key::AbstractString; context::AbstractString="")
    haskey(tbl, key) || return nothing
    value = tbl[key]
    value isa Integer || error(_prefix(context) * "Key \"$key\" must be an integer")
    return Int(value)
end

function required_bool(tbl, key::AbstractString; context::AbstractString="")
    haskey(tbl, key) || error(_prefix(context) * "Missing required key \"$key\"")
    value = tbl[key]
    value isa Bool || error(_prefix(context) * "Key \"$key\" must be a boolean")
    return value
end

function resolve_path(base_dir::AbstractString, value::Union{Nothing, AbstractString}; empty_value=nothing)
    isnothing(value) && return empty_value
    text = strip(String(value))
    isempty(text) && return empty_value
    return isabspath(text) ? text : normpath(joinpath(base_dir, text))
end

function parse_vec3_float(raw, key::AbstractString)::NTuple{3, Float64}
    raw isa AbstractVector || error("$key must be a 3-element array")
    length(raw) == 3 || error("$key must have length 3")
    return (Float64(raw[1]), Float64(raw[2]), Float64(raw[3]))
end

end
