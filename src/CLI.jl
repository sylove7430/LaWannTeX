module LatexCLI

include(joinpath(@__DIR__, "Core.jl"))

const InputIO = LatexCore.InputIO
const Model = LatexCore.Model
const Service = LatexCore.Service
const Table = LatexCore.Table

export run_main, parse_args, print_usage

function print_usage()
    println("""
Usage:
  julia main.jl input.toml
  julia main.jl --input input.toml
""")
end

function parse_args(args::Vector{String})
    input_path = nothing
    show_help = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--input"
            i += 1
            i <= length(args) || throw(ArgumentError("Missing value for --input"))
            input_path = String(args[i])
        elseif arg in ("-h", "--help")
            show_help = true
        elseif startswith(arg, "-")
            throw(ArgumentError("Unknown argument: $arg"))
        elseif isnothing(input_path)
            input_path = String(arg)
        else
            throw(ArgumentError("Unexpected positional argument: $arg"))
        end
        i += 1
    end

    return (
        input_path = input_path,
        show_help = show_help,
    )
end

function _summary_limit(result::Model.LatexRunResult)::Int
    return min(length(result.entries), Service.entry_limit(result.config))
end

function _render_summary(result::Model.LatexRunResult, input_path::AbstractString)
    println("Input file: ", abspath(input_path))
    println("Source hr: ", abspath(result.config.hr_path))
    println("Basis size: ", result.hr.num_wann)
    println("R blocks: ", result.hr.nrpts)
    println("Label source: ", result.label_source)
    println("Mode: ", result.config.selection isa Model.LegacySelection ? "legacy" : "orbital")
    println("Active basis: ", join(result.active_indices, ", "))
    println("Rendered entries: ", min(_summary_limit(result), length(result.entries)), " / ", length(result.entries))
    println("Output tex: ", abspath(result.output_path))
    Table.print_table(result.entries; top=_summary_limit(result))
    return nothing
end

function run_main(args::Vector{String})::Int
    opts = parse_args(args)
    if opts.show_help
        print_usage()
        return 0
    end
    isnothing(opts.input_path) && throw(ArgumentError("Missing input.toml path"))

    cfg = InputIO.read_input(opts.input_path)
    result = Service.run(cfg)
    _render_summary(result, opts.input_path)
    return 0
end

end
