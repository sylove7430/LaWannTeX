#!/usr/bin/env julia

import Pkg

Pkg.activate(@__DIR__; io=devnull)

include(joinpath(@__DIR__, "src", "CLI.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(LatexCLI.run_main(ARGS))
    catch err
        if err isa ArgumentError
            println(stderr, err.msg)
            LatexCLI.print_usage()
            exit(1)
        end
        rethrow()
    end
end
