import Pkg

Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

using Test

include(joinpath(@__DIR__, "..", "src", "CLI.jl"))

const LatexCLI = Main.LatexCLI
const GRAPHENE_INPUT = joinpath(@__DIR__, "..", "examples", "graphene", "input.toml")
const CSV3SB5_INPUT = abspath(joinpath(@__DIR__, "..", "examples", "csv3sb5", "input.toml"))
const CSV3SB5_HR = abspath(joinpath(@__DIR__, "..", "examples", "csv3sb5", "wannier90_hr.dat"))
const CSV3SB5_WIN = abspath(joinpath(@__DIR__, "..", "examples", "csv3sb5", "wannier90.win"))
const RUCL3_WIN = abspath(joinpath(@__DIR__, "..", "examples", "rucl3", "wannier90.win"))
const MINIMAL_HR = abspath(joinpath(@__DIR__, "fixtures", "minimal_hr.dat"))
const MINIMAL_WIN = abspath(joinpath(@__DIR__, "fixtures", "minimal.win"))

function _find_line(text::AbstractString, needle::AbstractString)
    for line in split(text, '\n')
        occursin(needle, line) && return line
    end
    return nothing
end

function _orbital_input(hr_path::AbstractString, out_path::AbstractString; max_nn::Int=1)
    return """
[files]
hr = "$(hr_path)"
out = "$(out_path)"

[render]
atol = 1.0e-6
include_second_quant = true

[selection]
mode = "orbital"
wannier_indices = [1, 2]
max_nn = $(max_nn)

[structure]
lattice = [
  [1.0, 0.0, 0.0],
  [0.0, 1.0, 0.0],
  [0.0, 0.0, 1.0],
]

[[orbitals]]
name = "A"
position_frac = [0.0, 0.0, 0.0]
orbitals = ["s"]

[[orbitals]]
name = "B"
position_frac = [0.5, 0.0, 0.0]
orbitals = ["s"]

[[orbitals]]
name = "C"
position_frac = [0.0, 0.5, 0.0]
orbitals = ["s"]
"""
end

@testset "latex module" begin
    @testset "input parsing reads top-level tables" begin
        cfg = LatexCLI.InputIO.read_input(GRAPHENE_INPUT)

        @test endswith(cfg.input_path, joinpath("examples", "graphene", "input.toml"))
        @test endswith(cfg.hr_path, joinpath("examples", "graphene", "graphene_hr.dat"))
        @test isnothing(cfg.win_path)
        @test isnothing(cfg.structure_path)
        @test endswith(cfg.out_path, joinpath("examples", "graphene", "reports", "graphene_hr_hamiltonian.tex"))
        @test cfg.orbital_labels == ["C:pz", "C:pz"]
        @test cfg.atol == 1.0e-6
        @test cfg.max_rblocks == 12
        @test cfg.include_second_quant
    end

    @testset "service renders a latex document from TOML config" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "input.latex.toml")
        out_path = joinpath(tmpdir, "graphene_report.tex")
        hr_path = abspath(joinpath(@__DIR__, "..", "examples", "graphene", "graphene_hr.dat"))

        write(input_path, """
[files]
hr = "$(hr_path)"
structure = ""
out = "$(out_path)"

[render]
atol = 1.0e-6
max_rblocks = 6
include_second_quant = true
""")

        cfg = LatexCLI.InputIO.read_input(input_path)
        result = LatexCLI.Service.run(cfg)
        tex = read(result.output_path, String)

        @test result.hr.num_wann == 2
        @test length(result.entries) >= 6
        @test result.label_source == "default labels"
        @test occursin(raw"\section*{Tight-Binding Hamiltonian}", tex)
        @test occursin(raw"\subsection*{Second-quantization form}", tex)
        @test occursin(raw"\caption{Hopping parameters extracted", tex)
        @test !occursin(raw"\mathbf{a}_1 &= ", tex)

        first_table_row = _find_line(tex, raw"$A$ & $B$ & 0 & 0 & 0")
        @test !isnothing(first_table_row)
        @test endswith(something(first_table_row, ""), " \\\\")

        @test occursin(raw"  &\quad +   f_{A \to B}(\mathbf{k})", tex)

        structure_tail = _find_line(tex, raw"- 2.7\,e^{i\mathbf{k}\cdot(-\mathbf{a}_2)}")
        @test !isnothing(structure_tail)
    end

    @testset "timestamp hr headers fall back to the hr filename" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "input.latex.toml")
        out_path = joinpath(tmpdir, "wannier90_report.tex")
        hr_path = abspath(joinpath(@__DIR__, "..", "examples", "cs2cu3snf12", "wannier90_hr.dat"))
        structure_path = abspath(joinpath(@__DIR__, "..", "examples", "cs2cu3snf12", "POSCAR"))

        if !all(isfile(path) for path in (hr_path, structure_path))
            @test_skip false
        else
            write(input_path, """
[files]
hr = "$(hr_path)"
structure = "$(structure_path)"
out = "$(out_path)"

[render]
atol = 1.0e-6
max_rblocks = 6
include_second_quant = false
""")

            cfg = LatexCLI.InputIO.read_input(input_path)
            result = LatexCLI.Service.run(cfg)
            tex = read(result.output_path, String)

            @test occursin(raw"\texttt{wannier90\_hr.dat}", tex)
            @test !occursin("written on 27Aug2023", tex)
            @test occursin(raw"\mathbf{a}_1 &= (7.0997, -3.5712, 0)", tex)
            @test occursin(raw"\mathbf{a}_3 &= (5.30336, 0, 5.91892)", tex)
        end
    end

    @testset "writer lattice renders a1 a2 a3 values" begin
        cfg = LatexCLI.InputIO.read_input(GRAPHENE_INPUT)
        result = LatexCLI.Service.run(cfg)
        tex = read(result.output_path, String)

        @test occursin(raw"The direct-lattice basis vectors used in $\mathbf{R} = R_1\mathbf{a}_1 + R_2\mathbf{a}_2 + R_3\mathbf{a}_3$ are", tex)
        @test occursin(raw"\mathbf{a}_1 &= (2.46, 0, 0)", tex)
        @test occursin(raw"\mathbf{a}_2 &= (1.23, 2.13042, 0)", tex)
        @test occursin(raw"\mathbf{a}_3 &= (0, 0, 20)", tex)
    end

    @testset "explicit matrix form is not rendered into latex output" begin
        hoppings = Dict((0, 0, 0) => reshape(ComplexF64[1.0 + 0im], 1, 1))
        hr = LatexCLI.LatexCore.WannierTypes.HrBlocks("synthetic", 1, 1, hoppings)
        entries = [
            LatexCLI.LatexCore.Model.HoppingEntry((0, 0, 0), 1, 1, 1.0 + 0im, nothing, nothing),
            LatexCLI.LatexCore.Model.HoppingEntry((1, 0, 0), 1, 1, 2.0 + 0im, nothing, nothing),
            LatexCLI.LatexCore.Model.HoppingEntry((0, 1, 0), 1, 1, 3.0 + 0im, nothing, nothing),
            LatexCLI.LatexCore.Model.HoppingEntry((0, 0, 1), 1, 1, 4.0 + 0im, nothing, nothing),
        ]

        tex = LatexCLI.LatexCore.Render.build_latex(
            hr,
            entries;
            orbital_labels=["A:s"],
            active_indices=[1],
            include_second_quant=false,
        )

        @test !occursin(raw"\subsection*{Explicit matrix form}", tex)
        @test !occursin(raw"Collecting the $R$-resolved contributions", tex)
    end

    @testset "structure factors wrap after three terms" begin
        hoppings = Dict((0, 0, 0) => reshape(ComplexF64[0.0 + 0im 0.0 + 0im; 0.0 + 0im 0.0 + 0im], 2, 2))
        hr = LatexCLI.LatexCore.WannierTypes.HrBlocks("synthetic", 2, 1, hoppings)
        entries = [
            LatexCLI.LatexCore.Model.HoppingEntry((0, 0, 0), 1, 2, 2.0 + 0im, nothing, nothing),
            LatexCLI.LatexCore.Model.HoppingEntry((1, 0, 0), 1, 2, 3.0 + 0im, nothing, nothing),
            LatexCLI.LatexCore.Model.HoppingEntry((0, 1, 0), 1, 2, 4.0 + 0im, nothing, nothing),
            LatexCLI.LatexCore.Model.HoppingEntry((0, 0, 1), 1, 2, 5.0 + 0im, nothing, nothing),
        ]

        tex = LatexCLI.LatexCore.Render.build_latex(
            hr,
            entries;
            orbital_labels=["A:s", "B:s"],
            active_indices=[1, 2],
            include_second_quant=true,
        )

        @test occursin("f_{\\text{A:}s \\to \\text{B:}s}(\\mathbf{k}) &= 2\\,", tex)
        @test occursin("+ 4\\,e^{i\\mathbf{k}\\cdot(\\mathbf{a}_2)} \\\\", tex)
        @test occursin("  &\\quad + 5\\,e^{i\\mathbf{k}\\cdot(\\mathbf{a}_3)}", tex)
    end

    @testset "win file can provide lattice and orbital labels" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "input.latex.toml")
        out_path = joinpath(tmpdir, "win_report.tex")
        win_path = joinpath(tmpdir, "wannier90.win")
        hr_path = abspath(joinpath(@__DIR__, "..", "examples", "graphene", "graphene_hr.dat"))

        write(win_path, """
begin unit_cell_cart
ang
2.46 0.0 0.0
1.23 2.130422493 0.0
0.0 0.0 20.0
end unit_cell_cart

begin projections
f=0.0,0.0,0.0:pz
f=0.3333333333,0.3333333333,0.0:pz
end projections
""")

        write(input_path, """
[files]
hr = "$(hr_path)"
win = "$(win_path)"
out = "$(out_path)"

[render]
atol = 1.0e-6
max_rblocks = 6
include_second_quant = true
""")

        cfg = LatexCLI.InputIO.read_input(input_path)
        result = LatexCLI.Service.run(cfg)
        tex = read(result.output_path, String)

        @test result.label_source == "from Wannier90 projections $(cfg.win_path)"
        @test result.orbital_labels == ["P1:pz", "P2:pz"]
        @test occursin(raw"\mathbf{a}_1 &= (2.46, 0, 0)", tex)
        @test occursin(raw"\mathbf{a}_2 &= (1.23, 2.13042, 0)", tex)
        @test occursin(raw"$\text{P1:}p_z$ & $\text{P2:}p_z$", tex)
    end

    @testset "spin labels render with arrows" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "input.latex.toml")
        out_path = joinpath(tmpdir, "spinful_report.tex")
        hr_path = abspath(joinpath(@__DIR__, "..", "examples", "graphene", "graphene_hr.dat"))

        write(input_path, """
[files]
hr = "$(hr_path)"
structure = ""
out = "$(out_path)"

[render]
atol = 1.0e-6
max_rblocks = 6
include_second_quant = true

[[orbitals]]
name = "A"
orbitals = ["pz"]
spins = ["up"]

[[orbitals]]
name = "B"
orbitals = ["pz"]
spins = ["down"]
""")

        cfg = LatexCLI.InputIO.read_input(input_path)
        result = LatexCLI.Service.run(cfg)
        tex = read(result.output_path, String)

        @test occursin(raw"$\text{A:}p_z\uparrow$ & $\text{B:}p_z\downarrow$", tex)
        @test occursin(raw"c^\dagger_{\text{A:}p_z\uparrow,\mathbf{k}}", tex)
        @test occursin(raw"orbital labels: $\text{A:}p_z\uparrow$, $\text{B:}p_z\downarrow$", tex)
        @test !occursin(raw"\text{:up}", tex)
        @test !occursin(raw"\text{:down}", tex)
    end

    @testset "cli argument parsing" begin
        opts = LatexCLI.parse_args(["input.toml"])
        @test opts.input_path == "input.toml"
        @test !opts.show_help

        opts_flag = LatexCLI.parse_args(["--input", "config.toml"])
        @test opts_flag.input_path == "config.toml"
    end

    @testset "unknown top-level keys are rejected" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "bad-top-level.toml")
        write(input_path, """
junk = "value"

[files]
hr = "graphene_hr.dat"
out = "report.tex"
""")

        @test_throws ErrorException LatexCLI.InputIO.read_input(input_path)
    end

    @testset "unknown file keys are rejected" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "bad-files.toml")
        write(input_path, """
[files]
hr = "graphene_hr.dat"
out = "report.tex"
extra = "unused"
""")

        @test_throws ErrorException LatexCLI.InputIO.read_input(input_path)
    end

    @testset "unknown render keys are rejected" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "bad-render.toml")
        write(input_path, """
[files]
hr = "graphene_hr.dat"
out = "report.tex"

[render]
atol = 1.0e-8
unknown = 3
""")

        @test_throws ErrorException LatexCLI.InputIO.read_input(input_path)
    end

    @testset "shorthand orbital expansion supports f orbitals in inline input" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "f-inline.toml")
        write(input_path, """
[files]
hr = "$(MINIMAL_HR)"
out = "$(joinpath(tmpdir, "unused.tex"))"

[[orbitals]]
name = "Fe"
position_frac = [0.0, 0.0, 0.0]
orbitals = ["f"]
""")

        cfg = LatexCLI.InputIO.read_input(input_path)
        @test cfg.orbital_labels == [
            "Fe:fz3",
            "Fe:fxz2",
            "Fe:fyz2",
            "Fe:fz(x2-y2)",
            "Fe:fxyz",
            "Fe:fx(x2-3y2)",
            "Fe:fy(3x2-y2)",
        ]
    end

    @testset "win site mapping expands d and f projections in Wannier90 order" begin
        specs = LatexCLI.LatexCore.WannierWinIO.read_win_orbital_specs(MINIMAL_WIN, 12)
        @test length(specs) == 12
        @test [spec.label for spec in specs[1:5]] == [
            "Fe1:dz2",
            "Fe1:dxz",
            "Fe1:dyz",
            "Fe1:dx2-y2",
            "Fe1:dxy",
        ]
        @test [spec.label for spec in specs[6:12]] == [
            "Fe2:fz3",
            "Fe2:fxz2",
            "Fe2:fyz2",
            "Fe2:fz(x2-y2)",
            "Fe2:fxyz",
            "Fe2:fx(x2-3y2)",
            "Fe2:fy(3x2-y2)",
        ]
        @test specs[1].position_frac == (0.0, 0.0, 0.0)
        @test specs[6].position_frac == (0.5, 0.0, 0.0)
    end

    @testset "win species projections apply to every matching atom when listed once" begin
        specs = LatexCLI.LatexCore.WannierWinIO.read_win_orbital_specs(CSV3SB5_WIN, 30)
        @test length(specs) == 30
        @test [spec.label for spec in specs[1:5]] == [
            "V1:dz2",
            "V1:dxz",
            "V1:dyz",
            "V1:dx2-y2",
            "V1:dxy",
        ]
        @test all(spec.position_frac == specs[1].position_frac for spec in specs[1:5])
        @test all(spec.position_frac == specs[6].position_frac for spec in specs[6:10])
        @test all(spec.position_frac == specs[11].position_frac for spec in specs[11:15])
        @test [spec.label for spec in specs[6:10]] == ["V2:dz2", "V2:dxz", "V2:dyz", "V2:dx2-y2", "V2:dxy"]
        @test [spec.label for spec in specs[11:15]] == ["V3:dz2", "V3:dxz", "V3:dyz", "V3:dx2-y2", "V3:dxy"]
        @test [spec.label for spec in specs[16:18]] == ["Sb1:pz", "Sb1:px", "Sb1:py"]
        @test all(spec.position_frac == specs[16].position_frac for spec in specs[16:18])
        @test [spec.label for spec in specs[28:30]] == ["Sb5:pz", "Sb5:px", "Sb5:py"]
        @test length(unique([spec.position_frac for spec in specs])) == 8
    end

    @testset "win inline fractional projections ignore axis metadata" begin
        specs = LatexCLI.LatexCore.WannierWinIO.read_win_orbital_specs(RUCL3_WIN, 6)
        @test length(specs) == 6
        @test [spec.label for spec in specs] == [
            "P1:dxz",
            "P1:dyz",
            "P1:dxy",
            "P2:dxz",
            "P2:dyz",
            "P2:dxy",
        ]
    end

    @testset "orbital mode keeps only active basis and adds NN metadata" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "orbital.toml")
        out_path = joinpath(tmpdir, "orbital.tex")
        write(input_path, _orbital_input(MINIMAL_HR, out_path))

        cfg = LatexCLI.InputIO.read_input(input_path)
        result = LatexCLI.Service.run(cfg)
        tex = read(result.output_path, String)

        @test result.active_indices == [1, 2]
        @test all(entry.m in (1, 2) && entry.n in (1, 2) for entry in result.entries)
        @test all(!isnothing(entry.nn_shell) for entry in result.entries)
        @test occursin("NN & distance", tex)
        @test !occursin(raw"\subsection*{Explicit matrix form}", tex)
        @test !occursin(raw"\text{C:}s", tex)
    end

    @testset "max_nn zero keeps only raw home-cell zero-distance terms" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "orbital-onsite.toml")
        out_path = joinpath(tmpdir, "orbital-onsite.tex")
        write(input_path, _orbital_input(MINIMAL_HR, out_path; max_nn=0))

        cfg = LatexCLI.InputIO.read_input(input_path)
        result = LatexCLI.Service.run(cfg)

        @test !isempty(result.entries)
        @test all(entry.R == (0, 0, 0) && entry.m == entry.n for entry in result.entries)
        @test Set(entry.m for entry in result.entries) == Set([1, 2])
    end

    @testset "orbital mode keeps raw R and uses physical bond distances" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "orbital-physical.toml")
        out_path = joinpath(tmpdir, "orbital-physical.tex")
        write(input_path, _orbital_input(MINIMAL_HR, out_path))

        cfg = LatexCLI.InputIO.read_input(input_path)
        result = LatexCLI.Service.run(cfg)

        onsite = filter(entry -> entry.R == (0, 0, 0) && entry.m == 1 && entry.n == 1, result.entries)
        @test length(onsite) == 1
        @test something(only(onsite).distance) == 0.0
        @test something(only(onsite).nn_shell) == 0

        left_bond = filter(entry -> entry.R == (-1, 0, 0) && entry.m == 1 && entry.n == 2, result.entries)
        @test length(left_bond) == 1
        @test something(only(left_bond).distance) == 0.5
        @test something(only(left_bond).nn_shell) == 1

        home_bond = filter(entry -> entry.R == (0, 0, 0) && entry.m == 1 && entry.n == 2, result.entries)
        @test length(home_bond) == 1
        @test something(only(home_bond).distance) == 0.5
        @test something(only(home_bond).nn_shell) == 1
    end

    @testset "orbital mode uses global distance shells like Kwant" begin
        tmpdir = mktempdir()
        input_path = joinpath(tmpdir, "csv3sb5-maxnn3.toml")
        out_path = joinpath(tmpdir, "csv3sb5-maxnn3.tex")
        input_text = """
[files]
hr = "$(CSV3SB5_HR)"
win = "$(CSV3SB5_WIN)"
out = "$(out_path)"

[render]
atol = 1.0e-6
include_second_quant = true

[selection]
mode = "orbital"
wannier_indices = [5, 10, 15]
max_nn = 3
"""
        write(input_path, input_text)

        cfg = LatexCLI.InputIO.read_input(input_path)
        result = LatexCLI.Service.run(cfg)

        first_shell = filter(entry -> something(entry.nn_shell) == 1, result.entries)
        third_shell = filter(entry -> something(entry.nn_shell) == 3, result.entries)

        @test !isempty(first_shell)
        @test all(isapprox(something(entry.distance), 2.74745; atol=1e-6) for entry in first_shell)
        @test any(isapprox(something(entry.distance), 5.4949; atol=1e-6) for entry in third_shell)
        @test !any(something(entry.nn_shell) == 1 && isapprox(something(entry.distance), 5.4949; atol=1e-6) for entry in result.entries)
    end
end
