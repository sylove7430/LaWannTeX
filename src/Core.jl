module LatexCore

include(joinpath(@__DIR__, "model", "Constants.jl"))
include(joinpath(@__DIR__, "model", "CrystalCell.jl"))
include(joinpath(@__DIR__, "model", "WannierTypes.jl"))
include(joinpath(@__DIR__, "model", "Model.jl"))

include(joinpath(@__DIR__, "io", "PoscarIO.jl"))
include(joinpath(@__DIR__, "io", "WannierHrIO.jl"))
include(joinpath(@__DIR__, "io", "InputIO.jl"))
include(joinpath(@__DIR__, "io", "WannierWinIO.jl"))
include(joinpath(@__DIR__, "io", "LabelIO.jl"))
include(joinpath(@__DIR__, "io", "LatticeIO.jl"))

include(joinpath(@__DIR__, "analysis", "Table.jl"))
include(joinpath(@__DIR__, "analysis", "OrbitalSelection.jl"))
include(joinpath(@__DIR__, "render", "Render.jl"))
include(joinpath(@__DIR__, "analysis", "Service.jl"))

end
