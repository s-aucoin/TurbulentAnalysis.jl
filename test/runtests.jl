using TurbulentAnalysis
using Test
using Dates

@testset "TurbulentAnalysis" begin
    include("velocity_structure_functions.jl")
    include("quality_control.jl")
    include("misc.jl")
end