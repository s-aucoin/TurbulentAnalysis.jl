using Random

@testset "velocity_structure_functions" begin
    @test DLLx(8.0) ≈ 8.0^(2/3)
    @test DLLLx(8.0) ≈ 8.0
    @test ϵ_DLL(2.0) ≈ 1.0
    @test ϵ_DLLL(2.0) ≈ -2.5
    @test σ_DLL(4.0) ≈ 2.0
    @test TurbulentAnalysis.TimeDuration([0.0, 2.0, 5.0]) == 5.0

    t = DateTime(2020, 1, 1):Second(1):DateTime(2020, 1, 1, 0, 0, 3)
    @test TurbulentAnalysis.TimeDuration(collect(t)) == 3.0

    x = collect(0.0:1.0:4.0)
    ur = hcat([0.0, 1.0, 0.0, -1.0, 0.0], [0.5, 1.5, 0.2, -1.2, 0.3])
    umom = umoments(ur, 1.0, x, 2.0; min_points = 3, calculate_sf3 = true)
    @test haskey(umom, :ML2)
    @test haskey(umom, :ML3)
    @test length(umom.x) == length(umom.ML2)
    @test all(length.(umom.ML2) .>= 1)
    @test all(length.(umom.ML3) .>= 1)

    rng = MersenneTwister(1234)
    Δr = 0.1
    rmax = 0.5
    x_long = collect(0.0:Δr:10.0)
    ur_long = randn(rng, length(x_long), 8)
    umom_long = umoments(ur_long, Δr, x_long, rmax)
    ϵ_long = ϵofx(umom_long)
    sf2_finite = filter(isfinite, ϵ_long.sf2_ϵ)
    sf3_finite = filter(isfinite, ϵ_long.sf3_ϵ)
    @test !isempty(sf2_finite)
    @test !isempty(sf3_finite)

    Δt = 0.1
    fs = 1/Δt
    Δt_new = 0.5
    t_series = collect(0.0:Δt:12.0)
    ur_ts = randn(rng, length(x_long), length(t_series))
    ϵ_ts = ϵofxandt(ur_ts, x_long, t_series, Δr, fs, Δt_new, rmax; xmax = minimum(x), xmin = maximum(x))
    @test length(ϵ_ts.t) == size(ϵ_ts.sf2, 2)
    @test size(ϵ_ts.sf2, 1) == length(ϵ_ts.x)
    @test !isempty(filter(isfinite, ϵ_ts.sf2))
    @test !isempty(filter(isfinite, ϵ_ts.sf3))
end
