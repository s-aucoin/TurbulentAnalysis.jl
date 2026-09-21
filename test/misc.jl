@testset "misc" begin
    x = collect(0.0:1.0:4.0)
    t = collect(0.0:0.1:1.0)
    u = randn(length(x), length(t))

    eddies = largest_eddies(u, x, t, 1.0, 0.4, 10; xmin = minimum(x), xmax = maximum(x))

    @test length(eddies.times) == length(eddies.ell)
    @test length(eddies.u_rms) == length(eddies.ϵ_rough)
    @test all(isfinite, filter(!isnan, eddies.ell))
    @test all(isfinite, filter(!isnan, eddies.u_rms))
end
