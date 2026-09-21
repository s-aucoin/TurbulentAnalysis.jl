@testset "quality_control" begin
    @test FailsQC(Int8(0b1010), 1) == true
    @test FailsQC(Int8(0b1010), 0) == false
    @test FailsQC(Int8(0b1010), 2) == false
    @test FailsQC(Int8(0b1010), 3) == true

    q = (;
        sf2 = [1.0, -0.5, 0.8, 0.4],
        vari = [0.1, 0.2, 0.1, 0.3],
        sf2σ = [0.5, 0.2, 0.01, 0.3],
        variσ = [0.05, 0.1, 0.01, 0.15],
        sf2n = [0.6, 0.65, 0.8, 0.7],
        sf3 = [0.4, 0.2, 0.5, 0.1],
        sf3σ = [0.1, 0.05, 0.01, 0.2],
        sf3n = [0.9, 1.1, 1.0, 0.4],
    )

    qc = QCϵts(q, 0.5, 0.1; exponent_correction = 0)
    @test haskey(qc, :qsf2)
    @test haskey(qc, :qvari)
    @test haskey(qc, :qsf3)
    @test all(length.(values(qc)) .== length(q.sf2))
end
