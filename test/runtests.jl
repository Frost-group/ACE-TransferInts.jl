using ACETransferInts
using Test
using LinearAlgebra

@testset "ACETransferInts.jl" begin

    @testset "rotation matrices" begin
        for R in (rot_x, rot_y, rot_z), θ in (0.0, π / 7, π / 2, 2.3)
            @test R(θ)' * R(θ) ≈ I(3)
            @test det(R(θ)) ≈ 1.0
        end
        @test rot_z(π / 2) * [1.0, 0.0, 0.0] ≈ [0.0, 1.0, 0.0]
        @test rot_x(π / 2) * [0.0, 1.0, 0.0] ≈ [0.0, 0.0, 1.0]
        @test rot_y(π / 2) * [0.0, 0.0, 1.0] ≈ [1.0, 0.0, 0.0]
        @test rot_y(0.4) * rot_y(-0.4) ≈ I(3)
    end

    @testset "centroid" begin
        X = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 3.0, 0.0]]
        @test centroid(X) ≈ [1.0, 1.0, 0.0]
    end

    @testset "metrics" begin
        ref  = [1.0, -2.0, 3.0, -4.0]
        pred = [1.5, -2.0, -3.0, -4.5]
        @test rmse(ref, ref) == 0.0
        @test mae(pred, ref) ≈ (0.5 + 0.0 + 6.0 + 0.5) / 4
        @test rmse(pred, ref) ≈ sqrt((0.25 + 0.0 + 36.0 + 0.25) / 4)
        @test sign_accuracy(pred, ref) ≈ 0.75
    end

    @testset "learning-curve CSV" begin
        csv = joinpath(@__DIR__, "..", "experiments", "ethylene", "results",
                       "LC_Ethylene_CSi_od2_td10_p19_q45_rcut75.csv")
        n, rmse_v, mae_v, sign_acc = read_lc_csv(csv)
        @test n == [100, 250, 500, 1000, 1500, 2000, 3000, 5000]
        @test length(rmse_v) == length(mae_v) == length(sign_acc) == length(n)
        @test all(mae_v .<= rmse_v)          # MAE ≤ RMSE always
        @test rmse_v[end] < rmse_v[1]        # more data, lower error
    end
end
