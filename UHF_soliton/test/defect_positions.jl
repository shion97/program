using Test
ENV["GKSwstype"] = "100"
include(joinpath(@__DIR__, "..", "experiment_defect_positions.jl"))

@testset "位置比較: 粒子数、収束、PBC平行移動" begin
    base_up, base_down = experiment_density(5, 24)
    shifted_up, shifted_down = experiment_density(15, 24)
    wrapped_up, wrapped_down = experiment_density(39, 24)
    @test base_up[4:6] ≈ [0.6, 0.6, 0.4]
    @test shifted_up ≈ circshift(base_up, 10)
    @test wrapped_up ≈ circshift(base_up, 34)
    for (up, down) in ((base_up, base_down), (shifted_up, shifted_down), (wrapped_up, wrapped_down))
        @test sum(up) ≈ 25
        @test sum(down) ≈ 25
        @test up + down ≈ ones(50)
    end
    @test wall_bonds(base_up - base_down) == [4, 28]
    base_result = redirect_stdout(devnull) do
        solve_uhf(base_up, base_down, N_up, N_down;
            t=t, U=U, boundary=boundary, mixing=mixing, tol=tol, maxiter=maxiter)
    end
    shifted_result = redirect_stdout(devnull) do
        solve_uhf(shifted_up, shifted_down, N_up, N_down;
            t=t, U=U, boundary=boundary, mixing=mixing, tol=tol, maxiter=maxiter)
    end
    for result in (base_result, shifted_result)
        @test result.converged
        @test result.residual < tol
        @test sum(result.n_up) ≈ 25 atol=tol
        @test sum(result.n_down) ≈ 25 atol=tol
    end
    @test shifted_result.energy ≈ base_result.energy atol=tol
    @test shifted_result.n_up ≈ circshift(base_result.n_up, 10) atol=tol
    @test shifted_result.n_down ≈ circshift(base_result.n_down, 10) atol=tol
end
