ENV["GKSwstype"] = "100"
include(joinpath(@__DIR__, "uhf.jl"))
using SHA

const POSITION_CASES = [
    (name="reference", first_site=1, segment_length=0),
    (name="gap02", first_site=5, segment_length=2),
    (name="gap04", first_site=5, segment_length=4),
    (name="gap08", first_site=5, segment_length=8),
    (name="gap16", first_site=5, segment_length=16),
    (name="gap24", first_site=5, segment_length=24),
    (name="shift24", first_site=15, segment_length=24),
    (name="wrap24", first_site=39, segment_length=24),
]
const COMPARISON_DIRECTORY = joinpath(@__DIR__, "results_vault", "defect_positions_20261106")

function experiment_density(first_site, segment_length)
    @assert iseven(L) && iseven(segment_length) && 0 <= segment_length < L
    up_base, down_base = N_up / L, N_down / L
    amplitude = 0.2 * min(up_base, 1-up_base, down_base, 1-down_base)
    staggered = Float64[(-1)^site for site in 1:L]
    staggered .-= sum(staggered) / L
    initial_up = up_base .+ amplitude .* staggered
    initial_down = down_base .- amplitude .* staggered
    swapped_sites = mod1.(first_site:first_site+segment_length-1, L)
    initial_up[swapped_sites], initial_down[swapped_sites] =
        copy(initial_down[swapped_sites]), copy(initial_up[swapped_sites])
    validate_density(initial_up, initial_down, L, N_up, N_down)
    return initial_up, initial_down
end

function wall_bonds(magnetization)
    staggered_magnetization = [(-1)^site * magnetization[site] for site in eachindex(magnetization)]
    return [site for site in eachindex(magnetization)
        if staggered_magnetization[site] * staggered_magnetization[mod1(site+1, length(magnetization))] < 0]
end

function save_experiment_density(path, up, down)
    open(path, "w") do density_io
        print(density_io, '\ufeff', "site,n_up,n_down\n")
        for site in eachindex(up)
            @printf(density_io, "%d,%.17g,%.17g\n", site, up[site], down[site])
        end
    end
end

function run_position_comparison()
    @assert (L, t, U, N_up, N_down, boundary, mixing, tol, maxiter) ==
        (50, 1.0, 4.0, 25, 25, "PBC", 0.25, 1e-8, 10000) "計画時の条件が変更されています"
    # 既存の結果を再実行で上書きしないため、ケースディレクトリの存在時は停止する。
    any(ispath(joinpath(COMPARISON_DIRECTORY, case.name)) for case in POSITION_CASES) &&
        error("比較結果が既に存在します。別の保存先を指定してください")
    mkpath(COMPARISON_DIRECTORY)
    source_path = joinpath(@__DIR__, "uhf.jl")
    cp(source_path, joinpath(COMPARISON_DIRECTORY, "uhf_source_snapshot.jl"); force=false)
    write(joinpath(COMPARISON_DIRECTORY, "source_sha256.txt"), bytes2hex(sha256(read(source_path))) * "  uhf.jl\n")
    results = Dict{String, Any}()
    initial_states = Dict{String, Any}()
    for case in POSITION_CASES
        case_directory = joinpath(COMPARISON_DIRECTORY, case.name)
        mkpath(case_directory)
        initial_up, initial_down = experiment_density(case.first_site, case.segment_length)
        initial_states[case.name] = (initial_up, initial_down)
        save_experiment_density(joinpath(case_directory, "initial_density.csv"), initial_up, initial_down)
        result = open(joinpath(case_directory, "uhf_results.txt"), "w") do report_io
            for setting in (:L, :t, :U, :N_up, :N_down, :boundary, :mixing, :tol, :maxiter)
                println(report_io, setting, " = ", repr(getfield(@__MODULE__, setting)))
            end
            println(report_io, "initial_mode = :manual\ncase = ", repr(case.name))
            println(report_io, "first_site = ", case.first_site, "\nsegment_length = ", case.segment_length)
            println(report_io, "initial_density = initial_density.csv\ninitial_amplitude = 0.1")
            println(report_io, "initial_wall_bonds = ", wall_bonds(initial_up - initial_down))
            calculation = solve_uhf(initial_up, initial_down, N_up, N_down;
                t=t, U=U, boundary=boundary, mixing=mixing, tol=tol, maxiter=maxiter, report_io=report_io)
            println(report_io, "\nconverged = ", calculation.converged)
            println(report_io, "iterations = ", calculation.iterations, "\nresidual = ", calculation.residual)
            println(report_io, "energy = ", calculation.energy, "\nE_kin = ", calculation.E_kin, "\nE_int = ", calculation.E_int)
            println(report_io, "sum_up = ", sum(calculation.n_up), "\nsum_down = ", sum(calculation.n_down))
            println(report_io, "final_wall_bonds = ", wall_bonds(calculation.magnetization))
            calculation
        end
        results[case.name] = result
        save_experiment_density(joinpath(case_directory, "density_checkpoint.csv"), result.n_up, result.n_down)
        save_plots(result, initial_up, initial_down;
            energy_path=joinpath(case_directory, "scf_energy.png"),
            residual_path=joinpath(case_directory, "scf_residual.png"),
            orbital_path=joinpath(case_directory, "final_orbital_energies.png"),
            initial_density_path=joinpath(case_directory, "initial_density_magnetization.png"),
            final_density_path=joinpath(case_directory, "final_density_magnetization.png"))
    end

    reference_energy = results["reference"].energy
    open(joinpath(COMPARISON_DIRECTORY, "comparison.csv"), "w") do comparison_io
        print(comparison_io, '\ufeff', "case,first_site,segment_length,min_wall_distance,converged,iterations,residual,energy,delta_energy_reference,sum_up,sum_down,max_abs_charge_deviation,final_wall_bonds,min_abs_magnetization\n")
        for case in POSITION_CASES
            result = results[case.name]
            println(comparison_io, join((case.name, case.first_site, case.segment_length,
                min(case.segment_length, L-case.segment_length), result.converged, result.iterations,
                result.residual, result.energy, result.energy-reference_energy,
                sum(result.n_up), sum(result.n_down), maximum(abs.(result.charge .- 1)),
                join(wall_bonds(result.magnetization), ";"), minimum(abs.(result.magnetization))), ','))
        end
    end
    open(joinpath(COMPARISON_DIRECTORY, "translation_checks.csv"), "w") do check_io
        print(check_io, '\ufeff', "case,shift,abs_delta_energy,max_abs_aligned_up,max_abs_aligned_down,max_abs_delta_orbitals\n")
        base_result = results["gap24"]
        for (case_name, shift) in (("shift24", 10), ("wrap24", 34))
            result = results[case_name]
            differences = (abs(result.energy-base_result.energy),
                maximum(abs.(result.n_up .- circshift(base_result.n_up, shift))),
                maximum(abs.(result.n_down .- circshift(base_result.n_down, shift))),
                max(maximum(abs.(result.orbital_energies_up .- base_result.orbital_energies_up)),
                    maximum(abs.(result.orbital_energies_down .- base_result.orbital_energies_down))))
            println(check_io, join((case_name, shift, differences...), ','))
            @assert all(difference < tol for difference in differences)
        end
    end

    staggered = Float64[(-1)^site for site in 1:L]
    panels = []
    for case in POSITION_CASES
        initial_up, initial_down = initial_states[case.name]
        result = results[case.name]
        panel = plot(1:L, staggered .* (initial_up-initial_down);
            label="initial", linestyle=:dash, color=:gray, linewidth=2,
            title="$(case.name): d=$(min(case.segment_length, L-case.segment_length)), walls=$(length(wall_bonds(result.magnetization)))",
            xlabel="Site i", ylabel="(-1)^i m_i", ylims=(-1, 1),
            legend=case.name == "reference" ? :bottomright : false)
        plot!(panel, 1:L, staggered .* result.magnetization; label="final", linewidth=2)
        push!(panels, panel)
    end
    savefig(plot(panels...; layout=(4,2), size=(1400,1600)), joinpath(COMPARISON_DIRECTORY, "staggered_comparison.png"))
    charge_plot = plot(; xlabel="Site i", ylabel="n_up + n_down - 1", title="Final charge modulation", legend=:outerright)
    for case in POSITION_CASES
        plot!(charge_plot, 1:L, results[case.name].charge .- 1; label=case.name, linewidth=2)
    end
    savefig(charge_plot, joinpath(COMPARISON_DIRECTORY, "charge_comparison.png"))
    separation_cases = POSITION_CASES[2:6]
    energy_plot = plot([case.segment_length for case in separation_cases],
        [results[case.name].energy-reference_energy for case in separation_cases];
        xlabel="Initial shortest wall distance", ylabel="E - E(reference)",
        title="Converged energy vs initial wall distance", marker=:circle, label="same first wall", linewidth=2)
    hline!(energy_plot, [0.0]; label="reference", linestyle=:dash)
    savefig(energy_plot, joinpath(COMPARISON_DIRECTORY, "energy_comparison.png"))
    @assert all(result.converged && result.residual < tol for result in values(results))
    @assert all(isapprox(sum(result.n_up), N_up; atol=tol, rtol=0) &&
        isapprox(sum(result.n_down), N_down; atol=tol, rtol=0) for result in values(results))
    println("COMPARISON_COMPLETE: ", COMPARISON_DIRECTORY)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_position_comparison()
end
