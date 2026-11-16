using Plots
using LaTeXStrings

function read_energy_differences(report_path)
    report = read(report_path, String)
    differences = map(("delta_E", "delta_E_kin", "delta_E_int")) do field
        entries = collect(eachmatch(Regex("(?m)(?:^|\\s)" * field * "=([^\\s]+)"), report))
        isempty(entries) && return nothing
        value = tryparse(Float64, last(entries).captures[1])
        (value === nothing || !isfinite(value)) && error("不正な$(field): $(report_path)")
        value
    end
    any(isnothing, differences) && return nothing
    convergence_entries = collect(eachmatch(r"(?m)^converged=(true|false)\s*$", report))
    isempty(convergence_entries) && error("収束状態がありません: $(report_path)")
    return (; differences, converged=last(convergence_entries).captures[1] == "true")
end

function plot_energy_differences(vault_directory=joinpath(@__DIR__, "results_vault"),
                                 output_path=joinpath(vault_directory, "energy_differences.png"))
    report_paths = String[]
    for (directory, _, filenames) in walkdir(vault_directory)
        "uhf_results.txt" in filenames && push!(report_paths, joinpath(directory, "uhf_results.txt"))
    end
    sort!(report_paths; by=path -> replace(relpath(path, vault_directory), '\\' => '/'))
    labels = String[]
    results = NTuple{3,Float64}[]
    for report_path in report_paths
        result = read_energy_differences(report_path)
        if result === nothing
            @warn "エネルギー差が3種類そろわないため除外: $(report_path)"
            continue
        end
        label = replace(relpath(dirname(report_path), vault_directory), '\\' => '/')
        push!(labels, label * (result.converged ? "" : " (unconverged)"))
        push!(results, result.differences)
    end
    isempty(results) && error("3種類のエネルギー差を持つレポートがありません: $(vault_directory)")

    positions = collect(eachindex(results))
    figure = plot(; title="Energy differences", xlabel="Energy difference",
        yticks=(positions, labels), yflip=true, legend=:outertopright,
        size=(1000, 600), left_margin=8Plots.mm, bottom_margin=10Plots.mm)
    group_offset = 0.18
    for (component, label) in enumerate((raw"$\Delta E$", raw"$\Delta E_{\mathrm{kin}}$", raw"$\Delta E_{\mathrm{int}}$"))
        scatter!(figure, [values[component] for values in results],
            positions .+ (component - 2) * group_offset; label,
            color=(:steelblue, :darkorange, :seagreen)[component],
            markershape=(:circle, :diamond, :utriangle)[component], markersize=7)
    end
    vline!(figure, [0.0]; color=:gray, linestyle=:dash, linewidth=1, label="")
    mkpath(dirname(abspath(output_path)))
    savefig(figure, output_path)
    return figure
end

if abspath(PROGRAM_FILE) == @__FILE__
    ENV["GKSwstype"] = "100"
    plot_energy_differences()
    println("保存先: ", joinpath(@__DIR__, "results_vault", "energy_differences.png"))
end
