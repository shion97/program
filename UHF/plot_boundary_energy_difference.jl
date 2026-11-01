ENV["GKSwstype"] = "100"
using Test, Printf, Plots

function report_match(text, pattern, source)
    matches = collect(eachmatch(pattern, text))
    length(matches) == 1 || throw(ArgumentError("$source: 必須項目が欠けているか重複しています: $pattern"))
    return only(matches)
end

function parse_uhf_report(text; source="入力テキスト")
    text = replace(text, "\r\n" => "\n")
    parameter(name) = report_match(text, Regex("^" * name * " = ([^\\n]+)\$", "m"), source).captures[1]
    parameters = (
        L=parse(Int, parameter("L")), t=parse(Float64, parameter("t")),
        U=parse(Float64, parameter("U")), N_up=parse(Int, parameter("N_up")),
        N_down=parse(Int, parameter("N_down")), mixing=parse(Float64, parameter("mixing")),
        tol=parse(Float64, parameter("tol")), maxiter=parse(Int, parameter("maxiter")),
    )
    boundary = report_match(text, r"^boundary = \"(OBC|PBC)\"$"m, source).captures[1]
    converged = report_match(text, r"^converged=(true|false)$"m, source).captures[1] == "true"
    final = report_match(text, r"^iterations=(\d+)\s+residual=(\S+)\s+E=(\S+)$"m, source)
    iterations = parse(Int, final.captures[1])
    residual, energy = parse.(Float64, final.captures[2:3])
    all(isfinite, (parameters..., residual, energy)) || throw(ArgumentError("$source: 非有限値があります"))
    (parameters.L > 0 && parameters.t > 0 && parameters.tol > 0 &&
     0 <= parameters.N_up <= parameters.L && 0 <= parameters.N_down <= parameters.L &&
     0 < parameters.mixing <= 1 && 1 <= iterations <= parameters.maxiter && residual >= 0) ||
        throw(ArgumentError("$source: 計算設定または最終残差が不正です"))
    return (; parameters, boundary, converged, iterations, residual, energy)
end

function compare_pair(obc, pbc)
    (obc.boundary == "OBC" && pbc.boundary == "PBC") || throw(ArgumentError("OBC・PBCの指定が不正です"))
    obc.parameters == pbc.parameters || throw(ArgumentError("OBC・PBCの計算パラメータが一致しません"))
    invalid = String[]
    for report in (obc, pbc)
        (!report.converged || report.residual >= report.parameters.tol) && push!(invalid, report.boundary)
    end
    eligible = isempty(invalid)
    return (; eligible, difference=eligible ? pbc.energy - obc.energy : NaN, reason=join(invalid, "/"))
end

function experiment_specs()
    return [
        (number=1, xlabel="Condition", values=[1.0], labels=["reference"],
         groups=[(label="", folders=[""])]),
        (number=2, xlabel="Condition", values=[1.0], labels=["reference"],
         groups=[(label="", folders=[""])]),
        (number=3, xlabel="N_e", values=[46.,48.,50.,52.,54.], labels=["46","48","50","52","54"],
         groups=[(label="", folders=["Ne46","Ne48","Ne50","Ne52","Ne54"])]),
        (number=4, xlabel="Delta N", values=[0.,2.,4.,12.], labels=["0","2","4","12"],
         groups=[(label="", folders=["dN00","dN02","dN04","dN12"])]),
        (number=5, xlabel="N_e", values=[46.,54.], labels=["46","54"],
         groups=[(label="Delta N=0", folders=["Ne46_dN00","Ne54_dN00"]),
                 (label="Delta N=4", folders=["Ne46_dN04","Ne54_dN04"])]),
        (number=6, xlabel="Initial density", values=[1.,2.,3.,4.], labels=["uniform","AF","domain_wall","random"],
         groups=[(label="", folders=["uniform","AF","domain_wall","random"])]),
        (number=7, xlabel="Mixing", values=[0.10,0.25,0.50], labels=["0.10","0.25","0.50"],
         groups=[(label="", folders=["mixing_010","mixing_025","mixing_050"])]),
        (number=8, xlabel="U/t", values=[0.,1.,2.,4.,6.,8.], labels=["0","1","2","4","6","8"],
         groups=[(label="", folders=["U00","U01","U02","U04","U06","U08"])]),
    ]
end

function collect_experiments(root=joinpath(@__DIR__, "results_vault"))
    experiments = NamedTuple[]
    for spec in experiment_specs()
        folder = joinpath(root, @sprintf("experiment_%02d", spec.number))
        groups = NamedTuple[]
        for group in spec.groups
            length(group.folders) == length(spec.values) == length(spec.labels) ||
                throw(ArgumentError("実験$(spec.number): 条件名と横軸の数が一致しません"))
            comparisons = NamedTuple[]
            for condition in group.folders
                reports = map(("OBC", "PBC")) do boundary
                    path = joinpath(folder, condition, boundary, "uhf_results.txt")
                    isfile(path) || throw(ArgumentError("入力がありません: $path"))
                    parse_uhf_report(read(path, String); source=path)
                end
                comparison = compare_pair(reports...)
                push!(comparisons, comparison)
            end
            push!(groups, (;label=group.label, comparisons))
        end
        push!(experiments, (;spec, folder, groups))
    end
    return experiments
end

function draw_energy_difference(experiment)
    spec = experiment.spec
    differences = [comparison.difference for group in experiment.groups for comparison in group.comparisons if comparison.eligible]
    lower, upper = extrema(vcat(0.0, differences))
    span = upper > lower ? upper - lower : 1.0
    spacing = length(spec.values) > 1 ? minimum(diff(spec.values)) : 1.0
    chart = plot(;
        title="Experiment $(spec.number): boundary energy difference\nNC = not converged; no energy difference plotted",
        xlabel=spec.xlabel, ylabel="E_PBC - E_OBC", size=(1000,600),
        xticks=(spec.values, spec.labels),
        xlims=(first(spec.values)-spacing/2, last(spec.values)+spacing/2),
        ylims=(lower-span/4, upper+span/10), legend=length(experiment.groups)>1 ? :topright : false,
        titlefontsize=12, markerstrokewidth=1,
        left_margin=10Plots.mm, bottom_margin=6Plots.mm, top_margin=6Plots.mm,
    )
    hline!(chart, [0.0]; label="", color=:gray, linestyle=:dash)
    for (group_index, group) in enumerate(experiment.groups)
        eligible = [comparison.eligible for comparison in group.comparisons]
        color = group_index == 1 ? :steelblue : :darkorange
        if any(eligible)
            scatter!(chart, spec.values[eligible], [comparison.difference for comparison in group.comparisons if comparison.eligible];
                label=group.label, color, markersize=6)
        end
        for (condition_index, comparison) in enumerate(group.comparisons)
            if !comparison.eligible
                note = isempty(group.label) ? "NC" : "NC\n$(group.label)"
                annotate!(chart, spec.values[condition_index], lower-span/8, text(note, 9, color, :center))
            end
        end
    end
    return chart
end

function main()
    # 入力不備で既存図が一部だけ更新されないよう、全ログと全図を先に評価する。
    experiments = collect_experiments()
    charts = draw_energy_difference.(experiments)
    for (experiment, chart) in zip(experiments, charts)
        destination = joinpath(experiment.folder, "boundary_energy_difference.png")
        savefig(chart, destination)
        println("Saved: ", destination)
        for group in experiment.groups, (index, comparison) in enumerate(group.comparisons)
            if !comparison.eligible
                println("  NC: ", experiment.spec.labels[index], " ", group.label, " (", comparison.reason, ")")
            end
        end
    end
end

function example_report(; boundary="OBC", energy=-3.0, converged=true, residual=1e-11)
    return """
    # パラメータ
    L = 50
    t = 1.0
    U = 4.0
    N_up = 25
    N_down = 25
    boundary = "$boundary"
    # 計算の設定
    mixing = 0.25
    tol = 1.0e-10
    maxiter = 10000
    # 初期電子密度
    initial_mode = :manual
    # 計算結果
    iter=    1  residual=1.000e-01  E=-99.000000000000  E_kin=-100.0  E_int=1.0
    converged=$converged
    iterations=126  residual=$residual  E=$energy
    """
end

function run_self_tests()
    @testset "最終結果の読み取り" begin
        report = parse_uhf_report(replace(example_report(), "\n" => "\r\n"))
        @test report.energy == -3.0
        @test report.boundary == "OBC"
        @test report.converged
        @test report.parameters.tol == 1e-10
        @test_throws ArgumentError parse_uhf_report(replace(example_report(), "iterations=126" => "missing=126"))
        @test_throws ArgumentError parse_uhf_report(example_report(energy=NaN))
    end
    @testset "境界条件の比較と未収束" begin
        obc = parse_uhf_report(example_report())
        pbc = parse_uhf_report(example_report(boundary="PBC", energy=-3.5))
        @test compare_pair(obc, pbc).difference == -0.5
        @test compare_pair(obc, pbc).eligible
        for invalid in (merge(pbc, (;converged=false)), merge(pbc, (;residual=1e-10)))
            comparison = compare_pair(obc, invalid)
            @test !comparison.eligible
            @test isnan(comparison.difference)
            @test comparison.reason == "PBC"
        end
        @test_throws ArgumentError compare_pair(pbc, obc)
        @test_throws ArgumentError compare_pair(obc, merge(pbc, (;parameters=merge(pbc.parameters, (;U=8.0)))))
    end
    @testset "再利用を含む実験の配置" begin
        specs = experiment_specs()
        @test length(specs) == 8
        @test sum(length(group.folders) for spec in specs for group in spec.groups) == 28
        @test specs[3].values == [46, 48, 50, 52, 54]
        @test specs[7].labels == ["0.10", "0.25", "0.50"]
        @test specs[8].groups[1].folders == ["U00", "U01", "U02", "U04", "U06", "U08"]
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if ARGS == ["--test"]
        run_self_tests()
    elseif isempty(ARGS)
        main()
    else
        error("実行方法: julia plot_boundary_energy_difference.jl [--test]")
    end
end
