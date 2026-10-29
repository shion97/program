using LinearAlgebra
using DelimitedFiles
using Printf
using Plots

# =============== 計算設定（ここを編集） ===============
# パラメータ
L = 50                # 電子数
t = 1.0                # ホッピング積分
U = 4.0               # オンサイト相互作用
m = 25                # N_up=m, N_down=L-m （0 <= m <= L）
N_up = m
N_down = L - m

# ファイルの扱い
initial_mode = :load  # :manual か :load を選ぶ
input_csv = "results_vault/反強磁性/density_checkpoint.csv"    # pathは相対pathで入力
output_csv = "density_checkpoint.csv"
result_dir_name = nothing

# 計算の設定
mixing = 0.25   # 密度混合割合
tol = 1e-10   # 収束判定 
maxiter = 10_000    # 最大反復回数

# png出力ファイル名
energy_plot_file = "scf_energy.png"
residual_plot_file = "scf_residual.png"
orbital_plot_file = "final_orbital_energies.png"
mag_plot_file = "final_magnetization.png"
report_file = "uhf_results.txt"

# ===========================================

# =============== 入出力用補助 ===============
function result_directory()
    result_dir_name === nothing && return joinpath(pwd(), "outputs")
    result_dir_name isa AbstractString || error("result_dir_nameはnothingまたはフォルダ名にしてください")
    (!isempty(result_dir_name) && result_dir_name ∉ (".", "..") &&
     !occursin(r"[\\/:]", result_dir_name)) || error("result_dir_nameには区切り文字を含まないフォルダ名を指定してください")
    return joinpath(pwd(), "results_vault", result_dir_name)
end

abs_path(name::AbstractString) = joinpath(result_directory(), basename(name))

function print_result(message::AbstractString, report_io)
    print(stdout, message)
    if report_io !== nothing
        print(report_io, message)
        flush(report_io)
    end
end

function write_initial_conditions(report_io)
    println(report_io, "# パラメータ")
    for (name, value) in (("L", L), ("t", t), ("U", U), ("m", m),
                          ("N_up", N_up), ("N_down", N_down))
        println(report_io, name, " = ", repr(value))
    end
    println(report_io, "\n# 計算の設定")
    for (name, value) in (("mixing", mixing), ("tol", tol), ("maxiter", maxiter))
        println(report_io, name, " = ", repr(value))
    end
    println(report_io, "\n# 初期電子密度")
    println(report_io, "initial_mode = ", repr(initial_mode))
    if initial_mode == :manual
        source_lines = readlines(@__FILE__)
        manual_start = findlast(==("    if initial_mode == :manual"), source_lines)
        manual_start === nothing && error("記録するmanual分岐が見つかりません")
        manual_end = findnext(==("    elseif initial_mode == :load"), source_lines, manual_start + 1)
        manual_end === nothing && error("記録するmanual分岐の終端が見つかりません")
        for line in source_lines[manual_start + 1:manual_end - 1]
            println(report_io, replace(line, r"^        " => ""))
        end
    elseif initial_mode == :load
        println(report_io, "input_csv = ", repr(relpath(abspath(input_csv), pwd())))
    end
    println(report_io, "\n# 計算結果")
    flush(report_io)
end

# =============== 入力・計算用関数 ===============
function validate_density(up, down, L, N_up, N_down)
    length(up) == L && length(down) == L || error("密度の配列長はL=$Lにしてください")
    all(isfinite, up) && all(isfinite, down) || error("密度には有限の数値を指定してください")
    all(x -> 0 <= x <= 1, up) || error("各n_up[i]は0～1である必要があります")
    all(x -> 0 <= x <= 1, down) || error("各n_down[i]は0～1である必要があります")
    isapprox(sum(up), N_up; atol=1e-8, rtol=0) || error("sum(n_up)=$(sum(up)), 指定N_up=$N_up")
    isapprox(sum(down), N_down; atol=1e-8, rtol=0) || error("sum(n_down)=$(sum(down)), 指定N_down=$N_down")
end

function save_density(path, up, down)
    open(path, "w") do io
        println(io, "site,n_up,n_down")
        for i in eachindex(up)
            @printf(io, "%d,%.17g,%.17g\n", i, up[i], down[i])
        end
    end
end

function load_density(path, L)
    raw = readdlm(path, ',', Float64; skipstart=1)
    size(raw) == (L, 3) || error("CSVは\$L行、3列(site,n_up,n_down)にしてください")
    all(Int.(raw[:, 1]) .== collect(1:L)) || error("CSVのsite列は1からLの順にしてください")
    return Vector{Float64}(raw[:, 2]), Vector{Float64}(raw[:, 3])
end

function site_density(vectors, N_occ)
    N_occ == 0 && return zeros(size(vectors, 1))
    return vec(sum(abs2, @view(vectors[:, 1:N_occ]); dims=2))
end

function output_density(h0, U, up, down, N_up, N_down)
    eig_up = eigen(Symmetric(h0 + Diagonal(U .* down)))
    eig_down = eigen(Symmetric(h0 + Diagonal(U .* up)))
    next_up = site_density(eig_up.vectors, N_up)
    next_down = site_density(eig_down.vectors, N_down)
    return next_up, next_down, eig_up, eig_down
end

function hf_energy(h0, U, eig_up, eig_down, dens_up, dens_down, input_up, input_down, N_up, N_down)
    psi_up = N_up == 0 ? zeros(eltype(eig_up.vectors), size(eig_up.vectors,1), 0) : @view eig_up.vectors[:, 1:N_up]
    psi_down = N_down == 0 ? zeros(eltype(eig_down.vectors), size(eig_down.vectors,1), 0) : @view eig_down.vectors[:, 1:N_down]
    E_kin = sum(psi_up .* (h0 * psi_up)) + sum(psi_down .* (h0 * psi_down))
    E_int = U * sum(dens_up .* dens_down)
    E_total = E_kin + E_int
    # 固有値は input_up / input_down の平均場から計算されている。
    # 非収束時も成り立つ恒等式を使ってエネルギーを検算する。
    E_check = sum(eig_up.values[1:N_up]) + sum(eig_down.values[1:N_down])
    E_check -= U * sum(input_down .* dens_up .+ input_up .* dens_down)
    E_check += E_int
    @assert abs(E_total - E_check) < 1e-7
    return E_total, E_kin, E_int
end

function solve_uhf(up_init, down_init, N_up, N_down;
                   t=1.0, U=4.0, mixing=0.25, tol=1e-10, maxiter=10_000,
                   report_io=nothing)
    L = length(up_init)
    0 <= N_up <= L && 0 <= N_down <= L || error("各スピンの粒子数は0～L")
    N_up + N_down == L || error("half-fillingはN_up+N_down=L")
    0 < mixing <= 1 || error("mixingは0より大きく1以下")
    tol > 0 && maxiter >= 1 || error("tol、maxiterを確認してください")
    up = Float64.(up_init)
    down = Float64.(down_init)
    validate_density(up, down, L, N_up, N_down)

    h0 = zeros(Float64, L, L)
    for i in 1:L-1
        h0[i, i+1] = -t
        h0[i+1, i] = -t
    end

    iter_hist = Int[]
    residual_hist = Float64[]
    energy_hist = Float64[]

    converged = false
    iterations = 0
    for k in 1:maxiter
        next_up, next_down, eig_up, eig_down = output_density(h0, U, up, down, N_up, N_down)
        residual = max(maximum(abs.(next_up .- up)), maximum(abs.(next_down .- down)))
        E_total, E_kin, E_int = hf_energy(h0, U, eig_up, eig_down, next_up, next_down, up, down, N_up, N_down)

        push!(iter_hist, k)
        push!(residual_hist, residual)
        push!(energy_hist, E_total)

        iterations = k
        print_result(@sprintf("iter=%5d  residual=%.3e  E=%.12f  E_kin=%.12f  E_int=%.12f\n",
                              k, residual, E_total, E_kin, E_int), report_io)

        if residual < tol
            up .= next_up
            down .= next_down
            converged = true
            break
        end

        up .= (1-mixing) .* up .+ mixing .* next_up
        down .= (1-mixing) .* down .+ mixing .* next_down
    end
    final_up, final_down, eig_up, eig_down = output_density(h0, U, up, down, N_up, N_down)
    final_residual = max(maximum(abs.(final_up .- up)), maximum(abs.(final_down .- down)))
    E_total, E_kin, E_int = hf_energy(h0, U, eig_up, eig_down, final_up, final_down, up, down, N_up, N_down)

    return (
        converged=converged,
        n_up=final_up,
        n_down=final_down,
        charge=final_up .+ final_down,
        magnetization=final_up .- final_down,
        orbital_energies_up=eig_up.values,
        orbital_energies_down=eig_down.values,
        energy=E_total,
        E_kin=E_kin,
        E_int=E_int,
        residual=final_residual,
        iterations=iterations,
        iter_hist=iter_hist,
        residual_hist=residual_hist,
        energy_hist=energy_hist,
    )
end

function save_plots(result; energy_path, residual_path, orbital_path, mag_path)
    sites = 1:length(result.n_up)

    p1 = plot(result.iter_hist, result.energy_hist;
        xlabel="Iteration", ylabel="Energy", label="E_total",
        title="SCF energy history", linewidth=2, marker=:circle)
    savefig(p1, energy_path)

    p2 = plot(result.iter_hist, result.residual_hist;
        xlabel="Iteration", ylabel="Residual", label="residual",
        title="SCF residual history", linewidth=2, marker=:circle,
        yscale=:log10)
    savefig(p2, residual_path)

    orbitals = 1:length(result.orbital_energies_up)
    p3 = plot(orbitals, result.orbital_energies_up;
        xlabel="Orbital index (ascending energy)", ylabel="Orbital energy", label="Spin up",
        title="Final orbital energies", linewidth=2, marker=:circle)
    plot!(p3, orbitals, result.orbital_energies_down;
        label="Spin down", linewidth=2, linestyle=:dash, marker=:square,
        markercolor=:transparent)
    savefig(p3, orbital_path)

    p4 = plot(sites, result.magnetization;
        xlabel="Site i", ylabel="m_i", label="m_i",
        title="Final local magnetization", linewidth=2, marker=:circle)
    hline!(p4, [0.0]; label="")
    savefig(p4, mag_path)
end

# =============== 実行部分：最初にif文で入力方法を選ぶ ===============
function run_calculation(report_io)
    0 <= m <= L || error("mは0～Lにしてください")
    write_initial_conditions(report_io)

    if initial_mode == :manual
        # 手動指定：以下の2配列を任意の形に変更できます。
        
        up_base = N_up / L
        down_base = N_down / L
        amplitude = 0.2 * min(up_base, 1-up_base, down_base, 1-down_base)
        up0 = [up_base + amplitude * (-1)^i for i in 1:L]
        down0 = [down_base - amplitude * (-1)^i for i in 1:L]
        
        #up0[5:30], down0[5:30] =copy(down0[5:30]), copy(up0[5:30])

        #up0 = [0.5 for _ in 1:L]
        #down0 = [0.5 for _ in 1:L]

    elseif initial_mode == :load
        up0, down0 = load_density(abspath(input_csv), L)
    else
        error("initial_modeは:manualまたは:loadにしてください")
    end

    validate_density(up0, down0, L, N_up, N_down)
    result = solve_uhf(up0, down0, N_up, N_down;
                       t=t, U=U, mixing=mixing, tol=tol, maxiter=maxiter,
                       report_io=report_io)

    print_result(@sprintf("\nN_up=%d N_down=%d  N_total=%d\n", N_up, N_down, N_up+N_down), report_io)
    print_result("converged=$(result.converged)\n", report_io)
    if !result.converged
        print_result("SCFは最大反復回数に達しても未収束です。最後の混合密度から再評価した結果を保存します。\n", report_io)
    end
    print_result(@sprintf("iterations=%d  residual=%.3e  E=%.12f\n",
                         result.iterations, result.residual, result.energy), report_io)
    print_result(@sprintf("E_kin=%.12f  E_int=%.12f\n", result.E_kin, result.E_int), report_io)
    print_result(@sprintf("sum(up)=%.10f  sum(down)=%.10f\n", sum(result.n_up), sum(result.n_down)), report_io)

    save_density(abs_path(output_csv), result.n_up, result.n_down)
    save_plots(result;
        energy_path=abs_path(energy_plot_file),
        residual_path=abs_path(residual_plot_file),
        orbital_path=abs_path(orbital_plot_file),
        mag_path=abs_path(mag_plot_file),
    )

    for filename in (output_csv, energy_plot_file, residual_plot_file, orbital_plot_file, mag_plot_file, report_file)
        println("Saved: ", abs_path(filename))
    end
    return result
end

function main()
    mkpath(dirname(abs_path(report_file)))
    open(abs_path(report_file), "w") do report_io
        run_calculation(report_io)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
