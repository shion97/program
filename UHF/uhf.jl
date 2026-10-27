using LinearAlgebra
using DelimitedFiles
using Printf
using Plots

# =============== 計算設定（ここを編集） ===============
L = 40
t = 1.0
U = 4.0
m = 20                # N_up=m, N_down=L-m （0 <= m <= L）
N_up = m
N_down = L - m

initial_mode = :manual  # :manual か :load を選ぶ
input_csv = "density_checkpoint.csv"
output_csv = "density_checkpoint.csv"

mixing = 0.25
tol = 1e-10
maxiter = 10_000

# png出力ファイル名
energy_plot_file = "scf_energy.png"
residual_plot_file = "scf_residual.png"
density_plot_file = "final_density.png"
mag_plot_file = "final_magnetization.png"

# =============== 入出力用補助 ===============
# Julia実行時のカレントディレクトリ配下の outputs に読込・保存を統一
output_dir = joinpath(pwd(), "outputs")
mkpath(output_dir)  # outputs が存在しなければ作成する
abs_path(name::AbstractString) = joinpath(output_dir, basename(name))

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
                   t=1.0, U=4.0, mixing=0.25, tol=1e-10, maxiter=10_000)
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
        @printf("iter=%5d  residual=%.3e  E=%.12f  E_kin=%.12f  E_int=%.12f\n",
                k, residual, E_total, E_kin, E_int)

        if residual < tol
            up .= next_up
            down .= next_down
            converged = true
            break
        end

        up .= (1-mixing) .* up .+ mixing .* next_up
        down .= (1-mixing) .* down .+ mixing .* next_down
    end
    converged || error("SCFが収束しませんでした。mixingや初期密度を変更してください")

    final_up, final_down, eig_up, eig_down = output_density(h0, U, up, down, N_up, N_down)
    final_residual = max(maximum(abs.(final_up .- up)), maximum(abs.(final_down .- down)))
    E_total, E_kin, E_int = hf_energy(h0, U, eig_up, eig_down, final_up, final_down, up, down, N_up, N_down)

    return (
        n_up=final_up,
        n_down=final_down,
        charge=final_up .+ final_down,
        magnetization=final_up .- final_down,
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

function save_plots(result; energy_path, residual_path, density_path, mag_path)
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

    p3 = plot(sites, result.n_up;
        xlabel="Site i", ylabel="Density", label="n_up",
        title="Final electron density", linewidth=2, marker=:circle)
    plot!(p3, sites, result.n_down; label="n_down", linewidth=2, marker=:square)
    plot!(p3, sites, result.charge; label="charge", linewidth=2, marker=:diamond)
    savefig(p3, density_path)

    p4 = plot(sites, result.magnetization;
        xlabel="Site i", ylabel="m_i", label="m_i",
        title="Final local magnetization", linewidth=2, marker=:circle)
    hline!(p4, [0.0]; label="")
    savefig(p4, mag_path)
end

# =============== 実行部分：最初にif文で入力方法を選ぶ ===============
0 <= m <= L || error("mは0～Lにしてください")

if initial_mode == :manual
    # 手動指定：以下の2配列を任意の形に変更できます。
    # デフォルトはN_up=m、N_down=L-mを満たす交代磁化の微小な種です。
    
    up_base = N_up / L
    down_base = N_down / L
    amplitude = 0.2 * min(up_base, 1-up_base, down_base, 1-down_base)
    up0 = [up_base + amplitude * (-1)^i for i in 1:L]
    down0 = [down_base - amplitude * (-1)^i for i in 1:L]
    
    up0[5:30], down0[5:30] =copy(down0[5:30]), copy(up0[5:30])

    #up0 = [0.5 for _ in 1:L]
    #down0 = [0.5 for _ in 1:L]

elseif initial_mode == :load
    up0, down0 = load_density(abs_path(input_csv), L)
else
    error("initial_modeは:manualまたは:loadにしてください")
end

validate_density(up0, down0, L, N_up, N_down)
result = solve_uhf(up0, down0, N_up, N_down;
                   t=t, U=U, mixing=mixing, tol=tol, maxiter=maxiter)

@printf("\nN_up=%d N_down=%d  N_total=%d\n", N_up, N_down, N_up+N_down)
@printf("iterations=%d  residual=%.3e  E=%.12f\n",
        result.iterations, result.residual, result.energy)
@printf("E_kin=%.12f  E_int=%.12f\n", result.E_kin, result.E_int)
@printf("sum(up)=%.10f  sum(down)=%.10f\n", sum(result.n_up), sum(result.n_down))

save_density(abs_path(output_csv), result.n_up, result.n_down)
save_plots(result;
    energy_path=abs_path(energy_plot_file),
    residual_path=abs_path(residual_plot_file),
    density_path=abs_path(density_plot_file),
    mag_path=abs_path(mag_plot_file),
)

println("Saved: ", abs_path(output_csv))
println("Saved: ", abs_path(energy_plot_file))
println("Saved: ", abs_path(residual_plot_file))
println("Saved: ", abs_path(density_plot_file))
println("Saved: ", abs_path(mag_plot_file))
