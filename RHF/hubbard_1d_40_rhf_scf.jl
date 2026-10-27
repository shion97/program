# Integer indices and loop counters remain Int, as appropriate.
using LinearAlgebra
using Printf
using Plots

# 1D Hubbard model, restricted Hartree-Fock, open boundary conditions.
# Half filling: L electrons (L/2 up and L/2 down).

# Kinetic (nearest-neighbor hopping) matrix.
function hopping_matrix(L::Int, t::Float64)
    T = zeros(Float64, L, L)
    for i in 1:(L - 1)
        T[i, i + 1] = -t
        T[i + 1, i] = -t
    end
    return T
end

# Site density n_i = 2 * sum_{a=1}^{L/2} |phi_a(i)|^2.
# Each spatial orbital contains one up- and one down-spin electron.
function site_density(Cocc::AbstractMatrix{Float64})
    return 2.0 .* vec(sum(abs2, Cocc; dims=2))
end

# RHF total energy for the determinant defined by occupied orbitals.
function energy_components(T::AbstractMatrix{Float64}, Cocc::AbstractMatrix{Float64}, n::AbstractVector{Float64}, U::Float64)
    Ekin = 2.0 * sum(Cocc .* (T * Cocc))
    Eint = (U / 4.0) * sum(abs2, n)
    return Ekin + Eint, Ekin, Eint
end

function solve_rhf(L::Int, t::Float64, U::Float64; alpha::Float64=0.2, tol::Float64=1e-10, maxiter::Int=1000)
    iseven(L) || error("L must be even at half filling in closed-shell RHF.")
    t > 0 || error("t must be positive.")
    U >= 0 || error("This example assumes repulsive U >= 0.")
    0 < alpha <= 1 || error("alpha must satisfy 0 < alpha <= 1.")

    nocc = L ÷ 2
    T = hopping_matrix(L, t)

    # Deliberately start away from n_i=1, without breaking spin symmetry.
    n = Float64[1.0 + 0.1 * (isodd(i) ? 1.0 : -1.0) for i in 1:L]

    history = Float64[]
    converged = false
    iterations = 0

    for iter in 1:maxiter
        # F_ij = T_ij + (U/2) * n_i * delta_ij.
        F = T + Diagonal((U / 2.0) .* n)

        # Eigenvalues are in ascending order. Columns are eigenvectors.
        eig = eigen(Symmetric(F))
        Cocc = eig.vectors[:, 1:nocc]

        # Density calculated from the newly occupied orbitals.
        n_new = site_density(Cocc)
        residual = maximum(abs.(n_new .- n))
        push!(history, residual)

        if iter/10 < 1 || iter % 10 == 0 || residual < tol
            E, _, _ = energy_components(T, Cocc, n_new, U)
            @printf("iter = %4d   residual = %.3e   E = %.12f\n",
                    iter, residual, E)
        end

        iterations = iter
        if residual < tol
            n = n_new
            converged = true
            break
        end

        # Linear density mixing stabilizes the self-consistency iteration.
        n = (1.0 - alpha) .* n .+ alpha .* n_new
    end

    converged || error("SCF did not converge: try a smaller alpha.")

    # Rebuild and diagonalize once more using the final mixed density.
    F = T + Diagonal((U / 2.0) .* n)
    eig = eigen(Symmetric(F))
    Cocc = eig.vectors[:, 1:nocc]
    n_final = site_density(Cocc)
    E, Ekin, Eint = energy_components(T, Cocc, n_final, U)
    final_residual = maximum(abs.(n_final .- n))

    return (
        density=n_final,
        eigenvalues=eig.values,
        orbitals=Cocc,
        energy=E,
        kinetic=Ekin,
        interaction=Eint,
        iterations=iterations,
        residual=final_residual,
        history=history,
    )
end

# Parameters: 40 sites, 40 electrons, t = 1, U = 4.
L = 40
t = Float64(1.0)
U = Float64(4.0)
result = solve_rhf(L, t, U; alpha=0.2, tol=1e-10)

println("\n---- Converged RHF result ----")
@printf("Sites                : %d\n", L)
@printf("Electrons            : %d (up = %d, down = %d)\n", L, L ÷ 2, L ÷ 2)
@printf("SCF iterations       : %d\n", result.iterations)
@printf("Final residual       : %.3e\n", result.residual)
@printf("Total electron number: %.12f\n", sum(result.density))
@printf("Max |n_i - 1|        : %.3e\n", maximum(abs.(result.density .- 1.0)))
@printf("Kinetic energy       : %.12f\n", result.kinetic)
@printf("Interaction energy   : %.12f\n", result.interaction)
@printf("Total RHF energy     : %.12f\n", result.energy)

# Independent analytic check for the OPEN chain at half filling.
E_reference = -4.0 * t * sum(cos(m * pi / (L + 1)) for m in 1:(L ÷ 2)) + U * L / 4.0
@printf("Analytic RHF energy  : %.12f\n", E_reference)
@printf("Energy difference    : %.3e\n", abs(result.energy - E_reference))

# HOMO-LUMO one-particle level separation (not a many-body charge gap).
gap = result.eigenvalues[L ÷ 2 + 1] - result.eigenvalues[L ÷ 2]
@printf("HOMO-LUMO separation : %.12f\n", gap)

println("\nSelected site densities:")
for i in 1:40
    @printf("site %2d : n_i = %.12f\n", i, result.density[i])
end

# 各サイトの電子密度
p1 = plot(
    1:L,
    result.density,
    marker=:circle,
    xlabel="Site i",
    ylabel="Electron density n_i",
    ylim=(0.95, 1.05),
    label="RHF density",
    title="Site density"
)

# SCF 残差の推移
p2 = plot(
    1:length(result.history),
    result.history,
    marker=:circle,
    yscale=:log10,
    xlabel="SCF iteration",
    ylabel="Density residual",
    label="SCF",
    title="SCF convergence"
)

# 2つのグラフを上下に並べる
fig = plot(p1, p2; layout=(2, 1), size=(800, 600))

display(fig)
savefig(fig, "rhf_scf.png")