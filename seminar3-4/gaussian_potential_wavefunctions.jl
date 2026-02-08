#=
これは輪講第３・４回目の時間非依存のシュレディンガー方程式における
ガウス型ポテンシャルに対する固有値問題を解き、
波動関数をプロットするプログラムである。
=#

using LinearAlgebra
using Plots

#v(q)の計算
function calc_vq(q, ξ)
    return sqrt(π * ξ^2) * exp(-q^2 * ξ^2 / 4)
end

#Vmnの計算
function calc_Vmn(k, kp, L, ξ, x0, V0)
    q1 = k - kp
    q2 = k + kp
    vq1 = calc_vq(q1, ξ)
    vq2 = calc_vq(q2, ξ)
    Vkkp = (V0 / L) * (cos(q1 * x0) * vq1 + cos(q2 * x0) * vq2)
    return Vkkp
end

#ハミルトニアン行列の作成
function make_H(N, L, ξ, x0, V0)
    H = zeros(Float64, N, N)
    for n in 1:N
        k = n * π / L
        for m in 1:N
            kp = m * π / L
            H[n, m] = (n == m ? k^2 : 0.0) +
                      calc_Vmn(k, kp, L, ξ, x0, V0)
        end
    end
    return H
end

#波動関数ψ(x)の計算
function calc_psi(cn, x, L)
    psi = 0.0
    for n in eachindex(cn)
        kn = n * π / L
        psi += cn[n] * sin(kn * x)
    end
    return psi * sqrt(2 / L)
end

#波動関数のプロット
function momentumspace(N, L, ξ, x0, V0)
    Hk = make_H(N, L, ξ, x0, V0)

    eig = eigen(Hk)
    ep = eig.values
    bn = eig.vectors

    xs = range(0, L, length=N)
    psi = zeros(Float64, N)

    n = 1  # 基底状態
    for (i, x) in enumerate(xs)
        psi[i] = calc_psi(bn[:, n], x, L)
    end

    plot(xs, psi,
        label = "result\n(n = 1)",
        xlabel = "x",
        ylabel = "ψ(x)"
    )

    savefig(joinpath(@__DIR__, "momu.png"))
end

function main()
    N = 1000      
    L = 10
    ξ = 1
    x0 = L / 2
    V0 = 10

    momentumspace(N, L, ξ, x0, V0)
end

main()