#=
これは輪講第３・４回目の時間非依存のシュレディンガー方程式における
ガウス型ポテンシャルに対する固有値問題を解き、
数値解とフーリエ変換を用いた解法の結果を比較するプログラムである。
=#

using LinearAlgebra
using Plots

#H行列の作成
function make_H(N,L,V)
    Δx = L/(N + 1)
    H = zeros(N,N)
    for i = 1:N 
        x = i * Δx
        H[i,i] = V(x)

        j = i + 1
        dij = -1/Δx^2
        if 1 <= j <= N
            H[i,j] += dij
        end
        j = i
        dij = 2/Δx^2
        if 1 <= j <= N
            H[i,j] += dij
        end
        j = i - 1
        dij = -1/Δx^2
        if 1 <= j <= N
            H[i,j] += dij
        end

    end
    return H
end

#ガウス型ポテンシャルの計算
function gaussV(x,ξ,x0,V0)
    return V0 * exp(-(x-x0)^2/ξ^2 )
end

#------------------------------------

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
function make_Hk(N, L, ξ, x0, V0)
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

#フーリエ変換での固有エネルギーの計算
function momentumspace(N, L, ξ, x0, V0)
    Hk = make_Hk(N, L, ξ, x0, V0)

    eig = eigen(Hk)
    ep = eig.values
    bn = eig.vectors
    return ep
end

function main()
    N = 1000
    L = 10
    V(x) = gaussV(x,1,L/2,10)
    H = make_H(N,L,V)
    e = eigvals(H)
    plot(1:N,e,labels="Numerial result",xlabel="n",ylabel="energy")

    e = momentumspace(N, L, 1, L/2, 10)
    plot!(1:N,e,labels="Momentum space result",xlabel="n",ylabel="energy")
    savefig(joinpath(@__DIR__,"eigen_comparision.png"))

end

main()