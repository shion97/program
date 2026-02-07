#=
これは輪講第５回目の時間依存のシュレディンガー方程式における
自由粒子の時間発展を数値的に求め、GIFとして保存するプログラムである。
=#

using LinearAlgebra
using Plots

function make_H!(H,N,L,V)
    @.H = 0
    Δx = L/(N + 1)
    
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
    return 
end

function timevolv(ψ,N,Nt,Δt,H)
    ψs = zeros(ComplexF64,N,Nt)
    U = exp(-im * Δt * H)

    for i in 1:Nt
        ψ = U * ψ
        ψs[:,i] = ψ
    end
    return ψs
end

function timedep_simple()
    anim = Animation()
    N = 4000
    L = 40.0
    xs = range(0,L,length = N)
    σ = 1
    k0 = 10
    ψ0 = zeros(ComplexF64,N)
    x0 = 5
    @.ψ0 = (π*σ^2)^(-1/4)*exp(-(xs-x0)^2 / (2σ^2) + im*k0*(xs-x0))
    dx = (xs[2] - xs[1])
    V(x) = 0
    H = zeros(Float64,N,N)
    make_H!(H,N,L,V)
    Δt = 0.02

    Nt = 100
    #規格化定数
    c = sqrt(norm(ψ0)^2 * dx)
    ψ = ψ0 / c
    ψ2max = maximum(abs.(ψ).^2)

    println("norm = $(norm(ψ)^2 * dx)")
    println("ψ2max = $(ψ2max)")

    @time ψs = timevolv(ψ,N,Nt,Δt,H)

    for i in 1:Nt
        plt = plot(xs,abs.(ψs[:,i].^2),ylim=(0,ψ2max),xlabel="x",ylabel="|ψ|^2")
        println("$i-th: norm = $(norm(ψs[:,i])^2 * dx)")
        frame(anim,plt)
    end
    gif(anim,"simple.gif",fps = 30)
    
end

#実行関数
#timedep_simple()
