#=
これは輪講第２回目の時間非依存のシュレディンガー方程式における
ポテンシャル０の場合の固有値を求め、数値解と解析解を比較するプログラムである。
main() -> 固有値比較
main2() -> 固有関数比較
=#

using LinearAlgebra
using Plots

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

function main()
    V(x) = 0
    N = 1000
    L = 1
    H = make_H(N,L,V)
    e = eigvals(H)
    
    
    e0 = zeros(Float64,N)
    for n=1:N
        e0[n] = (n^2*π^2)/L^2
    end
    plot(1:N,e,labels="Numerial result",xlabel="n",ylabel="energy")
    plot!(1:N,e0,labels="Analytical result",xlabel="n",ylabel="energy")
    savefig(joinpath(@__DIR__,"eigen.png"))
    println(e0[1],"\t",e[1])
    for n=1:N
        if e0[n] > 10^6
            println("10^6を超えるn=",n)
            println("解析解",e0[n],"\t","数値解",e[n]) 
            break
        end
    end
end

function main2()
    V(x) = 0
    N = 1000
    L = 1
    
    H = make_H(N,L,V)
    v = eigvecs(H)

    Δx = L/(N + 1)
    xs = zeros(Float64,N)
    ψ0 = zeros(Float64,N)
    ψ250 = similar(ψ0)
    n , m = 1 , 250
    
    for i = 1:N
        x = i*Δx
        xs[i] = x
        ψ0[i] = sqrt(2/L)*sin(x*n*π/L)
        ψ250[i] = sqrt(2/L)*sin(x*m*π/L)
    end
    coeff = 1/sqrt(Δx)
    plot(xs,coeff*v[:,n],label = "Numerial result(n=1)",xlabel="x",ylabel="ψ(x)")
    plot!(xs,ψ0,label = "Analytical result(n=1)",xlabel="x",ylabel="ψ(x)")
    savefig(joinpath(@__DIR__,"psi(n=1).png"))

    plot(xs,coeff*v[:,m],label = "Numerial result(n=250)",xlabel="x",ylabel="ψ(x)")
    plot!(xs,ψ250,label = "Analytical result(n=250)",xlabel="x",ylabel="ψ(x)")
    savefig(joinpath(@__DIR__,"psi(n=250).png"))
end 
main()
#main2()