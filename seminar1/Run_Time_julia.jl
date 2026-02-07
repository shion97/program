#=
これは輪講第1回目の実行速度の比較用プログラムです
=#
function loop()
    s = 1
    for i in 1:1000000
        s += i
    end
    return s
end

#1回実行の計測(JIT含む)
t1 = time()
loop()
t2 = time()
println("Julia (1 call): ", t2 - t1)

#50回実行の計測(JIT除く)
t3 = time()
for _ in 1:50
    loop()
end
t4 = time()
println("Julia (50 calls): ", t4 - t3)