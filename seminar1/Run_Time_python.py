"""
これは輪講第1回目の実行速度の比較用プログラムです
"""

import time

def loop():
    s = 1
    for i in range(1000000):
        s += i
    return s

# 1回実行
t1 = time.time()
loop()
t2 = time.time()
print("Python (1 call):", t2 - t1)

# 50回実行
t3 = time.time()
for _ in range(50):
    loop()
t4 = time.time()
print("Python (50 calls):", t4 - t3,)