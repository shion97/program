/*
これは輪講第1回目の実行速度の比較用プログラムです
*/
#include <windows.h>
#include <stdio.h>

long loop() {
    long s = 1;
    for (int i = 0; i < 100000; i++) s += i;
    return s;
}

double now() {
    LARGE_INTEGER freq, counter;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / freq.QuadPart;
}

int main() {
    double t1 = now();
    loop();
    double t2 = now();
    printf("C (1 call): %f \n", t2 - t1);

    t1 = now();
    for (int i = 0; i < 50; i++) loop();
    t2 = now();
    printf("C (50 calls): %f \n", t2 - t1);

    return 0;
}
