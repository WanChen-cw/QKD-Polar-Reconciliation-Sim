#include <iostream>
#include "PolarTest.hpp"
#include <chrono>
#include <thread>
#include <omp.h>
int main() {

    double sum = 0;
    double add = 1;
    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();
    const int max_threads = omp_get_max_threads(); // 获取最大线程数

    // 使用 OpenMP 并行创建线程
    int num_threads = 20;
#pragma omp parallel num_threads(num_threads)
   {
        // 获取当前线程的编号
        int thread_id = omp_get_thread_num();
        /*std::cout << "Thread " << thread_id << " is active." << std::endl;*/

        // 调用测试函数
        fastpolartest();

    }


    // Stop measuring time and calculate the elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

    printf("Result: %.20f\n", sum);
    printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
    return 0;

}