#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <bitset>
#include <chrono>
#include <random>
#include <cuda_runtime.h>

__device__ int convertToGrayCode(int n) {
    return n ^ (n >> 1);
}


__global__ void computeGrayCodeAndRowSum(double* d_M, double* d_r, double* d_pm, int n, int nc) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= (1 << n)) return;

    int gray_i = convertToGrayCode(i);
    int gray_prev = convertToGrayCode(i - 1);

    int idx = __ffs(gray_i ^ gray_prev) - 1;
    int sign = (gray_i >> idx) & 1 ? 1 : -1;

    for (int k = 0; k < nc; k++) {
        d_r[i * nc + k] += d_M[idx * nc + k] * sign;
    }

    double s = (i % 2 == 0) ? 1.0 : -1.0;
    double product = s;
    for (int k = 0; k < nc; k++) {
        product *= d_r[i * nc + k];
    }

    atomicAdd(d_pm, product);
}

long double perm3(const std::vector<std::vector<double>>& M) {
    int n = M.size();
    int nc = M[0].size();

    if (n != nc) {
        throw std::runtime_error("Matrix must be square");
    }

    std::vector<double> r = M[0];
    long double pm = -std::accumulate(r.begin(), r.end(), 1.0L, std::multiplies<long double>());

    double* d_M;
    double* d_r;
    double* d_pm;
    size_t matrixSize = n * nc * sizeof(double);
    size_t rowSize = (1 << n) * nc * sizeof(double);

    cudaMalloc(&d_M, matrixSize);
    cudaMalloc(&d_r, rowSize);
    cudaMalloc(&d_pm, sizeof(double));

    cudaMemcpy(d_M, M.data()->data(), matrixSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, r.data(), nc * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_pm, &pm, sizeof(double), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = ((1 << n) + blockSize - 1) / blockSize;

    computeGrayCodeAndRowSum<<<numBlocks, blockSize>>>(d_M, d_r, d_pm, n, nc);
    cudaDeviceSynchronize();

    cudaMemcpy(&pm, d_pm, sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_M);
    cudaFree(d_r);
    cudaFree(d_pm);

    return pm * std::pow(-1, n);
}

template<typename T>
std::vector<std::vector<T>> getMatrix(int n) {
    std::vector<std::vector<T>> matrix(n, std::vector<T>(n, 0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dist(0.0, 1.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double num = dist(gen);
            if (num < 0.6)
                matrix[i][j] = num;
            else
                matrix[i][j] = 0;
        }
    return matrix;
}

int main(int argc, char** argv) {
    int size = atoi(argv[1]);
    std::vector<std::vector<double>> matrix = getMatrix<double>(size);

    try {
        auto start = std::chrono::high_resolution_clock::now();
        auto result = perm3(matrix);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end - start;
        std::cout << "Elapsed time: " << elapsed_time.count() << " seconds" << std::endl;

        std::cout << "Result for Perm3: " << result << std::endl;
        // result = perm1(matrix);
        // std::cout << "Result1: " << result << std::endl;
    } catch (const std::exception& e) {
        // Handle any errors
        std::cout << "Error: " << e.what() << std::endl;
    }

    return 0;
}
