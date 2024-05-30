#include <vector>
#include <cmath>
#include <iostream>
#include <chrono>
#include <cuda_runtime.h>
#include "../crcs.h"

__device__ double atomicAddExp(double* address, double val)
{
    unsigned long long int* address_as_ull =
                             (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

__device__ int convertToGrayCode(int n) {
    return n ^ (n >> 1);
}

__global__ void initialPopulation(const int* d_rptrs, const double* d_rvals, double* d_x, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        double sum = 0;
        for (int ptr = d_rptrs[i]; ptr < d_rptrs[i + 1]; ++ptr) {
            sum += d_rvals[ptr];
        }
        d_x[i] = d_rvals[d_rptrs[i + 1] - 1] - (sum / 2);
    }
}

__global__ void mainLoop(const int* d_cptrs, const int* d_rows, const double* d_cvals, double* d_x, int n, double* p, long int g, int* d_nzeros) {
    long int gray_i = convertToGrayCode(g);
    long int gray_prev = convertToGrayCode(g - 1);
    long int j = __ffs(gray_i ^ gray_prev) - 1;
    int s = ((gray_i >> j) & 1) ? 1 : -1;

    for (int ptr = d_cptrs[j]; ptr < d_cptrs[j + 1]; ++ptr) {
        int row = d_rows[ptr];
        double val = d_cvals[ptr];
        if (d_x[row] == 0) {
            atomicAdd(d_nzeros, -1);
        }
        atomicAddExp(&d_x[row], s * val);
        if (d_x[row] == 0) {
            atomicAdd(d_nzeros, 1);
        }
    }

    if (*d_nzeros == 0) {
        double prod = 1;
        for (int i = 0; i < n; ++i) {
            prod *= d_x[i];
        }
        atomicAddExp(p, pow(-1, g) * prod);
    }
}

__global__ void mainLoop(const int* d_cptrs, const int* d_rows, const double* d_cvals, float* d_x, int n, float* p, long int g, int* d_nzeros) {
    long int gray_i = convertToGrayCode(g);
    long int gray_prev = convertToGrayCode(g - 1);
    int j = __ffs(gray_i ^ gray_prev) - 1;
    int s = ((gray_i >> j) & 1) ? 1 : -1;

    for (int ptr = d_cptrs[j]; ptr < d_cptrs[j + 1]; ++ptr) {
        int row = d_rows[ptr];
        double val = d_cvals[ptr];
        if (d_x[row] == 0) {
            atomicAdd(d_nzeros, -1);
        }
        atomicAdd(&d_x[row], s * val);
        if (d_x[row] == 0) {
            atomicAdd(d_nzeros, 1);
        }
    }

    if (*d_nzeros == 0) {
        double prod = 1;
        for (int i = 0; i < n; ++i) {
            prod *= d_x[i];
        }
        atomicAdd(p, pow(-1, g) * prod);
    }
}

template<typename T>
T SpaRyser(const CRS& crs, const CCS& ccs) {
    int n = crs.rptrs.size() - 1;

    // Allocate device memory
    int* d_rptrs, * d_columns, * d_cptrs, * d_rows, * d_nzeros;
    double* d_rvals, * d_cvals, * d_x, * d_p;
    cudaMalloc((void**)&d_rptrs, crs.rptrs.size() * sizeof(int));
    cudaMalloc((void**)&d_columns, crs.columns.size() * sizeof(int));
    cudaMalloc((void**)&d_rvals, crs.rvals.size() * sizeof(double));
    cudaMalloc((void**)&d_cptrs, ccs.cptrs.size() * sizeof(int));
    cudaMalloc((void**)&d_rows, ccs.rows.size() * sizeof(int));
    cudaMalloc((void**)&d_cvals, ccs.cvals.size() * sizeof(double));
    cudaMalloc((void**)&d_x, n * sizeof(T));
    cudaMalloc((void**)&d_p, sizeof(T));
    cudaMalloc((void**)&d_nzeros, sizeof(int));

    // Copy data to device
    cudaMemcpy(d_rptrs, crs.rptrs.data(), crs.rptrs.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_columns, crs.columns.data(), crs.columns.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rvals, crs.rvals.data(), crs.rvals.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cptrs, ccs.cptrs.data(), ccs.cptrs.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rows, ccs.rows.data(), ccs.rows.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cvals, ccs.cvals.data(), ccs.cvals.size() * sizeof(double), cudaMemcpyHostToDevice);

    // Initialize x vector
    int threadsPerBlock = 256;
    int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
    initialPopulation<<<blocksPerGrid, threadsPerBlock>>>(d_rptrs, d_rvals, d_x, n);

    // Copy x back to host to check zeros
    std::vector<double> x(n);
    cudaMemcpy(x.data(), d_x, n * sizeof(T), cudaMemcpyDeviceToHost);

    T p = 0;
    int nzeros = 0;
    for (int i = 0; i < n; ++i) {
        if (x[i] == 0) {
            ++nzeros;
        }
    }

    if (nzeros > 0) {
        p = 1;
        for (int i = 0; i < n; ++i) {
            p *= x[i];
        }
    } else {
        p = 0;
    }

    cudaMemcpy(d_p, &p, sizeof(T), cudaMemcpyHostToDevice);
    cudaMemcpy(d_nzeros, &nzeros, sizeof(int), cudaMemcpyHostToDevice);

    // Main loop
    blocksPerGrid = 32;
    threadsPerBlock = 2048;
    for (long int g = 0; g < pow(2, n-1) / (blocksPerGrid * threadsPerBlock); ++g) {                        // HANDLE THE FOR LOOP AND THREAD INDEXING
        if (g % 1000000 == 0) std::cout << g << std::endl;
        mainLoop<<<blocksPerGrid, threadsPerBlock>>>(d_cptrs, d_rows, d_cvals, d_x, n, d_p, g, d_nzeros);
    }

    // Copy result back to host
    cudaMemcpy(&p, d_p, sizeof(T), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_rptrs);
    cudaFree(d_columns);
    cudaFree(d_rvals);
    cudaFree(d_cptrs);
    cudaFree(d_rows);
    cudaFree(d_cvals);
    cudaFree(d_x);
    cudaFree(d_p);
    cudaFree(d_nzeros);

    return p * ((4 * (n % 2)) - 2);
}