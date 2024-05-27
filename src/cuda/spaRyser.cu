#include <vector>
#include <cmath>
#include <iostream>
#include <chrono>
#include <cuda_runtime.h>

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

class CRS {
public:
    std::vector<int> rptrs;
    std::vector<int> columns;
    std::vector<double> rvals;

    CRS(const std::vector<std::vector<double>>& matrix) {
        int row_index = 0;
        for (const auto& row : matrix) {
            rptrs.push_back(rvals.size());
            for (int col_index = 0; col_index < row.size(); ++col_index) {
                if (row[col_index] != 0) {
                    columns.push_back(col_index);
                    rvals.push_back(row[col_index]);
                }
            }
            ++row_index;
        }
        rptrs.push_back(rvals.size());
    }
};

class CCS {
public:
    std::vector<int> cptrs;
    std::vector<int> rows;
    std::vector<double> cvals;

    CCS(const std::vector<std::vector<double>>& matrix) {
        int num_cols = matrix[0].size();
        cptrs.resize(num_cols + 1, 0);

        for (int col_index = 0; col_index < num_cols; ++col_index) {
            for (int row_index = 0; row_index < matrix.size(); ++row_index) {
                if (matrix[row_index][col_index] != 0) {
                    rows.push_back(row_index);
                    cvals.push_back(matrix[row_index][col_index]);
                }
            }
            cptrs[col_index + 1] = cvals.size();
        }
    }
};

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

__global__ void mainLoop(const int* d_cptrs, const int* d_rows, const double* d_cvals, double* d_x, int n, double* p, int g, int* d_nzeros) {
    int gray_i = convertToGrayCode(g);
    int gray_prev = convertToGrayCode(g - 1);
    int j = __ffs(gray_i ^ gray_prev) - 1;
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

double SpaRyser(const CRS& crs, const CCS& ccs) {
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
    cudaMalloc((void**)&d_x, n * sizeof(double));
    cudaMalloc((void**)&d_p, sizeof(double));
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
    cudaMemcpy(x.data(), d_x, n * sizeof(double), cudaMemcpyDeviceToHost);

    double p = 0;
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

    cudaMemcpy(d_p, &p, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_nzeros, &nzeros, sizeof(int), cudaMemcpyHostToDevice);

    // Main loop
    for (int g = 0; g < (1 << (n - 1)); ++g) {
        mainLoop<<<1, 1>>>(d_cptrs, d_rows, d_cvals, d_x, n, d_p, g, d_nzeros);
    }

    // Copy result back to host
    cudaMemcpy(&p, d_p, sizeof(double), cudaMemcpyDeviceToHost);

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

int main() {
    std::vector<std::vector<double>> matrix = {
        {7, 9, 18},
        {41, 53, 16},
        {76, 82, 9}
    };

    CRS crs(matrix);
    CCS ccs(matrix);

    std::cout << "CRS Representation:" << std::endl;
    std::cout << "rptrs: ";
    for (auto val : crs.rptrs) std::cout << val << " ";
    std::cout << std::endl << "columns: ";
    for (auto val : crs.columns) std::cout << val << " ";
    std::cout << std::endl << "rvals: ";
    for (auto val : crs.rvals) std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "CCS Representation:" << std::endl;
    std::cout << "cptrs: ";
    for (auto val : ccs.cptrs) std::cout << val << " ";
    std::cout << std::endl << "rows: ";
    for (auto val : ccs.rows) std::cout << val << " ";
    std::cout << std::endl << "cvals: ";
    for (auto val : ccs.cvals) std::cout << val << " ";
    std::cout << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    double permanent = SpaRyser(crs, ccs);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Permanent of the matrix is: " << permanent << std::endl;
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << "ns" << std::endl;

    return 0;
}
