#include <vector>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <chrono>
#include "../crcs.cpp"

auto convertToGrayCode(int n) {
    return n ^ (n >> 1);
}
// Function to compute the permanent using the SpaRyser algorithm
double SpaRyser(const CRS& crs, const CCS& ccs) {
    int n = crs.rptrs.size() - 1; // TODO: CHECK IF THIS IS CORRECT
    std::vector<double> x(n, 0.0);
    double p = 0;
    int nzeros = 0;

    // Initial population of x vector
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int ptr = crs.rptrs[i]; ptr < crs.rptrs[i + 1]; ++ptr) {
            sum += crs.rvals[ptr];
        }
        x[i] = crs.rvals[crs.rptrs[i + 1] - 1] - (sum / 2);
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

    // Main loop
    for (int g = 0; g < (1 << (n - 1)); ++g) {
        auto gray_i = convertToGrayCode(g);
        auto gray_prev = convertToGrayCode(g - 1);

        auto j = std::log2(gray_i ^ gray_prev);
        auto s = (gray_i >> static_cast<int>(j)) & 1 ? 1 : -1;

        // #pragma omp parallel for num_threads(16)
        for (int ptr = ccs.cptrs[j]; ptr < ccs.cptrs[j + 1]; ++ptr) {
            int row = ccs.rows[ptr];
            double val = ccs.cvals[ptr];
            if (x[row] == 0) {
                --nzeros;
            }
            x[row] += s * val;
            if (x[row] == 0) {
                ++nzeros;
            }
        }

        if (nzeros == 0) {
            double prod = 1;
            for (int i = 0; i < n; ++i) {
                prod *= x[i];
            }
            p += std::pow(-1, g) * prod;
        }
    }

    return p * ((4 * (n % 2)) - 2);
}

// Main function or other appropriate place in your code
int main() {
    // Define a matrix
    // std::vector<std::vector<double>> matrix = {
    //     {1, 2, 3},
    //     {4, 5, 6},
    //     {7, 8, 9}
    // };

    std::vector<std::vector<double>> matrix = {
        {7, 9, 18},
        {41, 53, 16},
        {76, 82, 9}
    };

    // Create CRS and CCS objects
    CRS crs(matrix);
    CCS ccs(matrix);

    // Output CRS representation
    std::cout << "CRS Representation:" << std::endl;
    std::cout << "rptrs: ";
    for (auto val : crs.rptrs) std::cout << val << " ";
    std::cout << std::endl << "columns: ";
    for (auto val : crs.columns) std::cout << val << " ";
    std::cout << std::endl << "rvals: ";
    for (auto val : crs.rvals) std::cout << val << " ";
    std::cout << std::endl;

    // Output CCS representation
    std::cout << "CCS Representation:" << std::endl;
    std::cout << "cptrs: ";
    for (auto val : ccs.cptrs) std::cout << val << " ";
    std::cout << std::endl << "rows: ";
    for (auto val : ccs.rows) std::cout << val << " ";
    std::cout << std::endl << "cvals: ";
    for (auto val : ccs.cvals) std::cout << val << " ";
    std::cout << std::endl;

    // Compute the permanent using the SpaRyser algorithm
    auto start = std::chrono::high_resolution_clock::now();
    double permanent = SpaRyser(crs, ccs);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Permanent of the matrix is: " << permanent << std::endl;
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << "ns" << std::endl;

    return 0;
}