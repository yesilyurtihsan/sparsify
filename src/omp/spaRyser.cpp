#include <vector>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <chrono>
#include "../crcs.h"

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
