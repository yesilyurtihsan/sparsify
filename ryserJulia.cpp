#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <numeric>

// Function to calculate the permanent of a square matrix using Ryser's formula.
// This is an efficient algorithm for computing the permanent of a matrix, which is a concept similar to the determinant.
long long perm3(const std::vector<std::vector<int>>& M) {
    int n = M.size(); // Number of rows in the matrix
    int nc = M[0].size(); // Number of columns in the matrix
    
    // Ensure the matrix is square
    if (n != nc) {
        throw std::runtime_error("Matrix must be square");
    }
    
    // Gray code sequence initialization
    std::vector<long long> gd(1 << n, 0);
    gd[1] = 1; // Base case for gray code

    // Lookup table to map gray code to matrix row and sign
    std::unordered_map<long long, std::pair<int, int>> lut;
    for (int i = 0; i < n; ++i) {
        lut[1LL << i] = {i + 1, 1}; // Mapping for positive sign
        lut[-(1LL << i)] = {i + 1, -1}; // Mapping for negative sign
    }

    // Initialize row sum vector with the first row of the matrix
    std::vector<int> r = M[0];
    int s = -1; // Sign variable
    // Calculate the initial product of row sums
    long long pm = s * std::accumulate(r.begin(), r.end(), 1LL, std::multiplies<long long>());
    
    // Iterate over all subsets of rows using gray code
    for (int i = 2; i < (1 << n); ++i) {
        // Update gray code sequence
        if (i % 2 == 0) {
            gd[i] = 2 * gd[i / 2];
        } else {
            gd[i] = std::pow(-1, (i - 1) / 2);
        }
        
        // Use the lookup table to find the corresponding row and sign
        auto [idx, sign] = lut[gd[i]];
        // Update row sums
        for (int j = 0; j < nc; ++j) {
            r[j] += M[idx - 1][j] * sign;
        }
        // Update sign for the next iteration
        s *= -1;
        // Add the product of updated row sums to the permanent
        pm += s * std::accumulate(r.begin(), r.end(), 1LL, std::multiplies<long long>());
    }
    // Final adjustment of the sign based on the number of rows
    return pm * std::pow(-1, n);
}

int main() {
    // Example matrix
    std::vector<std::vector<int>> matrix = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    try {
        // Calculate and display the permanent of the matrix
        long long result = perm3(matrix);
        std::cout << "Result: " << result << std::endl;
    } catch (const std::exception& e) {
        // Handle any errors
        std::cout << "Error: " << e.what() << std::endl;
    }

    return 0;
}