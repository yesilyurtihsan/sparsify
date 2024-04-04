#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <bitset>
#include <chrono>

// long long perm1(const std::vector<std::vector<int>>& M) {
//     int n = M.size(); // Number of rows in the matrix
//     int nc = M[0].size(); // Number of columns in the matrix
//     long long P = 0; // Initialize permanent to 0

//     // Iterate over all subsets of rows except the empty set
//     for (int subset = 1; subset < (1 << n); ++subset) {
//         long long productSum = 1; // Initialize product of sums to 1
//         for (int col = 0; col < nc; ++col) {
//             long long colSum = 0; // Initialize column sum to 0
//             for (int row = 0; row < n; ++row) {
//                 if (subset & (1 << row)) { // Check if row is in subset
//                     colSum += M[row][col];
//                 }
//             }
//             productSum *= colSum; // Multiply by the sum of the current column
//         }
//         // Add or subtract the productSum based on the parity of the subset size
//         P += (__builtin_popcount(subset) % 2 == 0 ? -productSum : productSum);
//     }
//     // Multiply by (-1)^n
//     // return P * std::pow(-1, n); // TODO: bi hata var, bunu commentleyince calisiyor
//     // ama normalde o satirin kalmasi gerekir.bak
//     return P;
// }

auto convertToGrayCode(int n) {
    return n ^ (n >> 1);
}

// Function to calculate the permanent of a square matrix using Ryser's formula.
// This is an efficient algorithm for computing the permanent of a matrix, which is a concept similar to the determinant.
long long perm3(const std::vector<std::vector<int>>& M) {
    int n = M.size(); // Number of rows in the matrix
    int nc = M[0].size(); // Number of columns in the matrix
    
    // Ensure the matrix is square
    if (n != nc) {
        throw std::runtime_error("Matrix must be square");
    }
    
    // Initialize row sum vector with the first row of the matrix
    std::vector<int> r = M[0];
    int s = -1; // Sign variable
    // Calculate the initial product of row sums
    long long pm = s * std::accumulate(r.begin(), r.end(), 1LL, std::multiplies<long long>());
    
    // Iterate over all subsets of rows using gray code
    for (int i = 2; i < (1 << n); ++i) {
       auto gray_i = convertToGrayCode(i);
       auto gray_prev = convertToGrayCode(i - 1);
       
       auto idx = std::log2(gray_i ^ gray_prev);
       auto sign = (gray_i >> static_cast<int>(idx)) & 1 ? 1 : -1;

    //    std::cout << "gray_i: " << std::bitset<8>(gray_i) << " gray_prev: " << std::bitset<8>(gray_prev) << " idx: " << idx << " sign: " << sign << std::endl;

        // Update row sums
        for (int j = 0; j < nc; ++j) {
            r[j] += M[idx][j] * sign;
        }
        // Update sign for the next iteration
        s *= -1;
        // Add the product of updated row sums to the permanent
        pm += s * std::accumulate(r.begin(), r.end(), 1LL, std::multiplies<long long>());
    }
    // Final adjustment of the sign based on the number of rows
    return pm * std::pow(-1, n);
}

std::vector<std::vector<int>> getMatrix(int n){
    std::vector<std::vector<int>> matrix(n, std::vector<int>(n, 0));
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            matrix[i][j] = j+1;
    return matrix;
}


int main(int argc, char ** argv) {
    // std::vector<std::vector<int>> matrix = {
    //     {1, 2, 3},
    //     {4, 5, 6},
    //     {7, 8, 9}
    // };
    int size = atoi(argv[1]);
    std::vector<std::vector<int>> matrix = getMatrix(size);
   

    try {
        // Calculate and display the permanent of the matrix
        auto start = std::chrono::high_resolution_clock::now();
        long long result = perm3(matrix);
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