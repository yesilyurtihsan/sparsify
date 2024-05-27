#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <bitset>
#include <chrono>
#include <omp.h>
#include <random>

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

class GrayCodeSign {
private:
    int n;
    int* grayCodeValues;
    int* signValues;

public:
    GrayCodeSign(int size) {
        n = size;
        grayCodeValues = new int[n];
        signValues = new int[n];
        grayCodeValues[0] = 0;
        grayCodeValues[1] = 1;
        signValues[0] = 0;
        signValues[1] = 1;
    }

    void addValue(int idx, int grayCode, int sign) {
        grayCodeValues[idx] = grayCode;
        signValues[idx] = sign;
    }

    int getGrayCodeValue(int idx) {
        return grayCodeValues[idx];
    }

    int getSignValue(int idx) {
        return signValues[idx];
    }

    int size() {
        return n;
    }

    ~GrayCodeSign() {
        delete[] grayCodeValues;
        delete[] signValues;
    }
};



// Function to calculate the permanent of a square matrix using Ryser's formula.
// This is an efficient algorithm for computing the permanent of a matrix, which is a concept similar to the determinant.
long double perm3(const std::vector<std::vector<double>>& M) {
    int n = M.size(); // Number of rows in the matrix
    int nc = M[0].size(); // Number of columns in the matrix
    
    // Ensure the matrix is square
    if (n != nc) {
        throw std::runtime_error("Matrix must be square");
    }
    
    // Initialize row sum vector with the first row of the matrix
    std::vector<double> r = M[0];
    double s = -1; // Sign variable
    // Calculate the initial product of row sums

    auto pm = s * std::accumulate(r.begin(), r.end(), 1.0L, std::multiplies<long double>());
    
    #pragma omp parallel for //reduction(+: pm)
    for (int i = 2; i < (1 << n); i++){
        auto gray_i = convertToGrayCode(i);
        auto gray_prev = convertToGrayCode(i - 1);
        
        auto idx = std::log2(gray_i ^ gray_prev);
        auto sign = (gray_i >> static_cast<int>(idx)) & 1 ? 1 : -1;

        // #pragma omp parallel for
        for (int k = 0; k < nc; k++)
            r[k] += M[idx][k] * sign;

        s *= -1;
        // Add the product of updated row sums to the permanent
        // #pragma omp critical
        pm += s * std::accumulate(r.begin(), r.end(), 1.0L, std::multiplies<long double>());
    }

    // Final adjustment of the sign based on the number of rows
    return pm * std::pow(-1, n);
}


template<typename T>
std::vector<std::vector<T>> getMatrix(int n){
    std::vector<std::vector<T>> matrix(n, std::vector<T>(n, 0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dist(0.0, 1.0);
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++){
            double num = dist(gen);
            if (num < 0.6)
                matrix[i][j] = num;
            else
                matrix[i][j] = 0;
        }
    return matrix;
}


int main(int argc, char ** argv) {
    // std::vector<std::vector<double>> matrix = {
    //     {1, 2, 3},
    //     {4, 5, 6},
    //     {7, 8, 9}
    // };
    int size = atoi(argv[1]);
    std::vector<std::vector<double>> matrix = getMatrix<double>(size);
   
    // for (int i=0; i<size; i++){
    //     for(int j=0; j<size; j++)
    //         std::cout << matrix[i][j] << " ";
    //     std::cout << std::endl;
    // }

    std::cout << std::endl;
  

    try {
        // Calculate and display the permanent of the matrix
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