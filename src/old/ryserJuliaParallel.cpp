#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <numeric>
#include <omp.h>
#include <chrono>



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

    // for (int i=0; i < n; ++i) {     
    //     lut[1LL << i] = {i + 1, 1}; // Mapping for positive sign
    //     lut[-(1LL << i)] = {i + 1, -1}; // Mapping for negative sign
    // }

    #pragma omp parallel num_threads(n)
    {
        int ind = omp_get_thread_num();
        #pragma omp critical
        {
            lut[1LL << ind] = {ind + 1, 1}; // Mapping for positive sign
            lut[-(1LL << ind)] = {ind + 1, -1}; // Mapping for negative sign
        }
    }
    

    // Initialize row sum vector with the first row of the matrix
    std::vector<int> r = M[0];
    int s = -1; // Sign variable
    // Calculate the initial product of row sums
    long long pm = s * std::accumulate(r.begin(), r.end(), 1LL, std::multiplies<long long>());
    
    // Iterate over all subsets of rows using gray code

    for (int j=1; j < n; j++){
        for (int i = (1 << j); i < (1 << (j+1)); i++){
            if (i % 2 == 0) {
                gd[i] = 2 * gd[i / 2];
            } else {
                gd[i] = std::pow(-1, (i - 1) / 2);
            }
        }

        // for (int i = (1 << j); i < (1 << (j+1)); i++){
        //     // Use the lookup table to find the corresponding row and sign
            
        //         auto [idx, sign] = lut[gd[i]];
        //         // Update row sums
        //         for (int k = 0; k < nc; ++k) {
        //             r[k] += M[idx - 1][k] * sign;
        //         }

        //     // Update sign for the next iteration
        //     s *= -1;
        //     // Add the product of updated row sums to the permanent
        //     #pragma omp atomic
        //     pm += s * std::accumulate(r.begin(), r.end(), 1LL, std::multiplies<long long>());
        // }

        // int n_threads = pow(2, j);
        // #pragma omp parallel num_threads(n_threads)
        // {
        //     int i = omp_get_thread_num()+n_threads;
        //     // std::cout << omp_get_thread_num() << std::endl;
        //     auto [idx, sign] = lut[gd[i]];
        //     #pragma omp critical
        //     std::cout << idx << " " << sign << std::endl;
        //     // Update row sums
        //     #pragma omp critical
        //     {
        //         for (int k = 0; k < nc; ++k) {
        //         r[k] += M[idx - 1][k] * sign;
        //         }
        //         // Update sign for the next iteration
        //         s *= -1;
        //         // Add the product of updated row sums to the permanent
        //         pm += s * std::accumulate(r.begin(), r.end(), 1LL, std::multiplies<long long>());
        //     }
        // }

        for (int i = (1 << j); i < (1 << (j+1)); i++){
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

int main(int argc, char** argv) {
    // Example matrix
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

        std::cout << "Result: " << result << std::endl;
    } catch (const std::exception& e) {
        // Handle any errors
        std::cout << "Error: " << e.what() << std::endl;
    }

    return 0;
}