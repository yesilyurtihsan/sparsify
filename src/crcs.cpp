#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

typedef int dtype;

template<typename T>
int numNonZeros(T ** A, int m, int n){
    int cnt = 0;
    for(int i=0; i<m; i++)
        for(int j=0; j<n; j++)
            if (A[i][j] != 0) cnt++;
    return cnt; 
}

// Function to return the Compressed Row Storage (CRS) of matrix A
void crs(int** A, int m, int n, int*& val, int*& col_ind, int*& row_ptr) {
    int nnz_count = 0;
    row_ptr = new int[m + 1];
    row_ptr[0] = 0;
    
    for (int i = 0; i < m; i++) {
        
        bool row_started = false;
        for (int j = 0; j < n; j++) {
            if (A[i][j] != 0) {
                if (!row_started) {
                    row_ptr[i] = nnz_count;
                    row_started = true;
                }
                val[nnz_count] = A[i][j];
                col_ind[nnz_count] = j;
                nnz_count++;
            }
            if (!row_started)
                row_ptr[j + 1] = nnz_count;
        }
    }
}

// Function to return the Compressed Column Storage (CCS) of matrix A
void ccs(int** A, int m, int n, int*& val, int*& row_ind, int*& col_ptr) {
    int nnz_count = 0;
    col_ptr = new int[n + 1];
    col_ptr[0] = 0;
    
    for (int j = 0; j < n; j++) {
        bool col_started = false;
        for (int i = 0; i < m; i++) {
            if (A[i][j] != 0) {
                if (!col_started) {
                    col_ptr[j] = nnz_count;
                    col_started = true;
                }
                val[nnz_count] = A[i][j];
                row_ind[nnz_count] = i;
                nnz_count++;
            }
        }
        if (!col_started)
            col_ptr[j + 1] = nnz_count;
    }
}

// int main() {
//     srand(time(0));
//     const int m = 3;
//     const int n = 3;
//     int** A = new int*[m];
//     for (int i = 0; i < m; i++) {
//         A[i] = new int[n];
//     }

//     for (int i = 0; i < m; i++)
//         for (int j = 0; j < n; j++){
//             int num = rand() % 100;
//             if (num < 70) A[i][j] = 0;
//             else A[i][j] = num;
//             }
//     // Initialize A with values
//     // Add this


//     int nnz = numNonZeros(A, m, n);

//     int* val_crs = new int[nnz];
//     int* col_ind_crs = new int[nnz];
//     int* row_ptr_crs;
//     crs(A, m, n, val_crs, col_ind_crs, row_ptr_crs);
    

//     int* val_ccs = new int[nnz];
//     int* row_ind_ccs = new int[nnz];
//     int* col_ptr_ccs;
//     ccs(A, m, n, val_ccs, row_ind_ccs, col_ptr_ccs);

//     cout << nnz << endl;

//     // Use val_crs, col_ind_crs, row_ptr_crs, val_ccs, row_ind_ccs, col_ptr_ccs for further processing

//     for (int i = 0; i < m; i++){
//         for (int j = 0; j < n; j++)
//            cout << A[i][j] << " ";
//         cout << endl;
//     }

//     cout << endl << nnz << endl << endl;

//     for (int i=0; i<nnz; i++)
//         cout << col_ind_crs[i] << " ";

//     cout << endl << endl;

//     for (int i=0; i<nnz; i++)
//         cout << val_crs[i] << " ";

//     cout << endl << endl;

//     for(int i=0; i<m; i++)
//         cout << row_ptr_crs[i]<<endl;

//     return 0;
// }
