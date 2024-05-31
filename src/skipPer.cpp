#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <omp.h>
using namespace std;

//A structure that stores the i:row_id and j:col_id of nonzero:val
typedef struct Matrix {
    int i;
    int j;
    double val;
} Matrix;

Matrix* X;
int NumRows, NumCols;
int NumNonzeros;

typedef struct SparseMatrix {
    int* crs_ptrs;    // Row pointers for CRS, will be of size numrows+1
    int* crs_colids;  // Column indices for CRS, will be of size numnonzeros
    double* crs_values; // Non-zero values for CRS, will be of size numnonzeros

    int* ccs_ptrs;    // Column pointers for CCS, will be of size numcols+1
    int* ccs_rowids;     // Row indices for CCS, will be of size numnonzeros
    double* ccs_values; // Non-zero values for CCS, will be of size numnonzeros
} SparseMatrix;

//returns the number of nonzeros
int loadDataSparse(char* filename, Matrix** mat, int* numrows, int* numcols)
{
    FILE* myfile;
    if ((myfile = fopen(filename, "r")) == NULL) {
        printf("Error: file cannot be found\n");
        exit(1);
    }

    //read the first line of the file: should have numrows numnzs
    int n1, n2, numnzs, result;
    if ((result = fscanf(myfile, "%d %d", &n1,&numnzs)) <= 0) {
        printf("error while reading file %d\n", result);
        exit(1);
    }
    n2 = n1;

    printf("The number of rows, columns, and nonzeros are %d, %d, and %d respectively\n", n1, n2, numnzs);

    //will now read the nonzeros into the Matrix structure
    *mat = (Matrix*)malloc(numnzs * sizeof(Matrix));
    for (int i = 0; i < numnzs; i++) {
        int tempi, tempj;
        double tempval;
        if ((result = fscanf(myfile, "%d %d %lf", &tempi, &tempj, &tempval)) <= 0) {
            printf("error while reading file - %d\n", result);
            exit(1);
        }

        (*mat)[i].i = tempi;
        (*mat)[i].j = tempj;
        (*mat)[i].val = tempval;
    }

    fclose(myfile);

    *numrows = n1;
    *numcols = n2;

    return numnzs;
}

void compressedMatrixSetup(int*& rptrs, int*& colids, double*& rvalues, int*& cptrs, int*& rowids, double*& cvalues)
{
    rptrs = (int*)malloc((NumRows + 1) * sizeof(int));
    colids = (int*)malloc(NumNonzeros * sizeof(int));
    rvalues = (double*)malloc(NumNonzeros * sizeof(double));

    cptrs = (int*)malloc((NumCols + 1) * sizeof(int));
    rowids = (int*)malloc(NumNonzeros * sizeof(int));
    cvalues = (double*)malloc(NumNonzeros * sizeof(double));

    //fill the rptrs array
    memset(rptrs, 0, (NumRows + 1) * sizeof(int));
    for (int i = 0; i < NumNonzeros; i++) {
        int rowid = X[i].i;
        if (rowid < 0 || rowid >= NumRows) {
            printf("problem in X, quitting - %d\n", rowid);
            exit(1);
        }
        rptrs[rowid + 1]++;
    }

    //now we have cumulative ordering of rptrs.
    for (int i = 1; i <= NumRows; i++) {
        rptrs[i] += rptrs[i - 1];
    }
    printf("This number should be equal to the number of nonzeros %d\n", rptrs[NumRows]);

    // we set colids such that for each element, it holds the related column of that element
    for (int i = 0; i < NumNonzeros; i++) {
        int rowid = X[i].i;
        int index = rptrs[rowid];

        colids[index] = X[i].j;
        rvalues[index] = X[i].val;

        rptrs[rowid] = rptrs[rowid] + 1;
    }

    for (int i = NumRows; i > 0; i--) {
        rptrs[i] = rptrs[i - 1];
    }
    rptrs[0] = 0;
    printf("This number should be equal to the number of nonzeros %d\n", rptrs[NumRows]);

    //fill the cptrs array
    memset(cptrs, 0, (NumCols + 1) * sizeof(int));
    for (int i = 0; i < NumNonzeros; i++) {
        int colid = X[i].j;
        if (colid < 0 || colid >= NumCols) {
            printf("problem in X, quitting - %d\n", colid);
            exit(1);
        }
        cptrs[colid + 1]++;
    }

    //now we have cumulative ordering of cptrs.
    for (int i = 1; i <= NumCols; i++) {
        cptrs[i] += cptrs[i - 1];
    }
    printf("This number should be equal to the number of nonzeros %d\n", cptrs[NumCols]);

    // we set rowids such that for each element, it holds the related row of that element
    for (int i = 0; i < NumNonzeros; i++) {
        int colid = X[i].j;
        int index = cptrs[colid];

        rowids[index] = X[i].i;
        cvalues[index] = X[i].val;

        cptrs[colid] = cptrs[colid] + 1;
    }

    for (int i = NumCols; i > 0; i--) {
        cptrs[i] = cptrs[i - 1];
    }
    cptrs[0] = 0;
    printf("This number should be equal to the number of nonzeros %d\n", cptrs[NumCols]);
}


//int next(int g, int j) {
//    if (g < (1 << j)) {
//        return 1 << j;
//    } else {
//        return g + (1 << (j + 1)) - ((g - (1 << j)) % (1 << (j + 1)));
//    }
//}
//
//double SkipPer(const SparseMatrix& sparseMatrix) {
//    int n = NumRows;
//    double* x=(double*)malloc(NumRows*sizeof(double));
//    double p = 1.0;
//
//    // Initialize x and p
//    for (int i = 0; i < n; i++) {
//        double sum = 0.0;
//        for (int ptr = sparseMatrix.crs_ptrs[i]; ptr < sparseMatrix.crs_ptrs[i + 1]; ptr++) {
//            sum += sparseMatrix.crs_values[ptr];
//        }
//        x[i] = sum;
//        p *= sum;
//    }
//
//    int g_prime = 0;
//    int g = 1;
//    int max_g = (1 << (n - 1)) - 1;
//
//    while (g < max_g) {
//        int gr_diff = g ^ g_prime;
//
//        for (int j = 0; j < n; j++) {
//            if ((gr_diff & (1 << j)) != 0) {
//                gr_diff &= ~(1 << j);
//                int s = (g & (1 << j)) ? 1 : 0;
//                s = 2 * s - 1;
//
//                for (int ptr = sparseMatrix.ccs_ptrs[j]; ptr < sparseMatrix.ccs_ptrs[j + 1]; ptr++) {
//                    x[sparseMatrix.ccs_rowids[ptr]] += s * sparseMatrix.ccs_values[ptr];
//                }
//            }
//        }
//
//        double prod = 1.0;
//        for (int i = 0; i < n; i++) {
//            prod *= x[i];
//        }
//        p += pow(-1, g) * prod;
//
//        g_prime = g;
//        g = next(g, __builtin_ctz(g_prime));
//    }
//    free(x);
//    return p * (4 * (n % 2) - 2);
//}


double parallel_skipper(int*& rptrs, int*& cols, double*& rvals, int*& cptrs, int*& rows, double*& cvals, int threads) {
    int nov = NumRows;
    double* x = (double*)malloc(NumRows * sizeof(double));
    double rs, p;
    int j, ptr;
    unsigned long long ci, start, end, chunk_size, change_j;

    //initialize the vector entries                                                                                        
    for (j = 0; j < nov; j++) {
        rs = .0f;
        for (ptr = rptrs[j]; ptr < rptrs[j + 1]; ptr++)
            rs += rvals[ptr];
        x[j] = -rs / (2.0f);
    }

    for (ptr = cptrs[nov - 1]; ptr < cptrs[nov]; ptr++) {
        x[rows[ptr]] += cvals[ptr];
    }

    //update perman with initial x
    double prod = 1;
    for (j = 0; j < nov; j++) {
        prod *= x[j];
    }
    p = prod;

    //find start location        
    start = 1;
    for (int j = 0; j < nov; j++) {
        if (x[j] == 0) {
            change_j = -1;
            for (ptr = rptrs[j]; ptr < rptrs[j + 1]; ptr++) {
                ci = 1ULL << cols[ptr];
                if (ci < change_j) {
                    change_j = ci;
                }
            }
            if (change_j > start) {
                start = change_j;
            }
        }
    }

    end = (1ULL << (nov - 1));

    chunk_size = (end - start + 1) / threads + 1;

#pragma omp parallel num_threads(threads) private(j, ci, change_j) 
    {
        double my_x[64];
        memcpy(my_x, x, sizeof(double) * 64);

        int tid = omp_get_thread_num();
        unsigned long long my_start = start + tid * chunk_size;
        unsigned long long my_end = min(start + ((tid + 1) * chunk_size), end);

        //update if neccessary
        double my_p = 0;

        unsigned long long my_gray;
        unsigned long long my_prev_gray = 0;

        int ptr, last_zero;
        unsigned long long period, steps, step_start;

        unsigned long long i = my_start;

        while (i < my_end) {
            //k = __builtin_ctzll(i + 1);
            my_gray = i ^ (i >> 1);

            unsigned long long gray_diff = my_prev_gray ^ my_gray;
            //cout << "X: " << gray_diff << endl;
            j = 0;
            while (gray_diff > 0) { // this contains the bit to be updated
                unsigned long long onej = 1ULL << j;
                if (gray_diff & onej) { // if bit l is changed 
                    gray_diff ^= onej;   // unset bit
                    if (my_gray & onej) {    // do the update
                        for (ptr = cptrs[j]; ptr < cptrs[j + 1]; ptr++) {
                            my_x[rows[ptr]] += cvals[ptr];
                        }
                    }
                    else {
                        for (ptr = cptrs[j]; ptr < cptrs[j + 1]; ptr++) {
                            my_x[rows[ptr]] -= cvals[ptr];
                        }
                    }
                }
                j++;
            }
            //counter++;
            my_prev_gray = my_gray;
            last_zero = -1;
            double my_prod = 1;
            for (j = nov - 1; j >= 0; j--) {
                my_prod *= my_x[j];
                if (my_x[j] == 0) {
                    last_zero = j;
                    break;
                }
            }

            if (my_prod != 0) {
                my_p += ((i & 1ULL) ? -1.0 : 1.0) * my_prod;
                i++;
            }
            else {
                change_j = -1;
                for (ptr = rptrs[last_zero]; ptr < rptrs[last_zero + 1]; ptr++) {
                    step_start = 1ULL << cols[ptr];
                    period = step_start << 1;
                    ci = step_start;
                    if (i >= step_start) {
                        steps = (i - step_start) / period;
                        ci = step_start + ((steps + 1) * period);
                    }
                    if (ci < change_j) {
                        change_j = ci;
                    }
                }

                i++;
                if (change_j > i) {
                    i = change_j;
                }
            }
        }

#pragma omp critical
        p += my_p;
    }
    return ((4 * (nov & 1) - 2) * p);
}

int main(int argc, char* argv[]) {

    char* fileName;
    if (argc == 2)
        fileName = argv[1];  // Replace with your file name
    else {
        char tempfname[16] = "ey35_02.mat";
        //char tempfname[16] = "sample.mat";
        fileName = tempfname;
    }

    NumNonzeros = loadDataSparse(fileName, &X, &NumRows, &NumCols);

    printf("Matrix is read. There are %d rows, %d columns, and %d nonzeros.\n", NumRows, NumCols, NumNonzeros);

    SparseMatrix mat;
    compressedMatrixSetup(mat.crs_ptrs, mat.crs_colids, mat.crs_values, mat.ccs_ptrs, mat.ccs_rowids, mat.ccs_values);

    double start_time = omp_get_wtime();
    double permanent_slow = parallel_skipper(mat.crs_ptrs, mat.crs_colids, mat.crs_values, mat.ccs_ptrs, mat.ccs_rowids, mat.ccs_values, 1);
    double end_time = omp_get_wtime();

    cout << "Permanent with 1 thread: " << permanent_slow << endl << "Time: " << end_time - start_time << endl;

    start_time = omp_get_wtime();
    double permanent_fast = parallel_skipper(mat.crs_ptrs, mat.crs_colids, mat.crs_values, mat.ccs_ptrs, mat.ccs_rowids, mat.ccs_values, 16);
    end_time = omp_get_wtime();

    cout << "Permanent with 16 threads: " << permanent_fast << endl << "Time: " << end_time - start_time << endl;

    free(X);
    free(mat.crs_ptrs);
    free(mat.crs_colids);
    free(mat.crs_values);
    free(mat.ccs_ptrs);
    free(mat.ccs_rowids);
    free(mat.ccs_values);

    return 0;
}
