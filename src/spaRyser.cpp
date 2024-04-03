#include<iostream>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<unordered_map>
#include<random>
#include"crcs.cpp"

using namespace std;

typedef unsigned int dtype;

template<typename T>
T spaRyser(CRS<T>& crs, CCS<T>& ccs){       // Assumes square matrix
    int n = crs.getM();
    int nnz = 0;
    double* x = new double[n];
    for(int i=0; i<n; i++){
        double sum = 0;
        for(int ptr = crs.getRowPtr()[i]; ptr < crs.getRowPtr()[i+1]; ptr++)
            sum += crs.getVal()[ptr];
        x[i] = crs.getVal()[crs.getRowPtr()[i+1]-1] - sum/2;
        if(x[i]==0)
            nnz++;
    }
    
    double p=1;
    if (nnz > 0)
        for(int i=0; i<n; i++)
            p *= x[i];
    else p = 0;

    // Gray code sequence initialization
    long long gd[1 << n] = {0};
    gd[1] = 1; // Base case for gray code

    // Lookup table to map gray code to matrix row and sign
    std::unordered_map<long long, std::pair<int, int>> lut;
    for (int i = 0; i < n; ++i) {
        lut[1LL << i] = {i + 1, 1}; // Mapping for positive sign
        lut[-(1LL << i)] = {i + 1, -1}; // Mapping for negative sign
    }

    for(int g = 1; g < (1 << n); g++){
        if (g % 2 == 0) {
            gd[g] = 2 * gd[g / 2];
        } else {
            gd[g] = std::pow(-1, (g - 1) / 2);
        }
        
        // Use the lookup table to find the corresponding row and sign
        auto [idx, sign] = lut[gd[g]];
        // Update row sums
        for(int ptr = ccs.getColPtr()[idx]; ptr < ccs.getColPtr()[idx+1]-1;g++){
            auto row = ccs.getRowInd()[ptr];
            auto val = ccs.getVal()[ptr];
            if (x[row] == 0) nnz--;
            x[row] += sign*val;
            if (x[row] == 0) nnz++;
        }
        if (nnz == 0){
            double prod = 1;
            for(int i=0; i<n; i++)
                prod *= x[i];
            p += pow(-1, g) * prod;
        }
    }
    return p * (4 * (n % 2) - 2);
}


int main(int argc, char ** argv){
    int N = atoi(argv[1]);
    srand(time(0));

    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_real_distribution<dtype> dist(0.0, 1.0);

    const int m = N;
    const int n = N;
    
    dtype** A = new dtype*[m];
    for (int i = 0; i < m; i++) {
        A[i] = new dtype[n];
    }


    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++){
            int num = rand() % 100;
            if (num < 75) A[i][j] = 0;
            else A[i][j] = num;
            // double num = dist(gen);
            // if (num <0.5) A[i][j] = 0;
            // else A[i][j] = num;
            }
    
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
    
    cout << endl;

    CRS<dtype> crs(A,m,n);
    CCS<dtype> ccs(A,m,n);

    dtype perm = spaRyser(crs, ccs);
    cout << perm <<endl;

    return 0;
}