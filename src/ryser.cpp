#include<iostream>
#include<cmath> // This might be removed

// TODO: Add this functionality
unsigned int gray_g(int g, int n){
    return 0x0000;
}

template<typename T>
T ryser(T** A, int n){
    T* x = malloc(n*sizeof(T));
    for(int i=0; i<n; i++){
        T sum = 0;
        for(int j=0; j<n; j++)
            sum = sum + A[i][j];
        x[i] = A[i][n] - sum/2;
    }

    T p=1;
    for (int i=0; i<n;)
        p *= x[i++];

    for(int g=1; g<(pow(2, n-1)-1); g++){
        unsigned int j = log2(gray_g(g, n) ^ gray_g(g-1, n));
        unsigned int s = 2*gray_g(j, n) - 1;
        T prod = 1;
        for(int i=0; i<n; i++){
            x[i] = x[i] + s*A[i][j];
            prod *= x[i];
        }
        p += pow(-1,g)* prod;
    }
    
    return p * (4 * (n % 2) - 2);
}


int main(){

    return 0;
}