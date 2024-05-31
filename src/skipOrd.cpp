#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
using namespace std;

vector<vector<double>> sampleMat =//( 20, vector<double>(20, 0) );
{
    // 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9
      {0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0}, //0
      {1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1}, //1
      {0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0}, //2
      {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1}, //3
      {1,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0}, //4

      {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, //5
      {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0}, //6
      {0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0}, //7
      {0,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0}, //8
      {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0}, //9

      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},//10
      {0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1},//11
      {0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},//12
      {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},//13
      {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},//14

      {1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1},//15
      {0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0},//16
      {1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1},//17
      {1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0},//18
      {0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0},//19

}; //54 non-zeros in a 20x20 matrix, so 0.135 density

// Function to generate a sparse matrix with random values
vector<vector<double>> generateSparseMatrix(int n, double density) {
    vector<vector<double>> mat(n, vector<double>(n, 0.0)); // Initialize matrix with zeros

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(0.0, 1.0); // Uniform distribution between 0 and 1
    uniform_int_distribution<int> expDis(0, 5); // to randomly multiply the result from dis by some power of 10

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (dis(gen) <= density) {
                double value = dis(gen);
                int power = expDis(gen);
                mat[i][j] = value * pow(10, power); // Assign a random value between 0 and 1
            }
        }
    }

    return mat;
}

void printVectorMatrix(vector<vector<double>> mat,double density) {
    int n = mat.size(); //it's an nxn matrix
    cout << "Sparse Matrix of size " << n << "x" << n<< " with density " << density << ":\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            cout << mat[i][j] << "  ";
        cout << endl;
    }

}

void numNonzeros(vector<vector<double>> mat) {

    int nzcount = 0;
    for (int i = 0; i < mat.size(); i++)
        for (int j = 0; j < mat[i].size(); j++)
            if (mat[i][j] != 0)
                nzcount++;
    cout << nzcount << " nonzeros in the matrix\n";
}

// Function to permute rows and columns of a matrix
vector<vector<double>> permuteMatrix(const vector<vector<double>>& matrix, const vector<int>& rowPerm, const vector<int>& colPerm) {
    int numRows = matrix.size();
    int numCols = matrix[0].size();

    vector<vector<double>> permutedMatrix(numRows, vector<double>(numCols));

    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            permutedMatrix[rowPerm[i]][colPerm[j]] = matrix[i][j];
        }
    }

    return permutedMatrix;
}

void testingPermutes()
{
    vector<vector<double>> mat = {
        {1,2,3,4,5},
        {6,7,8,9,10},
        {11,12,13,14,15},
        {16,17,18,19,20},
        {21,22,23,24,25}
    };

    vector<vector<double>> matPermuted = permuteMatrix(mat,
        { 0,1,2,3,4 },
        { 0,1,2,3,4 });
    printVectorMatrix(matPermuted,1);
    matPermuted = permuteMatrix(mat,
        { 0,1,2,3,4 },
        { 1,0,2,3,4 });
    printVectorMatrix(matPermuted, 1);
    matPermuted = permuteMatrix(mat,
        { 0,1,2,4,3 },
        { 0,1,2,3,4 });
    printVectorMatrix(matPermuted, 1);
    matPermuted = permuteMatrix(mat,
        { 1,2,3,4,0 },
        { 1,2,3,4,0 });
    printVectorMatrix(matPermuted, 1);
}

void printVector(vector<int> vec)
{
    for (int i : vec)
        cout << i << " ";
    cout << endl;
}

// Function to perform the SkipOrd algorithm
void skipOrd(const vector<vector<double>>& matrix, vector<int>& rowPerm, vector<int>& colPerm) {
    int n = matrix.size();
    vector<bool> rowVisited(n, false);
    vector<int> degs(n, 0);

    // Calculate the degree of each column
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            if (matrix[i][j] != 0) {
                degs[j]++;
            }
        }
    }

    // Main loop of the SkipOrd algorithm
    for (int j = 0; j < n; ++j) {
        // Find the column with the minimum degree among the remaining columns
        int curCol = min_element(degs.begin(), degs.end()) - degs.begin();
        degs[curCol] = n; // Mark the column as visited by setting its degree to n
        colPerm.push_back(curCol);

        // Process rows with non-zero entries in the selected column
        for (int i = 0; i < n; ++i) {
            if (matrix[i][curCol] != 0 && !rowVisited[i]) {
                rowVisited[i] = true;
                rowPerm.push_back(i);

                // Update the degrees of columns connected to the current row
                for (int k = 0; k < n; ++k) {
                    if (matrix[i][k] != 0 && degs[k] != n) {
                        degs[k]--;
                    }
                }
            }
        }
    }

    //push back unvisited rows
    for (int i = 0; i < n; i++)
        if (!rowVisited[i])
            rowPerm.push_back(i);
}

//int main() {
//    int n = 50; // Set the size of the matrix (20, 50, 100, or 1000)
//    double density = 0.01; // Set the density of non-zero elements (adjust as needed)
//
//    vector<vector<double>> sparseMat = generateSparseMatrix(n, density);
//
//    printVectorMatrix(sparseMat,0.01);
//    numNonzeros(sparseMat);
//    
//    vector<int> rowPerm, colPerm;
//    skipOrd(sparseMat, rowPerm, colPerm);
//    printVector(colPerm); printVector(rowPerm);
//
//    vector<vector<double>> sparseMatPermuted = permuteMatrix(sparseMat, rowPerm, colPerm);
//    printVectorMatrix(sparseMatPermuted,0.01);
//
//
//    return 0;
//}




