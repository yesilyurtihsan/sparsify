#include<vector>

class CRS {
public:
    std::vector<int> rptrs;
    std::vector<int> columns;
    std::vector<double> rvals;

    CRS(const std::vector<std::vector<double>>& matrix) {
        int row_index = 0;
        for (const auto& row : matrix) {
            rptrs.push_back(rvals.size());
            for (int col_index = 0; col_index < row.size(); ++col_index) {
                if (row[col_index] != 0) {
                    columns.push_back(col_index);
                    rvals.push_back(row[col_index]);
                }
            }
            ++row_index;
        }
        rptrs.push_back(rvals.size()); // Add an extra entry to mark the end of the last row
    }
};

class CCS {
public:
    std::vector<int> cptrs;
    std::vector<int> rows;
    std::vector<double> cvals;

    CCS(const std::vector<std::vector<double>>& matrix) {
        int num_cols = matrix[0].size();
        cptrs.resize(num_cols + 1, 0);

        for (int col_index = 0; col_index < num_cols; ++col_index) {
            for (int row_index = 0; row_index < matrix.size(); ++row_index) {
                if (matrix[row_index][col_index] != 0) {
                    rows.push_back(row_index);
                    cvals.push_back(matrix[row_index][col_index]);
                }
            }
            cptrs[col_index + 1] = cvals.size();
        }
    }
};
