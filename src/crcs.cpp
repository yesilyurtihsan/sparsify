#include<vector>
#include"crcs.h"


CRS::CRS(const std::vector<std::vector<double>>& matrix) {
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

CRS::CRS(const std::vector<M3>& m3s, int n) {
    rptrs.resize(n + 1, 0);

    for (const auto& m : m3s) {
        rptrs[m.x + 1]++;
    }

    for (int i = 1; i <= n; ++i) {
        rptrs[i] += rptrs[i - 1];
    }

    columns.resize(m3s.size());
    rvals.resize(m3s.size());

    std::vector<int> row_count(n, 0);
    for (const auto& m : m3s) {
        int index = rptrs[m.x] + row_count[m.x]++;
        columns[index] = m.y;
        rvals[index] = m.v;
    }
}


CCS::CCS(const std::vector<std::vector<double>>& matrix) {
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

CCS::CCS(const std::vector<M3>& m3s, int n) {
    cptrs.resize(n + 1, 0);

    for (const auto& m : m3s) {
        cptrs[m.y + 1]++;
    }

    for (int i = 1; i <= n; ++i) {
        cptrs[i] += cptrs[i - 1];
    }

    rows.resize(m3s.size());
    cvals.resize(m3s.size());

    std::vector<int> col_count(n, 0);
    for (const auto& m : m3s) {
        int index = cptrs[m.y] + col_count[m.y]++;
        rows[index] = m.x;
        cvals[index] = m.v;
    }
}
