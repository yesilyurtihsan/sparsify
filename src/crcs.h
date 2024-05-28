#ifndef CRCS_H
#define CRCS_H

#include <vector>

class M3 {
    public:
    int x;
    int y;
    int v;

    M3(int x, int y, int v): x(x), y(y), v(v) {}
    int getX() const { return x; }
    int getY() const { return y; }
    int getV() const { return v; }
};

class CRS {
public:
    std::vector<int> rptrs;
    std::vector<int> columns;
    std::vector<double> rvals;

    CRS(const std::vector<std::vector<double>>& matrix);
    CRS(const std::vector<M3>& m3s, int n);
};

class CCS {
public:
    std::vector<int> cptrs;
    std::vector<int> rows;
    std::vector<double> cvals;

    CCS(const std::vector<std::vector<double>>& matrix);
    CCS(const std::vector<M3>& m3s, int n);
};

#endif
