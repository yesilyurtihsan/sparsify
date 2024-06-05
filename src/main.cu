#include <vector>
#include <cmath>
#include <iostream>
#include <chrono>
#include <cuda_runtime.h>
#include <fstream>
#include "crcs.h"
#include "cuda/spaRyser.cu"

int main(int argc, char ** argv) {
    cudaSetDevice(1);
    // int n, nonzeros;
    // std::vector<M3> m3s;
    // std::ifstream fin(argv[1]);
    // fin >> n >> nonzeros;

    // for(int i=0; i<nonzeros; i++){
    //     int x, y;
    //     double v;
    //     fin >> x >> y >> v;
    //     M3 m(x, y, v);
    //     m3s.push_back(m);
    // }

    // CRS crs(m3s, n);
    // CCS ccs(m3s, n);

    int n = 3, nonzeros = 9;
    std::vector<M3> m3s;
    m3s.push_back(M3(0, 0, 1));
    m3s.push_back(M3(0, 1, 2));
    m3s.push_back(M3(0, 2, 3));
    m3s.push_back(M3(1, 0, 4));
    m3s.push_back(M3(1, 1, 5));
    m3s.push_back(M3(1, 2, 6));
    m3s.push_back(M3(2, 0, 7));
    m3s.push_back(M3(2, 1, 8));
    m3s.push_back(M3(2, 2, 9)); // RESULT: 159808

    CRS crs(m3s, n);
    CCS ccs(m3s, n);

    std::cout << "CRS Representation:" << std::endl;
    std::cout << "rptrs: ";
    for (auto val : crs.rptrs) std::cout << val << " ";
    std::cout << std::endl << "columns: ";
    for (auto val : crs.columns) std::cout << val << " ";
    std::cout << std::endl << "rvals: ";
    for (auto val : crs.rvals) std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "CCS Representation:" << std::endl;
    std::cout << "cptrs: ";
    for (auto val : ccs.cptrs) std::cout << val << " ";
    std::cout << std::endl << "rows: ";
    for (auto val : ccs.rows) std::cout << val << " ";
    std::cout << std::endl << "cvals: ";
    for (auto val : ccs.cvals) std::cout << val << " ";
    std::cout << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    auto permanent = SpaRyser<float>(crs, ccs);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Permanent of the matrix is: " << permanent << std::endl;
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << "ns" << std::endl;

    return 0;
}
