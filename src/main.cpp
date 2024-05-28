#include <iostream>
#include <fstream>
#include <vector>
#include "crcs.h"
#include "spa_ryser.cpp"

using namespace std;


int main(int argc, char ** argv){
    int n, nonzeros;
    vector<M3> m3s;
    ifstream fin(argv[1]);
    fin >> n >> nonzeros;

    for(int i=0; i<nonzeros; i++){
        int x, y, v;
        fin >> x >> y >> v;
        M3 m(x, y, v);
        m3s.push_back(m);
    }

    // int n = 3, nonzeros = 9;
    // vector<M3> m3s;
    // m3s.push_back(M3(0, 0, 1));
    // m3s.push_back(M3(0, 1, 2));
    // m3s.push_back(M3(0, 2, 3));
    // m3s.push_back(M3(1, 0, 4));
    // m3s.push_back(M3(1, 1, 5));
    // m3s.push_back(M3(1, 2, 6));
    // m3s.push_back(M3(2, 0, 7));
    // m3s.push_back(M3(2, 1, 8));
    // m3s.push_back(M3(2, 2, 9));

    CRS crs(m3s, n);
    CCS ccs(m3s, n);

    double perm = SpaRyser(crs, ccs);
    cout << perm << endl;

    return 0;
}