#include "rbm.hpp"

#include <armadillo>
#include <boost/container/vector.hpp>
#include <iostream>

using namespace arma;
using namespace std;

int main() {
    cout << "Running ..." << endl;

    int n = 3;
    int kmax = 4;

    vector< vector<vec> > ss;

    for (int k=0; k <= kmax; ++k) {
        ss.push_back(sample_lattice(n, k));
    }

    for (int k=0; k <= kmax - 1; ++k) {
        double d = hausdorff(ss[k], ss[k+1]);
        cout << "Hausdorff distance between k=" << k 
             << " and k=" << (k+1) 
             << " is: " << d
             << endl;
    }
    return 0;
}
