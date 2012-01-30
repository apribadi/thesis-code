#include "rbm.hpp"

#include <armadillo>
#include <boost/container/vector.hpp>
#include <iostream>

using namespace arma;
using namespace std;

int main() {
    cout << "Running ..." << endl;

    int n = 2;
    int kmax = 3;
    int ntrials = 1 << 14;

    vector< vector<vec> > ss;
    vector< vector<vec> > tt;

    for (int k=0; k <= kmax; ++k) {
        ss.push_back(sample(n, k, ntrials));
        tt.push_back(sample(n, k, ntrials));
    }

    for (int k=0; k <= kmax - 1; ++k) {
        double d = hausdorff(ss[k], ss[k+1]);
        cout << "Hausdorff distance between k=" << k 
             << " and k=" << (k+1) 
             << " is: " << d
             << endl;
    }

    for (int k=0; k <= kmax; ++k) {
        double d = hausdorff(ss[k], tt[k]);
        cout << "Hausdorff distance between k=" << k 
             << " and k=" << (k) 
             << " is: " << d
             << endl;
    }
    return 0;
}
