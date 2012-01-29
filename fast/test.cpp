#include "rbm_mc.hpp"

#include <armadillo>
#include <boost/container/vector.hpp>
#include <iostream>

using namespace arma;
using namespace std;

int main() {
    cout << "Running ..." << endl;

    int n = 3;
    int k = 2;
    int ntrials = 10;
    vector<vec> ss = sample(n, k, ntrials);

    for (int t=0; t < ntrials; ++t) {
        cout << "Trial: " << t << endl;
        ss[t].t().print();
    }


    return 0;
}
