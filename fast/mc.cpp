#include <iostream>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>

using namespace std;
using namespace boost::random;



int main() {
    cout << "Hello, world!" << endl;

    boost::random::mt19937 eng;

    uniform_real_distribution<> dist = uniform_real_distribution<> (0, 10);

    printf("Random number is %f", dist(eng));
    return 0;
}


