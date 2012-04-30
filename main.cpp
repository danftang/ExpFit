#include <iostream>
#include <cmath>
#include <vector>
#include "AntisymmetricExpFit.cpp"

using namespace std;

int main() {
  const int		N = 16;	// half number of points to fit
  vector<double>	data(N);
  int			i;
  const double		PI = 3.1415926535897932;

  // --- setup data to fit
  // ---------------------
  for(i = 0; i<N; ++i) {
    data[i] = 1.0/(2*N*tan(PI*(i+0.5)/(2*N)));
  }

  // --- do fit
  // ----------
  bool 			provable = true; // set to false for alternative method 
  AntisymmetricExpFit	myFit(data, provable);
  myFit.fit(1e-5);

  // --- print out result
  // --------------------
  cout << "# y = ";
  myFit.print_approximant(); 
  cout << "# one norm = " << myFit.one_norm() << endl;
  cout << "# rms = " << sqrt(myFit.two_norm()/myFit.mu.size()) << endl;
  cout << endl;
  myFit.print();

  return(0);
}
