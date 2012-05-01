#include <iostream>
#include <cmath>
#include <vector>
#include "AntisymmetricExpFit.cpp"

using namespace std;

int main() {
  vector<double>	data(8);
  int			N;
  int			i,j;
  const double		PI = 3.1415926535897932;
  bool 			provable = false;
  AntisymmetricExpFit	myFit(data, provable);

  for(N = 8; N < (2<<15); N *= 2) {
    // --- setup data to fit
    // ---------------------
    data.resize(N);
    for(i = 0; i<N; ++i) {
      data[i] = 1.0/(2*N*tan(PI*(i+0.5)/(2*N)));
    }

    // --- do fit
    // ----------
    myFit.set_data(data, provable);
    
    for(i = 1; i<myFit.eigenvalues().size(); ++i) {
      myFit.fit(i);

      // --- print out result
      // --------------------
      cout.precision(5);
      cout << N << " " << i << " "
      	   << myFit.one_norm() << " " 
      	   << 2.0*myFit.eigenvalues()(i).real() << endl;
      cout.precision(21);
      for(j = 0; j<i; ++j) {
      	cout << scientific << myFit.A(j) << " "; 
	cout << scientific << myFit.k(j) << endl;
      }
      cout << endl;
      if(myFit.eigenvalues()(i).real() < 1e-16) break;
    }
  }
  return(0);
}
