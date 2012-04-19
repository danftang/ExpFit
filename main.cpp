#include <iostream>
#include <valarray>
#include "SVDFit.cpp"

int main() {
  SVDFit		myFit(32);
  int			n;
  int			i,j,k;
  const int		N=13;	// log max number of points
  const int		M=24;	// max number of exponents
  double		e;
  double		f;
  double		extra_comp;
  double		err[M][N];	// log10 error given terms and points
  double		ecomp[M][N];	// extra computation needed for decay

  for(i = 0; i<M; ++i) {
    for(n = 0; n<N; ++n) {
      err[i][n] = -12.0;
      ecomp[i][n] = 0.0;
    }
  }

  for(n = 3; n<N; ++n) {
    myFit.resize(1<<n);
    myFit.decompose();
    e = 0.0;
    for(i=1; 
	i<myFit.ESolver.eigenvalues().size() && 
	  myFit.ESolver.eigenvalues()(i).real() > 1e-15 &&
	  e > -11.0
	  ; ++i) {
     myFit.fit_exponents(i);
     myFit.fit_amplitudes();

     e = log10(myFit.one_norm());
     err[myFit.k.size()][n] = e;

     extra_comp = 0.0;
     for(j = 0; j<myFit.k.size(); ++j) {
       f = (e - log10(myFit.A(j)))/log10(myFit.k(j));
       if(f > (1<<n)) f = 1<<n;
       if(f < 0.5) f = 0.5;
       extra_comp += (2.0*f - 1)/(1<<n);
     }
     ecomp[myFit.k.size()][n] = extra_comp;
     std::cout << n << " " << e << " " 
	       << 7.0*myFit.k.size() + 2.0*extra_comp - 1.0 << " " 
	       << 7.0*myFit.k.size() - 1.0 << " " 
	       << 2.0*extra_comp << std::endl;
    }
  }

  /****
  for(i = 1; i<M; ++i) {
    for(n = 3; n<N; ++n) {
      std::cout << i << " " << n << " " << err[i][n] << " " << ecomp[i][n] << std::endl;
    }
    std::cout << std::endl;
  }  
  *****/

  //myFit.fit(1e-5);
  //myFit.test();
  //myFit.fullHfit();

  return(0);

  //std::cout << "----" << std::endl;
  //std::cout << myFit.k << std::endl;
  //std::cout << "----" << std::endl;
  //std::cout << myFit.A << std::endl;
  std::cout << "# y = ";
  myFit.print_approximant(); 
  std::cout << "# one norm = " << myFit.one_norm() << std::endl;
  std::cout << std::endl;
  myFit.print();

  return(0);
}
