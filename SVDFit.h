///////////////////////////////////////////////////////////////////////////////
//
// Finds the best sum-of-exponentials fit to a set of points
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/SVD>
#include "digamma.cpp"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class SVDFit {
public:
  typedef long double 						Real;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, 1>		RealVector;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> 	RealMatrix;

  SVDFit(int);

  Real	 	mu_approx(int);
  Real	 	mu_exact(double);
  Real	 	residual(int);
  Real	 	two_norm();
  Real		one_norm();
  void   	fit(double);
  void		fullHfit();
  void 		print();
  void 		print_approximant();
  void		resize(int);
  void		decompose();
  void		fit_exponents(int);
  void		fit_amplitudes();
  Real		eval_eigenpoly(RealVector &, Real);
  std::complex<double>	eval_eigenpoly(RealVector &, std::complex<double>);

  Eigen::JacobiSVD<RealMatrix>	SVD;	// decomposition solver
  Eigen::EigenSolver<RealMatrix> ESolver;	// eigenvalue solver
  RealVector		mu;	// values of data points
  RealMatrix		H;	// Reduced Hankel matrix of data points
  RealVector		A;	// Amplitudes of exponentials
  RealVector		k;	// exponants;
  static const Real	PI 	= 3.141592653589793238463;

protected:
  void		change_variable(RealVector &, const Real &, RealVector &);
  void		find_zeroes(RealVector &, const Real &, const Real &, RealVector &);

public:
  // testy stuff

  void		test();
  void		createHankel();
  RealVector	convolution(RealVector , RealVector );
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline void SVDFit::decompose() {
  //  SVD.compute(H,Eigen::ComputeFullU);
  ESolver.compute(H);
}
