///////////////////////////////////////////////////////////////////////////////
//
// Finds the best sum-of-exponentials fit to a set of points
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Core>

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class AntisymmetricExpFit {
public:
  typedef long double 						Real;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, 1>		RealVector;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> 	RealMatrix;

  AntisymmetricExpFit(const std::vector<double> &);

  void   	fit(double);
  void		set_data(const std::vector<double> &);

  Real	 	mu_approx(int);
  Real	 	residual(int);
  Real	 	two_norm();
  Real		one_norm();
  void 		print();
  void 		print_approximant();
  
  void		decompose();
  void		fit_exponents(int);
  void		fit_amplitudes();

  Eigen::EigenSolver<RealMatrix> ESolver;	// eigenvalue solver
  RealVector		mu;	// values of data points
  RealMatrix		H;	// Reduced Hankel matrix of data points
  RealVector		A;	// Amplitudes of exponentials
  RealVector		k;	// exponants;
  static const Real	PI 	= 3.141592653589793238463;

protected:
  void		change_variable(RealVector &, const Real &, RealVector &);
  void		find_zeroes(RealVector &, const Real &, const Real &, RealVector &);

};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline void AntisymmetricExpFit::decompose() {
  ESolver.compute(H);
}
