///////////////////////////////////////////////////////////////////////////////
//
// Finds the best sum-of-exponentials fit to a set of points
//
// Daniel F. Tang
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
  typedef const Eigen::EigenSolver<RealMatrix>::EigenvalueType  Eigenvalues;
  typedef const Eigen::EigenSolver<RealMatrix>::EigenvectorsType Eigenvectors;

  AntisymmetricExpFit(const std::vector<double> &, bool=true);

  void   	fit(double);
  void		set_data(const std::vector<double> &, bool=true);

  Real	 	mu_approx(int);
  Real	 	residual(int);
  Real	 	two_norm();
  Real		one_norm();
  void 		print();
  void 		print_approximant();

  void		decompose();
  void		fit_exponents(int);
  void		fit_amplitudes();
  Eigenvalues & eigenvalues();
  Eigenvectors &eigenvectors();

  RealVector		mu;	// values of data points
  RealVector		A;	// Amplitudes of exponentials
  RealVector		k;	// exponants;
  static const Real	PI 	= 3.141592653589793238463;

protected:
  void		change_variable(RealVector &, const Real &, RealVector &);
  void		find_zeroes(RealVector &, const Real &,
			    const Real &,RealVector &);

  RealMatrix			 H;	// Reduced Hankel matrix of data points
  Eigen::EigenSolver<RealMatrix> ESolver; // eigenvalue solver
  bool				 provable; // use the provable method if true
};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline void AntisymmetricExpFit::decompose() {
  ESolver.compute(H);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline AntisymmetricExpFit::Eigenvalues &AntisymmetricExpFit::eigenvalues() {
  ESolver.eigenvalues();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline AntisymmetricExpFit::Eigenvectors &AntisymmetricExpFit::eigenvectors() {
  ESolver.eigenvectors();
}

