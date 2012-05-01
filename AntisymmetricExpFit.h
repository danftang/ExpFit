///////////////////////////////////////////////////////////////////////////////
//
//
// Copyright (c) 2012 Daniel Tang.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing,
//   software distributed under the License is distributed on an "AS
//   IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
//   express or implied.  See the License for the specific language
//   governing permissions and limitations under the License.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Core>

///////////////////////////////////////////////////////////////////////////////
//
///
/// \brief Finds a sum-of-exponentials fit to a monotonic,
/// anitsymmetric function.
///
/// The function is sampled at a set of points. One half of the points
/// are provided on construction, these are reflected about the point
/// one half step above the last point.
///
/// \author Daniel Tang
//
///////////////////////////////////////////////////////////////////////////////
class AntisymmetricExpFit {
public:
  typedef long double 						Real;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, 1>		RealVector;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> 	RealMatrix;
  typedef const Eigen::EigenSolver<RealMatrix>::EigenvalueType 	Eigenvalues;
  typedef Eigen::EigenSolver<RealMatrix>::EigenvectorsType 	Eigenvectors;

  AntisymmetricExpFit(const std::vector<double> &, bool=true);

  void   	fit(double);
  void   	fit(int);
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
  Eigenvalues &	eigenvalues();
  Eigenvectors 	eigenvectors();

  RealVector		mu;	// values of data points
  RealVector		A;	// Amplitudes of exponentials
  RealVector		k;	// exponants;
  static const Real	PI 	= 3.141592653589793238463;

protected:
  void		change_variable(RealVector &, const Real &, RealVector &);
  void		find_zeroes(RealVector &, const Real &,
			    const Real &,RealVector &);

  Eigen::EigenSolver<RealMatrix> ESolver; // eigenvalue solver
  RealMatrix			 H;	// Reduced Hankel matrix of data points
  bool				 provable; // use the provable method if true
};


///////////////////////////////////////////////////////////////////////////////
//
/// Find the eigenvectors and eigenvalues of the Hankel matrix
//
///////////////////////////////////////////////////////////////////////////////
inline void AntisymmetricExpFit::decompose() {
  ESolver.compute(H);
}

///////////////////////////////////////////////////////////////////////////////
//
/// \return A vector containing the eigenvalues of the Hankel matrix
//
///////////////////////////////////////////////////////////////////////////////
inline AntisymmetricExpFit::Eigenvalues &AntisymmetricExpFit::eigenvalues() {
  return(ESolver.eigenvalues());
}

///////////////////////////////////////////////////////////////////////////////
//
/// \return A matrix containing the eigenvectors of the Hankel matrix
//
///////////////////////////////////////////////////////////////////////////////
inline AntisymmetricExpFit::Eigenvectors AntisymmetricExpFit::eigenvectors() {
  return(ESolver.eigenvectors());
}

