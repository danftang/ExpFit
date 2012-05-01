///////////////////////////////////////////////////////////////////////////////
//
/// \author Daniel Tang
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
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include "AntisymmetricExpFit.h"

using namespace Eigen;
using namespace std;


///////////////////////////////////////////////////////////////////////////////
//
/// \fn AntisymmetricExpFit::AntisymmetricExpFit(const std::vector<double> &y, bool p)
/// Constructs an instance of AsymmetricExpFit, forms the Hankel matrix
/// and solves for the eigenvectors.
///
/// \param y Half of the data we wish to find an antisymmetric fit of.
/// The point of antisymmetry lies at y.size()-0.5. y.size() should be
/// even.
///
/// \param p if true, the provable method of solving is used,
/// otherwise the alternative method is used. See accompanying paper for
/// further details.
//
///////////////////////////////////////////////////////////////////////////////
AntisymmetricExpFit::AntisymmetricExpFit(const vector<double> &y, bool p) {
  set_data(y,p);
}


///////////////////////////////////////////////////////////////////////////////
//
/// Fits an exponential sum to the data in the Hankel matrix, H,
/// whose error is bounded by max_error
//
///////////////////////////////////////////////////////////////////////////////
void AntisymmetricExpFit::fit(double max_error) {
  int j;

  // --- Identify eigenvector of
  // --- first eigenvalue below max_error
  // ------------------------------------
  j = 0;
  while(ESolver.eigenvalues().size()-1 > j && 
	fabs(ESolver.eigenvalues()(j).real()) > max_error/2.0) {
    ++j;
  }
  fit(j);
}


///////////////////////////////////////////////////////////////////////////////
//
/// Fits an exponential sum with 'n' terms to the data in the Hankel
/// matrix, H.
//
///////////////////////////////////////////////////////////////////////////////
void AntisymmetricExpFit::fit(int n) {
  fit_exponents(n);
  if(k.size() > 0) fit_amplitudes();
}


///////////////////////////////////////////////////////////////////////////////
//
/// \fn void AntisymmetricExpFit::set_data(const std::vector<double> &y, bool p)
/// Used to define the data points we wish to fit the approximant to
/// (if different from that sent to the constructor).
///
/// \param y Half of the data to fit. The point of antisymmetry lies
/// at y.size()-0.5. y.size() should be even
///
/// \param p if true, the provable method of solving is used,
/// otherwise the alternative method is used. See accompanying paper for
/// further details.
//
///////////////////////////////////////////////////////////////////////////////
void AntisymmetricExpFit::set_data(const vector<double> &y, bool p) {
  int i,j;
  int N = 2*y.size();
  provable = p;

  // --- set mu to exact values
  // --------------------------
  mu.resize(N);
  for(i = 0; i< N/2; ++i) {
    mu(i) = y[i];
    mu(N-1-i) = -y[i];
  }

  // --- musum = x*mu/(1-x)
  // ----------------------
  RealVector musum(N+1);
  musum(0) = 0.0;
  for(i = 1; i<musum.size(); ++i) {
    musum(i) = musum(i-1) + mu(i-1);
  }

  // --- form Hankel matrix from mu
  // ------------------------------
  H.resize(N/4,N/4);
  for(i = 0; i < N/4; ++i) {
    for(j = 0; j < N/4; ++j) {
      if(provable) {
	H(i,j) = musum(i+j) - musum(i+N/2-j);
      } else {
	H(i,j) = mu(i+j) + mu(i+N/2-1-j);
      }
    }
  }
  decompose();
}


///////////////////////////////////////////////////////////////////////////////
//
/// returns the value of the approximant at the j'th point.
//
///////////////////////////////////////////////////////////////////////////////
AntisymmetricExpFit::Real AntisymmetricExpFit::mu_approx(int j) {
  int 	m;
  Real 	y = 0.0;

  for(m = 0; m < A.size(); ++m) {
    y += A(m)*(pow(k(m),j) - pow(k(m),mu.size()-1-j));
  }

  return(y);
}


///////////////////////////////////////////////////////////////////////////////
//
/// returns the 2-norm error of the approximant
//
///////////////////////////////////////////////////////////////////////////////
AntisymmetricExpFit::Real AntisymmetricExpFit::two_norm() {
  int 	 i;
  Real 	 r;

  r = 0.0;
  for(i = 0; i < mu.size(); ++i) {
    r += pow(residual(i),2);
  }
  return(r);
}


///////////////////////////////////////////////////////////////////////////////
//
/// returns the approximant minus the exact value at point \b i
//
///////////////////////////////////////////////////////////////////////////////
AntisymmetricExpFit::Real AntisymmetricExpFit::residual(int i) {
    return(mu_approx(i) - mu(i));
}


///////////////////////////////////////////////////////////////////////////////
//
/// returns the 1-norm error of the approximant
//
///////////////////////////////////////////////////////////////////////////////
AntisymmetricExpFit::Real AntisymmetricExpFit::one_norm() {
  int 		i;
  AntisymmetricExpFit::Real 	err;
  
  err = 0.0;
  for(i = 0; i < mu.size(); ++i) {
    err += fabs(residual(i));
  }
  return(err);
}


///////////////////////////////////////////////////////////////////////////////
//
/// Prints out the values of the approximant against the data
//
///////////////////////////////////////////////////////////////////////////////
void AntisymmetricExpFit::print() {
  int i;

  cout << scientific;
  for(i = 0; i < mu.size(); ++i) {
    cout << i << " "
	      << mu(i) << " "
	      << mu_approx(i) << " "
	      << residual(i) << endl;
  }
}


///////////////////////////////////////////////////////////////////////////////
//
/// Prints the approximant to standard out
//
///////////////////////////////////////////////////////////////////////////////
void AntisymmetricExpFit::print_approximant() {
  int m;
  double x;

  cout << scientific;
  for(m = 0; m < A.size()-1; ++m) {
    cout << A(m) << "." << k(m) << "^x + ";
  }
  cout << A(m) << "." << k(m) << "^x";
  cout << endl;
}



///////////////////////////////////////////////////////////////////////////////
//
/// Find zeroes of the polynomial \f$y = sum_i e_ix^i\f$,
/// where \f$e\f$ is the \f$c^{th}\f$ eigenvector, and put any zeroes
/// that lie in the interval [0:1] into \c k.
///
/// \param c The identifier of the eigenvector to find zeroes of.
//
///////////////////////////////////////////////////////////////////////////////
void AntisymmetricExpFit::fit_exponents(int c) {
  const int 	N = 2.0*ESolver.eigenvectors().rows();
  int		i,j;

  // --- Construct the eigenpolynomial
  // --- from the eigenvalues
  // ---------------------------------
  RealVector P(N);	// eigen polynomial coefficients

  if(provable) {
    // --- PA = P_q(x) - x^{N/2} P_q(x^{-1})
    // --- P = PA/(1-x)
    RealVector PA(N+1);	// eigen polynomial coefficients
    PA(N/2) = 0.0;
    PA.topRows(N/2) = ESolver.eigenvectors().col(c).real();
    PA.bottomRows(N/2) = -ESolver.eigenvectors().col(c).real().reverse();
    P(0) = PA(0);
    for(i = 1; i<P.size(); ++i) {
      P(i) = P(i-1) + PA(i);
    }
  } else {
    P.topRows(N/2) = ESolver.eigenvectors().col(c).real();
    P.bottomRows(N/2) = ESolver.eigenvectors().col(c).real().reverse();
  }

  RealVector Q;		// eigenpolynomial with changed variable
  RealVector Z;		// scaled positions of zeroes
  long double alpha;    // scaling factor of polynomial
  long double s;	// scale
  k.resize(c);
  i = 0;
  alpha = 2.0/3.0;
  if(P.size() > 64) {
    while(i < c && alpha > 1e-8) {
      change_variable(P,alpha,Q);
      j = Q.size()-1;
      s = pow(2.0,j);
      while(j && fabs(Q(j))/s < 1e-18) {
	s /= 2.0;
	--j;
      }
      Q.conservativeResize(j+1);
      find_zeroes(Q, -0.5, 0.5, Z);
      for(j = Z.size()-1; j>=0; --j) {
	s = 1.0 - alpha*(1.0 -Z(j));
	k(i) = s;
	++i;
      }
      alpha /= 3.0;
    }
  } else {
    find_zeroes(P, 0.0, 1.0, k);
  }
}


///////////////////////////////////////////////////////////////////////////////
//
/// Solevs the Vandermonde system to find the least squares fit of the 
/// amplitudes for a given \c k.
//
///////////////////////////////////////////////////////////////////////////////
void AntisymmetricExpFit::fit_amplitudes() {
  const int 	M = k.size();
  RealMatrix	T(mu.size(),M);
  RealMatrix	D(M,M);
  int		i,j;

  for(i = 0; i<mu.size(); ++i) {
    for(j = 0; j<M; ++j) {
      T(i,j) = pow(k(j),i) - pow(k(j),mu.size()-1-i);
    }
  }

  JacobiSVD<RealMatrix> svd(T, ComputeFullU | ComputeFullV);
  D.fill(0.0);
  D.diagonal() = svd.singularValues();
  for(i = 0; i<M; ++i) {
    if(D(i,i) > 1e-12) {
      D(i,i) = 1.0/D(i,i);
    } else {
      D(i,i) = 0.0;
    }
  }  

  A.resize(M);
  A = svd.matrixV() * D * svd.matrixU().transpose().topRows(M) * mu;
}


///////////////////////////////////////////////////////////////////////////////
//
/// Changes the variable of a polynomial represented in P and stores the
/// result in Q so that \f$P(x) = Q(x')\f$ where
///
/// \f$x = \alpha x' + (1 - \alpha)\f$
//
///////////////////////////////////////////////////////////////////////////////
void AntisymmetricExpFit::change_variable(RealVector &P, const Real &alpha, RealVector &Q) {
  const int M = P.size(); 	// degree of P + 1
  const int SMAX = 96;		// maximum degree of Q
  int m;			// index into P
  int r;			// index into Q
  Real beta = 1.0-alpha;	// offset
  Real lgb;			// log of alpha^r beta^(m-r) (m choose r) 
  Q.resize(SMAX);
  Q.fill(0.0);

  for(r = 0; r<SMAX; ++r) {
    lgb = r*log(alpha);
    for(m = r; m<M; ++m) {
      Q(r) += exp(lgb)*P(m);
      lgb += log(beta) + log((m+1.0)/(m+1.0-r));
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
//
/// Finds the zeroes of P(x) in the domain \f$a < x <= b\f$, putting the
/// results into Z
//
///////////////////////////////////////////////////////////////////////////////
void AntisymmetricExpFit::find_zeroes(RealVector &P, const Real &a, const Real &b, RealVector &Z) {
  const int N = P.size();
  RealMatrix cm(N-1,N-1);
  int i,j;

  // --- Form the companion matrix
  // -----------------------------
  cm.fill(0.0);
  for(j=0; j<N-1; ++j) {
    cm(0,j) = -P(N-2-j)/P(N-1);
  }
  for(i=1; i<N-1; ++i) {	
    cm(i,i-1) = 1.0;
  }

  // --- find eigenvalues of
  // --- companion matrix and
  // --- extract roots in [0:1]
  // --------------------------
  EigenSolver<RealMatrix> roots(cm,false);

  Z.resize(N);
  i = 0;
  for(j = 0; j<roots.eigenvalues().size(); ++j) {
    if(fabs(roots.eigenvalues()(j).imag()) < 1e-15 &&
       roots.eigenvalues()(j).real() > a &&
       roots.eigenvalues()(j).real() <= b) {
      Z(i) = roots.eigenvalues()(j).real();
      ++i;
    }
  }
  Z.conservativeResize(i);
}
