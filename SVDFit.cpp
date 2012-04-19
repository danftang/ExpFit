#include <Eigen/Eigenvalues>
#include "SVDFit.h"

using namespace Eigen;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
// m = number of points to interpolate
//
///////////////////////////////////////////////////////////////////////////////
SVDFit::SVDFit(int m) {
  resize(m);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SVDFit::Real SVDFit::mu_approx(int j) {
  int 	m;
  Real 	y = 0.0;

  for(m = 0; m < A.size(); ++m) {
    y += A(m)*(pow(k(m),j) - pow(k(m),mu.size()-1-j));
  }

  return(y);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SVDFit::Real SVDFit::mu_exact(double j) {
  //return(digamma((0.5+j)/mu.size())/digamma(0.5/mu.size()));
  return(tan((0.5*PI/mu.size()))/tan(PI*(j+0.5)/mu.size()));
  //if(j < mu.size()/2) { 
  // return(tan((0.5*PI/mu.size()))/tan(PI*(j+0.5)/mu.size()));
  //} else {
  //  return(-tan((0.5*PI/mu.size()))/tan(PI*(j+0.5)/mu.size()));
  //}
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SVDFit::Real SVDFit::two_norm() {
  int 	 i;
  Real 	 r;

  r = 0.0;
  for(i = 0; i < mu.size(); ++i) {
    r += pow(residual(i),2);
  }
  return(r);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SVDFit::Real SVDFit::residual(int i) {
    return(mu_approx(i) - mu(i));
}


///////////////////////////////////////////////////////////////////////////////
//
// compare error between approximant and 
// tan(PI/4N)/tan(x) for x = (0.5+i)PI/2N
//
///////////////////////////////////////////////////////////////////////////////
SVDFit::Real SVDFit::one_norm() {
  int 		i;
  SVDFit::Real 	err;
  
  err = 0.0;
  for(i = 0; i < mu.size(); ++i) {
    err += fabs(residual(i));
  }
  return(err);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void SVDFit::print() {
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
///////////////////////////////////////////////////////////////////////////////
void SVDFit::print_approximant() {
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
///////////////////////////////////////////////////////////////////////////////
void SVDFit::fit(double max_error) {
  const int			N = (mu.size()+1)/2;
  int				i,j;

  // --- Calc SVD and identify eigenvector of
  // --- first singular value below max_error
  // ----------------------------------------
  decompose();
  
  //   cout << "Eigenvalues are " << endl;
  // cout << ESolver.eigenvalues() << endl;

  j = 0;
  while(ESolver.eigenvalues().size()-1 > j && 
	fabs(ESolver.eigenvalues()(j).real()) > max_error) {
    ++j;
  }
  fit_exponents(j);

  //cout << "k's are" << endl;
  //cout << k << endl;

  fit_amplitudes();

}


///////////////////////////////////////////////////////////////////////////////
//
// Find zeroes of the polynomial y = sum_i c(i)x^i and put them in X
//
///////////////////////////////////////////////////////////////////////////////
void SVDFit::fit_exponents(int c) {
  const int 	N = 2.0*ESolver.eigenvectors().rows();
  RealMatrix	cm(N-1,N-1);
  int		i,j;

  // --- Construct the eigenpolynomial
  // --- from the eigenvalues
  // ---------------------------------
  RealVector P(N);	// eigen polynomial coefficients
  P.topRows(N/2) = ESolver.eigenvectors().col(c).real();
  P.bottomRows(N/2) = ESolver.eigenvectors().col(c).real().reverse();

  //cout << "Eigenpolynomial is " << endl;
  //cout << P << endl;

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

  /****
  cout << "---- Poly zeroes" << endl;
  for(i = 0; i<k.size(); ++i) {
    cout << k(i) << " " << eval_eigenpoly(P,k(i)) << endl;
  }
  cout << "----" << endl;
  ****/

  /******
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
  //cout << "Roots are" << endl;
  //cout << roots.eigenvalues() << endl;

  k.resize(N);
  i = 0;
  for(j = 0; j<roots.eigenvalues().size(); ++j) {
    if(fabs(roots.eigenvalues()(j).imag()) < 1e-15 &&
       roots.eigenvalues()(j).real() > 0.0 &&
       roots.eigenvalues()(j).real() < 1.0) {
      k(i) = roots.eigenvalues()(j).real();
      ++i;
    }
  }
  k.conservativeResize(i);
  ****/
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void SVDFit::fit_amplitudes() {
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

  //cout << "Amplitude singular values are" << endl;
  //cout << svd.singularValues() << endl;
  //cout << "D = " << endl;
  //cout << D << endl;
  //cout << "U^T mu = " << endl;
  //cout << svd.matrixU().transpose().topRows(M) * mu << endl;

  A.resize(M);
  A = svd.matrixV() * D * svd.matrixU().transpose().topRows(M) * mu;

  //cout << "Amplitudes are " << endl;
  //cout << A << endl;
}


///////////////////////////////////////////////////////////////////////////////
//
// reset number of data points to N
//
///////////////////////////////////////////////////////////////////////////////
void SVDFit::resize(int N) {
  int i,j;

  // --- set mu to exact values
  // --------------------------
  mu.resize(N);
  for(i = 0; i< N; ++i) {
    mu(i) = mu_exact(i);
  }

  // --- form Hankel matrix from mu
  // ------------------------------
  H.resize(N/4,N/4);
  for(i = 0; i < N/4; ++i) {
    for(j = 0; j < N/4; ++j) {
      H(i,j) = mu(i+j) + mu(i+N/2-1-j);
    }
  }
  /****
  mu.resize(N+1);
  for(i = 0; i<= N; ++i) {
    mu(i) = mu_exact(i);
  }


  H.resize(N/4,N/4);
  for(i = 0; i < N/4; ++i) {
    for(j = 0; j < N/4; ++j) {
      H(i,j) = mu(i+j) + mu(i+N/2-j);
    }
  }
  *****/
}


///////////////////////////////////////////////////////////////////////////////
//
// evaluates the c'th eigenpolynomial at point x
//
///////////////////////////////////////////////////////////////////////////////
SVDFit::Real SVDFit::eval_eigenpoly(RealVector &P, Real x) {
  int i;
  Real y;

  y = 0.0;
  for(i = P.size()-1; i>= 0; --i) {
    y = y*x + P(i);
  }
  return(y);
}

complex<double> SVDFit::eval_eigenpoly(RealVector &P, complex<double> x) {
  int i;
  complex<double> y;
  complex<double> c;

  y = 0.0;
  c.imag() = 0.0;
  for(i = P.size()-1; i>= 0; --i) {
    c.real() = P(i);
    y = y*x + c;
  }
  return(y);
}


///////////////////////////////////////////////////////////////////////////////
//
// Changes the variable of a polynomial represented in P and stores the
// result in Q so that P(x) = Q(x') where
//
// x = alpha x' + (1 - alpha)
//
///////////////////////////////////////////////////////////////////////////////
void SVDFit::change_variable(RealVector &P, const Real &alpha, RealVector &Q) {
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

  /****
  for(m = M-1; m>=0; --m) {
    lgb = m*log(beta);
    for(r = 0; r<m; ++r) {
      Q(r) += exp(lgb)*P(m);
      lgb += log(alpha/beta) + log((m-r)/(r+1.0));
    }
    Q(m) += exp(lgb)*P(m);
  }
  ****/

  /*****
  cout << "P = " << P << endl;
  cout << "Q = " << Q << endl;
  double x, px, qx;
  for(x = -1.0; x<=1.0; x+= 0.01) {
    px = eval_eigenpoly(P,alpha*x + beta);
    qx = eval_eigenpoly(Q,x);
    cout << px << " " << qx << " " << px-qx << endl; 
  }
  ******/
}


///////////////////////////////////////////////////////////////////////////////
//
// finds the zeroes of P(x) in the domain a < x <= b, putting the
// results into Z
//
///////////////////////////////////////////////////////////////////////////////
void SVDFit::find_zeroes(RealVector &P, const Real &a, const Real &b, RealVector &Z) {
  const int N = P.size();
  RealMatrix cm(N-1,N-1);
  int i,j;

  // cout << "Finding zeroes of poly size " << P.size() << endl;

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
  
  //cout << "Finding roots of " << endl;
  //cout << P << endl;
  //cout << "Complex roots are" << endl;
  //cout << roots.eigenvalues() << endl;

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


///////////////////////////////////////////////////////////////////////////////
//
// Turns H into a full Hankel (odd N) or extended-Hankel (even) matrix
// (where an extra row of the identity matrix is added)
//
///////////////////////////////////////////////////////////////////////////////
void SVDFit::createHankel() {
  int i,j;
  int N,M;

  N = mu.size()/2 + 1;
  M = mu.size() - N + 1;

  // --- form Hankel matrix from mu
  // ------------------------------
  H.resize(N,N);
  for(i = 0; i < M; ++i) {
    for(j = 0; j < N; ++j) {
	H(i,j) = mu(i+j);
    }
  }
  if(M != N) {
    for(j = 0; j < N-1; ++j) {
	H(M,j) = 0.0;
    }
    H(M,M) = 1.0;
  }
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void SVDFit::fullHfit() {
  RealVector P;

  createHankel();
  cout << "H = " << endl;
  cout << H << endl;
  
  /***
  JacobiSVD<RealMatrix> Jac(H,ComputeFullU | ComputeFullV);
  cout << "Eigenvalues are" << endl;
  cout << Jac.singularValues() << endl;
  cout << "U is" << endl;
  cout << Jac.matrixU() << endl;
  cout << "V is" << endl;
  cout << Jac.matrixV() << endl;


  P = Jac.matrixV().col(4);
  find_zeroes(P,0.0,1.0,k);
  cout << "-----" << endl;
  P = Jac.matrixU().col(4);
  find_zeroes(P,0.0,1.0,k);
  ***/

  decompose();
  cout << "Eigenvalues are" << endl;
  cout << ESolver.eigenvalues() << endl;
  cout << "Eigenvectors are" << endl;
  cout << ESolver.eigenvectors() << endl;

  int C = 0;
  while(C<ESolver.eigenvalues().size() && fabs(ESolver.eigenvalues()(C).real()) > 1e-5) ++C;
  P = ESolver.eigenvectors().col(C).real(); 

  cout << "P is " << endl << P << endl;

  find_zeroes(P,0.0,1.0,k);
  cout << "Roots are" << endl;
  cout << k << endl;

  fit_amplitudes();

  cout << "Amplitudes are" << endl;
  cout << A << endl;

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SVDFit::RealVector SVDFit::convolution(RealVector a, RealVector b) {
  RealVector result(a.size()+b.size()-1);
  int i,j;

  result.fill(0.0);
  for(i = 0; i<a.size(); ++i) {
    for(j = 0; j<b.size(); ++j) {
      result(i+j) += a(i)*b(j);
    }
  }
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void SVDFit::test() {
  int		i,j;
  int	C = 1;


  decompose();

  // --- create full Hankel with mid-row missing
  H.resize(mu.size()/2, mu.size()/2);
  for(i = 0; i < mu.size()/4; ++i) {
    for(j = 0; j < mu.size()/2; ++j) {
      H(i,j) = mu(i+j);
    }
  }
  for(i = mu.size()/4+1; i <= mu.size()/2; ++i) {
    for(j = 0; j < mu.size()/2; ++j) {
      H(i-1,j) = mu(i+j);
    }
  }
  //cout << "H = " << endl;
  //cout << H << endl;


  while(C<ESolver.eigenvalues().size() && fabs(ESolver.eigenvalues()(C).real()) > 1e-5) ++C;

  //cout << "Chosen eigenvalue is" << endl;
  //cout << ESolver.eigenvalues()(C) << endl;
  //cout << "Eigenvalues are" << endl;
  //cout << ESolver.eigenvalues() << endl;
  //cout << "Eigenvector is" << endl;
  //cout << ESolver.eigenvectors() << endl;


  RealVector E = ESolver.eigenvectors().col(C).real(); 
  RealVector E2 = -ESolver.eigenvectors().col(C+1).real(); 

  //cout << "E = " << endl << E << endl;
  //cout << "E2 = " << endl << -ESolver.eigenvectors().col(C+1).real() << endl;

  RealVector P(2*E.size());	// eigen polynomial coefficients
  RealVector AP(2*E.size());	// eigen polynomial coefficients
  RealVector APt(2*E.size());	// eigen polynomial coefficients
  RealVector AP2(4*E.size());	// eigen polynomial coefficients
  RealVector x1(2);	// eigen polynomial coefficients
  x1.fill(1.0);

  P.fill(0.0);
  AP.fill(0.0);
  AP2.fill(0.0);
  APt.fill(0.0);
  P.topRows(E.size()) = E;
  P.bottomRows(E.size()) = E.reverse();
  AP.topRows(E.size()) = E;
  AP.bottomRows(E.size()) = -E.reverse();
  AP2.topRows(2*E.size()) = 0.5*P + 0.5*AP;
  AP2.bottomRows(2*E.size()) = 0.5*P + 0.5*AP;

  //AP2.bottomRows(2*E.size()) = -AP;

  //APt.topRows(E.size()) = E.reverse();
  //APt.bottomRows(E.size()) = -E;
  //AP2.bottomRows(2*E.size()) = -APt;

  //P = ESolver.eigenvectors().col(C).real() + ESolver.eigenvectors().col(C+1).real(); 
  //AP = ESolver.eigenvectors().col(C).real() - ESolver.eigenvectors().col(C+1).real(); 

  //cout << "Sym  P = " << endl << P << endl;
  //cout << "Asym P = " << endl << AP << endl;
  //cout << "AP2 = " << endl << AP2 << endl;
  //cout << "E = " << endl << E << endl;
  //cout << "E2 = " << endl << E2 << endl;
  //cout << "Sym/Asym P = " << endl << P.cwiseProduct(AP.cwiseInverse()) << endl;
  //createHankel();
  RealVector U;

  // ---- plot characteristic poly on unit circle
  double a;
  complex<double> e;
  complex<double> e2;
  for(a = 0; a<3.1415; a += 3.14159/1024) {
    e.real() = cos(a);
    e.imag() = -sin(a);
    e2 = eval_eigenpoly(AP2,e);
    e.real() = cos(a);
    e.imag() = sin(a);
    e = eval_eigenpoly(P,e);
    cout << a << " " << abs(e) << " " << abs(e2) << " " << abs(e2/e) << endl;
  }

  U = H*P/ESolver.eigenvalues()(C).real();
  //cout << "HP = " << endl << U << endl;
  //cout << "HE = " << endl << H*E/ESolver.eigenvalues()(C).real() << endl;
  //cout << "HE2 = " << endl << H*E2/ESolver.eigenvalues()(C).real() << endl;


  //  find_zeroes(U, 0.0000001,0.999999, k);
  find_zeroes(P, 0.0000001,0.999999, k);

  //cout << "Roots are" << endl;
  //cout << k << endl;

  fit_amplitudes();

  //cout << "Amplitudes are" << endl;
  //cout << A << endl;

  /****
  RealMatrix HE(2*H.cols(), H.cols()); // Extended Hankel
  for(i = 0; i < 2*H.cols(); ++i) {
    for(j = 0; j < H.cols(); ++j) {
      HE(i,j) = mu_approx((i+j)%mu.size()) - mu((i+j) % mu.size());
    }
  }

  cout << "H_E P = " << endl << HE*P/ESolver.eigenvalues()(C).real() << endl;
  ***/
  

  return;

  //  fit_exponents(C);

  /****
  RealVector PP(H.cols());
  RealVector P(2*H.cols());

  P.fill(0.0);
  PP = ESolver.eigenvectors().col(C).real();
  P.topRows(PP.size()) = PP;
  P.bottomRows(PP.size()) += PP.reverse();

  //cout << "Half EigenPoly is" << endl;
  //cout << PP << endl;
  //cout << "EigenPoly is" << endl;
  //cout << P << endl;

  find_zeroes(P, 0.0,1.0, k);
  fit_amplitudes();

  // --- Form the companion matrix
  // -----------------------------
  int N = P.size();
  RealMatrix cm(N-1,N-1);
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
  EigenSolver<RealMatrix> roots2(cm,false);
  ****/
  //  cout << roots2.eigenvalues() << endl;
}


