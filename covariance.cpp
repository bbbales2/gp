#include <Rcpp.h>
#include <cmath>
#include "stan/math/fwd/mat.hpp"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
List rbf_cov_chol(NumericVector x1, double l_) {
  using namespace Eigen;
  using namespace stan;
  
  math::fvar<double> l(l_, 1.0);

  Eigen::Matrix<math::fvar<double>, Dynamic, Dynamic> Sigma(x1.size(), x1.size());
  
  for(int i = 0; i < x1.size(); i++) {
    for(int j = 0; j < x1.size(); j++) {
      Sigma(i, j) = math::exp(-(x1[i] - x1[j]) * (x1[i] - x1[j]) / (2 * l * l));
    }
  }
  
  for(int i = 0; i < x1.size(); i++) {
    Sigma(i, i) += 1e-10;
  }
  
  Eigen::Matrix<math::fvar<double>, Dynamic, Dynamic> L_(x1.size(), x1.size());
  
  L_ = math::cholesky_decompose(Sigma);
  
  NumericMatrix L(x1.size(), x1.size());
  NumericMatrix dLdl(x1.size(), x1.size());
  
  for(int i = 0; i < x1.size(); i++) {
    for(int j = 0; j < x1.size(); j++) {
      L(i, j) = L_(i, j).val();
      dLdl(i, j) = L_(i, j).tangent();
    }
  }
  
  List output;
  
  output["L"] = L;
  output["dLdl"] = dLdl;
  
  return output;
}

NumericMatrix approx_L(double l,
         NumericVector lp,
         List Ls,
         List dLdls) {
  const NumericMatrix &m = Ls[1];
  int N = m.nrow();
  NumericMatrix output(N, N);
  
  int lidx = 0;
  for(; lidx < lp.size() - 1; lidx++) {
    if(lp[lidx + 1] >= l)
      break;
  }
  
  double x1 = lp[lidx];
  double x2 = lp[lidx + 1];
  
  double t = (l - x1) / (x2 - x1);
  
  //std::cout << lidx << std::endl;
  //std::cout << stan::math::value_of(alpha) << std::endl;
  
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      if(j <= i) {
        const NumericMatrix &Ls1 = Ls[lidx];
        const NumericMatrix &Ls2 = Ls[lidx + 1];
        const NumericMatrix &dLdls1 = dLdls[lidx];
        const NumericMatrix &dLdls2 = dLdls[lidx + 1];
        
        double y1 = Ls1(i, j);
        double y2 = Ls2(i, j);
          
        double k1 = dLdls1(i, j);
        double k2 = dLdls2(i, j);
            
        double a = k1 * (x2 - x1) - (y2 - y1);
        double b = -k2 * (x2 - x1) + (y2 - y1);
            
        output(i, j) = (1 - t) * y1 + t * y2 + t * (1 - t) * (a * (1 - t) + b * t);
      } else {
        output(i, j) = 0.0;
      }
    }
  }
  
  return output;
}
