using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Dynamic;

Matrix<var, Dynamic, Dynamic> build_output(const MatrixXd &v,
                                           const MatrixXd &dvdl,
                                           const var &l) {
  Matrix<var, Dynamic, Dynamic> output(v.rows(), v.cols());

  // For each output we need to build a special var
  for(int i = 0; i < v.rows(); i++) {
    for(int j = 0; j < v.cols(); j++) {
      // Create a new var which holds a vari which holds
      //   1. The value of the output
      //   2. A pointer to the parameter on which this value depends
      //   3. The partial derivative of this output with respect to that value
      //
      // The vari is allocated with new, but we don't need to worry about
      // cleaning up the memory. It is allocated in a special place that
      // Stan handles
      output(i, j) = var(new precomp_v_vari(v(i, j), l.vi_, dvdl(i, j)));
    }
  }
   
  return output;
}

MatrixXd build_output(const MatrixXd &v, const MatrixXd &dvdl, const double &l) {
  return v;  
}

template <typename T0__, typename T1__, typename T2__, typename T3__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type,
              Eigen::Dynamic,Eigen::Dynamic>
approx_L(const T0__& l,
         const std::vector<T1__>& lp,
         const std::vector<Eigen::Matrix<T2__, Eigen::Dynamic,Eigen::Dynamic> >& Ls,
         const std::vector<Eigen::Matrix<T3__, Eigen::Dynamic,Eigen::Dynamic> >& dLdls, std::ostream* pstream__) {
  using namespace Eigen;
  typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type var;
  int N = Ls.front().rows();
  Eigen::Matrix<var, Dynamic, Dynamic> output(N, N);
  
  int lidx = 0;
  for(; lidx < lp.size() - 1; lidx++) {
    if(lp[lidx + 1] >= l)
      break;
  }
  
  double x1 = lp[lidx],
    x2 = lp[lidx + 1];
  
  double t = (value_of(l) - x1) / (x2 - x1);
  double dtdl = 1 / (x2 - x1);
  
  //std::cout << lidx << std::endl;
  //std::cout << stan::math::value_of(alpha) << std::endl;
  
  const Eigen::MatrixXd &y1 = Ls[lidx];
  const Eigen::MatrixXd &y2 = Ls[lidx + 1];
  
  const Eigen::MatrixXd &k1 = dLdls[lidx];
  const Eigen::MatrixXd &k2 = dLdls[lidx + 1];
  
  Eigen::MatrixXd a = k1 * (x2 - x1) - (y2 - y1);
  Eigen::MatrixXd b = -k2 * (x2 - x1) + (y2 - y1);
  
  Eigen::MatrixXd v = (1 - t) * y1 + t * y2 + t * (1 - t) * (a * (1 - t) + b * t);
  Eigen::MatrixXd dvdl = (b * (2 - 3 * t) * t + a * (1 + t * (-4 + 3 * t)) - y1 + y2) * dtdl;
  
  return build_output(v, dvdl, l);
}