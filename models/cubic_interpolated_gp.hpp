using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Dynamic;

Matrix<var, Dynamic, 1> build_output(const MatrixXd &dvdl,
                                     const Matrix<var, Dynamic, 1> &z_,
                                     const var &l) {
  Matrix<var, Dynamic, 1> output(z_.rows());
  Matrix<double, Dynamic, 1> z(z_.rows());
  
  for(int i = 0; i < z_.rows(); i++) {
    z[i] = z_[i].val();
  }
  
  Matrix<double, Dynamic, 1> dfdl = dvdl * z;

  // For each output we need to build a special var
  for(int i = 0; i < z_.rows(); i++) {
    // Create a new var which holds a vari which holds
    //   1. The value of the output
    //   2. A pointer to the parameter on which this value depends
    //   3. The partial derivative of this output with respect to that value
    //
    // The vari is allocated with new, but we don't need to worry about
    // cleaning up the memory. It is allocated in a special place that
    // Stan handles
    output(i) = var(new precomp_v_vari(0.0, l.vi_, dfdl(i)));
  }

  return output;
}

VectorXd build_output(const MatrixXd &dvdl, const VectorXd &z, const double &l) {
  return VectorXd::Zero(z.size());  
}

template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__>::type>::type, Eigen::Dynamic,1>
approx_Lz(const T0__& l,
          const std::vector<T1__>& lp,
          const std::vector<Eigen::Matrix<T2__, Eigen::Dynamic,Eigen::Dynamic> >& Ls,
          const std::vector<Eigen::Matrix<T3__, Eigen::Dynamic,Eigen::Dynamic> >& dLdls,
          const Eigen::Matrix<T4__, Eigen::Dynamic,1>& z, std::ostream* pstream__) {
  using namespace Eigen;
  int N = Ls.front().rows();
  
  int lidx = 0;
  for(; lidx < lp.size() - 1; lidx++) {
    if(lp[lidx + 1] >= l)
      break;
  }
  
  double x1 = lp[lidx],
    x2 = lp[lidx + 1];
  
  double t = (value_of(l) - x1) / (x2 - x1);
  double dtdl = 1 / (x2 - x1);
  
  const Eigen::MatrixXd &y1 = Ls[lidx];
  const Eigen::MatrixXd &y2 = Ls[lidx + 1];
  
  const Eigen::MatrixXd &k1 = dLdls[lidx];
  const Eigen::MatrixXd &k2 = dLdls[lidx + 1];
  
  Eigen::MatrixXd a = k1 * (x2 - x1) - (y2 - y1);
  Eigen::MatrixXd b = -k2 * (x2 - x1) + (y2 - y1);
  
  Eigen::MatrixXd v = (1 - t) * y1 + t * y2 + t * (1 - t) * (a * (1 - t) + b * t);
  Eigen::MatrixXd dvdl = (b * (2 - 3 * t) * t + a * (1 + t * (-4 + 3 * t)) - y1 + y2) * dtdl;
  
  return multiply(v, z) + build_output(dvdl, z, l);
}