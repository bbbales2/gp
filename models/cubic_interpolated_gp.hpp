template <typename T0__, typename T1__, typename T2__, typename T3__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type,
              Eigen::Dynamic,Eigen::Dynamic>
approx_L(const T0__& l,
         const std::vector<T1__>& lp,
         const std::vector<Eigen::Matrix<T2__, Eigen::Dynamic,Eigen::Dynamic> >& Ls,
         const std::vector<Eigen::Matrix<T3__, Eigen::Dynamic,Eigen::Dynamic> >& dLdls, std::ostream* pstream__) {
  typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type var;
  int N = Ls.front().rows();
  Eigen::Matrix<var, Eigen::Dynamic,Eigen::Dynamic> output(N, N);
  
  int lidx = 0;
  for(; lidx < lp.size() - 1; lidx++) {
    if(lp[lidx + 1] >= l)
      break;
  }
  
  double x1 = lp[lidx],
    x2 = lp[lidx + 1];
  
  var t = (l - x1) / (x2 - x1);
  
  //std::cout << lidx << std::endl;
  //std::cout << stan::math::value_of(alpha) << std::endl;
  
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      if(j <= i) {
        double y1 = Ls[lidx](i, j),
          y2 = Ls[lidx + 1](i, j);
          
        double k1 = dLdls[lidx](i, j),
          k2 = dLdls[lidx + 1](i, j);
            
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