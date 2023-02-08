#include <Rcpp.h>
using namespace Rcpp;
// this cpp script is not necessary because I directly define the function with `cppFunction` in `allele_freq_weighted_similarity.R`

// [[Rcpp::export]]
NumericMatrix rcpp_freq_w_similarity(NumericMatrix mat, NumericVector freq1) {

  // allocate the matrix to return
  NumericMatrix rmat(mat.nrow(), mat.nrow());
  for (int i = 0; i < rmat.nrow(); i++) {
    for (int j = 0; j < i; j++) {
      // computation ver1
      //NumericMatrix::Row x1 = mat.row(i);
      //NumericMatrix::Row y1 = mat.row(j);
      //NumericVector x0 = 2 - x1;
      //NumericVector y0 = 2 - y1;
      //NumericVector p1 = (x1 * y1)*0.25;
      //NumericVector p0 = (x0 * y0)*0.25;
      //NumericVector freq0 = 1 - freq1;
      //double s = mean(freq0*p1 + freq1*p0); // freq0 = 1-freq1 and freq1 = 1-freq0

      // computation ver2
      NumericMatrix::Row x1 = mat.row(i);
      NumericMatrix::Row y1 = mat.row(j);
      std::vector<double> v(x1.size());
      for (int k = 0; k < x1.size() ;k++) {
        v[k] = (1 - freq1[k])*((x1[k] * y1[k])*0.25) + freq1[k]*(((2 - x1[k]) * (2 - y1[k]))*0.25);
      }
      double s = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
      // write to output matrix
      rmat(i,j) = s;
    }
  }
  return rmat;
}
