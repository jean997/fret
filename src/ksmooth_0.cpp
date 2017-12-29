#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//' @title Kernel smoother for integer locations
//' @description Assumes missing locations have observation of 0
//' @param x Input locations (integers), sorted
//' @param y Input observations, should be length of x and non-missing
//' @param xout Output locations
//' @param bandwidth Bandwidth. Should be an odd integer
// [[Rcpp::export]]
NumericVector ksmooth_0_cpp(NumericVector x,  NumericVector y,
                            NumericVector xout, double bandwidth){
  int k = xout.size();

  NumericVector yout(k);
  double s;
  double left;
  double right;
  int ixleft;
  int ixright;
  for(size_t i=0; i < k; i++){
    left =xout[i] - bandwidth/2.0;
    right = xout[i] + bandwidth/2.0;
    NumericVector nbhd = y[x >= left & x <= right];
    s = std::accumulate(nbhd.begin(), nbhd.end(), 0.0);
    yout[i] = s/bandwidth;
  }
  return yout;
}

//' @title Kernel smoother for integer locations
//' @description Assumes missing locations have observation of 0
//' @param x Input locations (integers), sorted
//' @param y Input observations, should be length of x and non-missing
//' @param xout Output locations
//' @param bandwidth Bandwidth. Should be an odd integer
//' @param chunksize Chunksize to work in
//' @export
// [[Rcpp::export]]
NumericVector ksmooth_0_stitch(NumericVector x, NumericVector y,
                              NumericVector xout, double bandwidth, double chunksize){
  int k = xout.size();
  int nchunks = ceil(k/chunksize);
  int margin = ceil(bandwidth/2.0);
  IntegerVector l = (seq_len(nchunks)-1)*chunksize;
  NumericVector left_xout = xout[l];

  IntegerVector r = pmin(seq_len(nchunks)*chunksize-1, xout.size()-1);
  NumericVector right_xout = xout[r];
  NumericVector yout(k);
  for(size_t i = 0; i< nchunks; i ++){
  //parallel_for(blocked_range<size_t>(0, nchunks),
    //           [&](const blocked_range<size_t>& m){
    //             for(size_t i=m.begin(); i!=m.end(); ++i){
    NumericVector xo = xout[xout >= left_xout[i] & xout <= right_xout[i]];
    NumericVector xx = x[x >= left_xout[i]-margin & x <= right_xout[i] + margin];
    NumericVector yy = y[x >= left_xout[i]-margin & x <= right_xout[i] + margin];
    if(xx.size()==0){
      for(size_t j = l[i]; j <= r[i]; j ++){
        yout[j] = 0.0;
      }
    }else{
      NumericVector yo = ksmooth_0_cpp(xx, yy, xo, bandwidth);
      for(size_t j = 0; j < yo.size(); j ++ ){
        yout[l[i] + j] = yo[j];
      }
    }
  }
 //              });
  return yout;
}
