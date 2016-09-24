#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void ksmooth_0_cpp(NumericVector x,  NumericVector y, double bandwidth,
                NumericVector xout, NumericVector yout){
  //y.out <- sapply(x, FUN=function(xx){
  //  sum(y[ x <= (xx+bandwidth/2) & x >= (xx - bandwidth/2)])/(bandwidth+1)
  //})
  int i;
  int k = xout.size();
  int l;
  int j = x.size();
  //NumericVector yout(k);
  double s;
  double left;
  double right;

  for(i=0; i <= k; i++){
    left =xout[i] - bandwidth/2.0;
    right = xout[i] + bandwidth/2.0;
    s = 0;
    for(l = 0; l <= j; l++){
      if(x[l] >= left & x[l] <= right){
        s += y[l];
      }
      if(x[l] > right){
        break;
      }
    }
    yout[i] = s/bandwidth;
  }
}
