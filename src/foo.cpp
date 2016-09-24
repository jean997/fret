#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
void double_me3(NumericVector x, NumericVector y) {
  int k = x.size();
  int i;
  // takes a numeric input and doubles it
  for(i=0; i <= k; i++){
    y[i] = 2*x[i];
  }
}
