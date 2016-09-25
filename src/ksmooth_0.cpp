#include <string>
#include <vector>
#include <algorithm>


using namespace std;

void ksmooth_0_cpp(double *x,  double *y, double bandwidth, double *xout, double *yout, int lx, int lxout){

  //y.out <- sapply(x, FUN=function(xx){
  //  sum(y[ x <= (xx+bandwidth/2) & x >= (xx - bandwidth/2)])/(bandwidth+1)
  //})
  int i;
  int l;


  double s;
  double left;
  double right;

  for(i=0; i <= lxout; i++){
    left =xout[i] - bandwidth/2.0;
    right = xout[i] + bandwidth/2.0;
    s = 0;
    for(l = 0; l <= lx; l++){
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

//
//
// R wrappers
extern "C" {
  void ksmooth_0_cpp_R(double *x, double *y,  double *bandwidth, double *xout, double *yout, int *lx, int *lxout) {
    ksmooth_0_cpp(x, y, *bandwidth, xout, yout, *lx, *lxout);
  }
}
