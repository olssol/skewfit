
#include <Rcpp.h>
#include <Rmath.h>

extern "C" {
  void pava_(double *y, double *w, double *kt, int *n); 
}

extern "C" {
  void ufit_(double *xk, double *wk, double *xmode,
             double *x,  double *w,  double *mse,
             double *x1, double *w1, double *x2,
             double *w2, double *ind, double *kt,
             int *n, int *goof);
}

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// // PAVA algorithm
// [[Rcpp::export]]
NumericVector cPava(NumericVector y, NumericVector w) {

  int n = y.size();
  std::vector<int> kt(n);

  int    i, j, k1, k2, same = 0;
  int    wnew;
  double ynew;

  if (1 == n)
    return(y);

  for (i = 0; i < n; i++) {
    kt[i] = i;
  }

  while (0 == same) {
    same = 1;
    for (i = 1; i < n; i++) {
      if (y[i-1] <= y[i]) 
        continue;

      k1 = kt[i];
      k2 = kt[i-1];

      for (j = 0; j < n; j++) {
        if (kt[j] == k1) {
          kt[j] = k2;
        }
      }

      wnew = w[i-1] + w[i];
      ynew = (w[i-1] * y[i-1] + w[i] * y[i])/wnew;

      for (j = 0; j < n; j++) {
        if (k2 == kt[j]) {
          y[j] = ynew;
          w[j]   = wnew;
        }
      }
      same = 0;
    }
  }

  return(y);
}

// [[Rcpp::export]]
NumericVector FPava(NumericVector y, NumericVector w) {
  int n = y.size();
  NumericVector kt(n);
  pava_(y.begin(), w.begin(), kt.begin(), &n);
  return(y);
}


// [[Rcpp::export]]
List FUfit(NumericVector y, NumericVector w) {
  int n = y.size();
  NumericVector x0(n), w0(n), x1(n), w1(n), x2(n), w2(n), ind(n), kt(n);
  double xmode = -1, mse;
  int goof = 1;

  ufit_(y.begin(),  w.begin(),  &xmode, x0.begin(), w0.begin(), &mse,
        x1.begin(), w1.begin(), x2.begin(), w2.begin(),
        ind.begin(), kt.begin(), &n, &goof);

  List rst = List::create(_["ai"] = x0,
                          _["mode"] = xmode);
  return(rst);
}

// Iso functions 
// [[Rcpp::export]]
List cIso(NumericVector y, int unimodal) {
  List rst(2);
  NumericVector ai;
  double mode;

  Function uf("get.isoreg.ufit");
  NumericVector w(y.size(), 1);

  if (0 == unimodal) {
    ai   = cPava(y, w);
    mode = NA_REAL;
    rst  = List::create(_["ai"] = ai, _["mode"] = mode);
  } else {
    rst  = FUfit(y, w);
  }

  return(rst);
}

