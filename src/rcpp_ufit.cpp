
#include <Rcpp.h>
#include <Rmath.h>

extern "C" {
  void hello_();
}

extern "C" {
  void pava_(double *y, double *w, double *kt, const int n); 
}

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// PAVA algorithm
// [[Rcpp::export]]
NumericVector cPava0(NumericVector y, NumericVector w) {

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
NumericVector cPava(NumericVector y, NumericVector w) {
  int n = y.size();
  NumericVector yy = clone(y), ww = clone(w);
  NumericVector kt(n);

  pava_(yy.begin(), ww.begin(), kt.begin(), n);
  return(yy);
}


// UFIT
// [[Rcpp::export]]
List cUfit(NumericVector y, NumericVector x, int unimodal) {
  List rst(2);
  NumericVector ai;
  double mode;

  Function uf("get.isoreg.ufit");

  if (0 == unimodal) {
    NumericVector w(y.size(), 1);
    ai   = cPava(y, w);
    mode = NA_REAL;
    rst  = List::create(_["ai"] = ai, _["mode"] = mode);
  } else {
    rst = uf(_["y"] = y, _["x"] = x);
  }

  return(rst);
}

