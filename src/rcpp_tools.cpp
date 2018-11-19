
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]


// E-Step
// [[Rcpp::export]]
NumericMatrix cGetEw(NumericVector x, double eta, double sig2) {
  NumericMatrix rst(x.size(), 2);

  double es  = pow(eta,2) + sig2;
  double tau = sqrt(sig2/es);
  double mu, ut, p, ew, ew2;
  int    i;
  
  for (i = 0; i < x.size(); i++) {
    mu  = eta * x[i] / es;
    ut  = - mu / tau;
    p   = R::dnorm(ut, 0, 1, 0) / (1 - R::pnorm(ut, 0, 1, 1, 0));
    ew  = mu + tau * p;
    ew2 = pow(tau, 2) * (1 + ut * p - pow(p,2)) + pow(ew, 2);

    rst(i,0) = ew;
    rst(i,1) = ew2;
  }

  return(rst);
}

// Simple Linear regression without intercept
// [[Rcpp::export]]
double cSlm(NumericVector y, NumericVector x) {
  int n = y.size(), i;
  double ysum = 0, xsum = 0;
  double xy = 0, xx = 0;
  double beta;

  for (i = 0; i < n; i++) {
    ysum += y[i];
    xsum += x[i];
  }

  ysum = ysum / n;
  xsum = xsum / n;

  for (i = 0; i < n; i++) {
    xy += (x[i] - xsum) * (y[i] - ysum);
    xx += pow(x[i] - xsum, 2);
  }

  beta = xy / xx;

  return(beta);
} 


// double cGetSnLpdf(NumericVector u, double eta, double mu, double sigma) {
//   double rst, e2, s2, e2ps2, r0, x1, x2;
//   int i;

//   e2    = pow(eta, 2);
//   s2    = pow(sigma, 2);
//   e2ps2 = sqrt(e2 + s2);
//   r0    = log(2) - log(e2ps2);

//   rst = 0;
//   for (i = 0; i < u.size(); i++) {
//     x1   = (u[i] - mu) / e2ps2;
//     x2   = (u[i] - mu) * eta / e2ps2 / sigma;
//     rst  += r0 + R::dnorm(x1, 0, 1, 1) + log(R::pnorm(x1, 0, 1, 1, 0));
//   }

//   return(rst);
// }

// EM Model 3 with Skewed Error
// [[Rcpp::export]]
List cEMMdl3(NumericVector init_pa, NumericVector init_ai, 
             NumericVector y, NumericVector x, NumericVector z, int unimodal, 
             int max_steps, double tol) {

  NumericVector pa = clone(init_pa);
  NumericVector ai = clone(init_ai);

  NumericVector last_pa(init_pa.size());
  NumericVector last_ai(init_pa.size());

  int           n = y.size();
  NumericVector ybz(n), yai(n);

  List          fitrst(2), rst(2);
  double        beta, mode, sig2, last_diff = 10000;
  int           inx = 0, i;
  double        tmp1;

  while(inx < max_steps & last_diff > tol) {
    last_pa = clone(pa);
    last_ai = clone(ai);

    for (i = 0; i < n; i++) {
      ybz[i] = y[i] - pa[0] * z[i];
    }

    fitrst = cUfit(ybz, x, unimodal);
    ai     = fitrst[0];
    if (1 == unimodal) {
      mode = as<double>(fitrst[1]);
    } else {
      mode = pa[2];
    }

    for (i = 0; i < n; i++) {
      yai[i] = y[i] - ai[i];
    }

    beta = cSlm(yai, z);

    tmp1 = 0;
    for (i = 0; i < n; i++) {
      tmp1 += pow(yai[i] - beta * z[i], 2);
    }
    sig2 = tmp1 / n;

    pa[0] = beta;
    pa[1] = sig2;
    pa[2] = mode;

    last_diff = 0;
    for (i = 0; i < pa.size(); i++) {
      last_diff = fmax(last_diff, fabs(pa[i] - last_pa[i]));
    }

    for (i = 0; i < ai.size(); i++) {
      last_diff = fmax(last_diff, fabs(ai[i] - last_ai[i]));
    }

    inx++;
  }

  if (last_diff > tol) {
    rst = List::create(NA_REAL, NA_REAL);
  } else {
    rst = List::create(_["mle.pa"] = pa, _["mle.ai"] = ai);
  }

  return(rst);
}

// EM Model 4 with Skewed Error
// [[Rcpp::export]]
List cEMMdl4(NumericVector init_pa, NumericVector init_ai, NumericVector init_ci,
             NumericVector y, NumericVector x, NumericVector z, int unimodal, 
             int max_steps, double tol) {

  NumericVector pa = clone(init_pa);
  NumericVector ai = clone(init_ai);
  NumericVector ci = clone(init_ci);

  NumericVector last_pa(init_pa.size());
  NumericVector last_ai(init_pa.size());

  int           n = y.size();
  NumericVector yaieta(n), yzeta(n);
  NumericMatrix ew(n, 2);

  List          fitrst(2), rst(2);
  double        beta, mode, eta, sig2, last_diff = 10000;
  double        tmp1, tmp2;
  int           inx = 0, i;

  while(inx < max_steps & last_diff > tol) {
    last_pa = clone(pa);
    last_ai = clone(ai);

    ew = cGetEw(ci, pa[1], pa[2]);

    for (i = 0; i < n; i++) {
      yaieta[i] = y[i] - ai[i] - pa[1] * ew(i,0);
    }

    beta = cSlm(yaieta, z);

    for (i = 0; i < n; i++) {
      yzeta[i] = y[i] - beta * z[i] - pa[1] * ew(i,0);
    }

    fitrst = cUfit(yzeta, x, unimodal);
    ai     = fitrst[0];
    if (1 == unimodal) {
      mode = as<double>(fitrst[1]);
    } else {
      mode = pa[3];
    }

    for (i = 0; i < n; i++) {
      ci[i] = y[i] - ai[i] - beta * z[i];
    }

    tmp1 = 0;
    tmp2 = 0;
    for (i = 0; i < n; i++) {
      tmp1 += ci[i] * ew(i,0);
      tmp2 += ew(i,1);
    }
    eta = tmp1 / tmp2;


    tmp1 = 0;
    for (i = 0; i < n; i++) {
      tmp1 += pow(ci[i], 2) + pow(eta, 2) * ew(i,1) - 2 * eta * ci[i] * ew(i,0);
    }
    sig2 = tmp1 / n;

    pa[0] = beta;
    pa[1] = eta;
    pa[2] = sig2;
    pa[3] = mode;


    last_diff = 0;
    for (i = 0; i < pa.size(); i++) {
      last_diff = fmax(last_diff, fabs(pa[i] - last_pa[i]));
    }

    for (i = 0; i < ai.size(); i++) {
      last_diff = fmax(last_diff, fabs(ai[i] - last_ai[i]));
    }

    inx++;
  }

  if (last_diff > tol) {
    rst = List::create(NA_REAL, NA_REAL);
  } else {
    rst = List::create(_["mle.pa"] = pa, _["mle.ai"] = ai);
  }

  return(rst);
}
