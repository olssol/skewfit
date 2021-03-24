#include <Rcpp.h>

using namespace Rcpp;

List cIso(NumericVector y, int unimodal);

// [[Rcpp::plugins("cpp11")]]

// rank x in a sorted ascending vector
// [[Rcpp::export]]
int rank_x(double x, NumericVector v) {
  int n = v.size();
  int l = 0, u = n - 1, cur;
  int rst;

  if (x <= v[l]) {
    rst = l;
  } else if (x >= v[u]) {
    rst = u;
  } else {
    while (l < u - 1) {
      cur = (l + u)/2;
      if (x < v[cur]) {
        u = cur;
      } else if (x >= v[cur]) {
        l = cur;
      }
    }
    rst = l;
  }

  return(rst);
}

// approximate pnorm from Aludaat and Alodat (2008)
// [[Rcpp::export]]
double apnorm(double u) {
  int    sgn;
  double rst;
  if (u >= 0) {
    sgn = 1;
  } else {
    sgn = -1;
  }

  rst = 0.5 + sgn * 0.5 * sqrt(1 - exp(-sqrt(M_PI / 8) * pow(u, 2)));
  return(rst);
}

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

    if (ut < -4) {
      ut = -4.0;
    } else if (ut > 7) {
      ut = 7.0;
    }

    p   = R::dnorm(ut, 0, 1, 0) / (1 - R::pnorm(ut, 0, 1, 1, 0));
    ew  = mu + tau * p;
    ew2 = pow(tau, 2) * (1 + ut * p - pow(p,2)) + pow(ew, 2);

    rst(i,0) = ew;
    rst(i,1) = ew2;
  }

  return(rst);
}

//Simple linear regression without intercept
// [[Rcpp::export]]
double cSlm0(NumericVector y, NumericVector x) {
  int    n = y.size(), i;
  double ysum = 0, xsum = 0;
  double xy = 0, xx = 0;

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

  return(xy / xx);
}

// [[Rcpp::export]]
double cSlm1(NumericVector y, NumericVector x) {
  NumericVector beta(1);
  Function coef("get_lm_coeff");

  beta = coef(y, x);
  return(beta(0));
}


// Simple Linear regression without intercept
// [[Rcpp::export]]
NumericVector cSlm(NumericVector y, NumericMatrix x) {

  int np = x.ncol();
  NumericVector beta(np);

  Function coef("get_lm_coeff");
  beta = coef(y, x);

  // if (1 == np) {
  //   beta[0] = cSlm0(y, x(_,1));
  // } else {
  //   Function coef("get.lm.coeff");
  //   beta = coef(y, x);
  // }

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

//' Parametric alpha(x) with Skewed Error
//'
//' alpha(x) =  alpha * x where alpha > 0
//' @export
// [[Rcpp::export]]
List fit_para_skew(NumericVector init_pa, NumericVector y, NumericVector x, NumericMatrix z,
                   int max_steps, double tol) {

  NumericVector pa = clone(init_pa);
  NumericVector last_pa(pa.size());
  int           n  = y.size(), nz = z.ncol();
  NumericVector yeta(n), yzeta(n), beta(1+nz), alpha(1);
  NumericMatrix ew(n, 2);
  NumericVector ai(n), ci(n);

  Function      coef("get_lm_coeff_2");
  Function      coef0("get_lm_coeff");

  List          rst(2);
  double        eta, sig2, last_diff = 10000;
  double        tmp1, tmp2;
  int           inx = 0, i, j;


  //initial ci
  for (i = 0; i < n; i++) {
    tmp1 = pa[0] + pa[1] * x[i];
    for (j = 0; j < nz; j++) {
      tmp1 += pa[2+j] * z(i,j);
    }

    ci[i] = y[i] - tmp1;
  }

  while (inx < max_steps & last_diff > tol) {
    last_pa = clone(pa);
    ew      = cGetEw(ci, pa[2+nz], pa[3+nz]);

    for (i = 0; i < n; i++) {
      yeta[i] = y[i] - pa[1] * x[i] - pa[2+nz] * ew(i,0);
    }

    // get intercept and coeff for z
    beta = coef(yeta, z);

    //get alpha
    for (i = 0; i < n; i++) {
      tmp1 = beta[0];
      for (j = 0; j < nz; j++) {
        tmp1 += beta[1+j] * z(i,j);
      }

      yzeta[i] = y[i] - tmp1 - pa[2+nz] * ew(i,0);
    }

    //Rcout << yzeta << "\n";
    alpha = coef0(yzeta, x);
    if (alpha[0] < 0)
      alpha[0] = 0;

    //update ci
    for (i = 0; i < n; i++) {
      ai[i] = beta[0] + alpha[0] * x[i];
      tmp1  = ai[i];
      for (j = 0; j < nz; j++) {
        tmp1 += beta[1+j] * z(i,j);
      }
      ci[i] = y[i] - tmp1;
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

    //return parameter
    pa[0] = beta[0];
    pa[1] = alpha[0];
    for (i = 1; i < beta.size(); i++) {
      pa[i+1] = beta[i];
    }
    pa[2+nz] = eta;
    pa[3+nz] = sig2;

    last_diff = 0;
    for (i = 0; i < pa.size(); i++) {
      last_diff = fmax(last_diff, fabs(pa[i] - last_pa[i]));
    }

    inx++;
  }

  if (last_diff > tol) {
    rst = List::create(NA_REAL, NA_REAL);
  } else {
    rst = List::create(_["mle_pa"]   = pa,
                       _["mle_ai"]   = ai,
                       _["residual"] = ci);
  }

  return(rst);
}

//' Isotonic regression with Normal Error
//'
//' @export
// [[Rcpp::export]]
List fit_iso_norm(NumericVector init_pa, NumericVector init_ai,
                  NumericVector y, NumericVector x, NumericMatrix z, int unimodal,
                  int max_steps, double tol) {

  NumericVector pa = clone(init_pa);
  NumericVector ai = clone(init_ai);

  NumericVector last_pa(init_pa.size());
  NumericVector last_ai(init_pa.size());

  int           n  = y.size();
  int           nz = z.ncol();
  NumericVector ybz(n), yai(n), beta(nz);
  NumericVector residual(n);

  List          fitrst(2), rst(2);
  double        mode, sig2, last_diff = 10000;
  int           inx = 0, i, j;
  double        tmp1, tmp2;

  while(inx < max_steps & last_diff > tol) {
    last_pa = clone(pa);
    last_ai = clone(ai);

    for (i = 0; i < n; i++) {
      tmp2 = 0;
      for (j = 0; j < nz; j++) {
        tmp2 += pa[j] * z(i,j);
      }
      ybz[i] = y[i] - tmp2;
    }

    fitrst = cIso(ybz, unimodal);
    ai     = fitrst[0];
    if (1 == unimodal) {
      mode = x[as<int>(fitrst[1])-1];
    } else {
      mode = pa[nz+1];
    }

    for (i = 0; i < n; i++) {
      yai[i] = y[i] - ai[i];
    }

    beta = cSlm(yai, z);
    tmp1 = 0;
    for (i = 0; i < n; i++) {
      tmp2 = 0;
      for (j = 0; j < nz; j++) {
        tmp2 += beta[j] * z(i,j);
      }
      residual[i] = yai[i] - tmp2;

      tmp1 += pow(residual[i], 2);
    }
    sig2 = tmp1 / n;

    //return parameter
    for (i = 0; i < nz; i++) {
      pa[i] = beta[i];
    }
    pa[nz]   = sig2;
    pa[nz+1] = mode;

    last_diff = 0;
    for (i = 0; i < nz+2; i++) {
      last_diff = fmax(last_diff, fabs(pa[i] - last_pa[i]));
    }

    for (i = 0; i < ai.size(); i++) {
      last_diff = fmax(last_diff, fabs(ai[i] - last_ai[i]));
    }

    inx++;
  }

  if (last_diff > tol) {
    rst = List::create(NA_REAL, NA_REAL, NA_REAL);
  } else {
    rst = List::create(_["mle_pa"]   = pa,
                       _["mle_ai"]   = ai,
                       _["residual"] = residual
                       );
  }

  return(rst);
}

// cEMMdl4
//' Isotonic regression with Skewed Error
//'
//' @export
// [[Rcpp::export]]
List fit_iso_skew(NumericVector pa, NumericVector ai,
                  NumericVector y, NumericVector x, NumericMatrix z, int unimodal, int usez,
                  int max_steps, double tol) {

  NumericVector last_pa(pa.size());
  NumericVector last_ai(ai.size());

  int           n  = y.size(), nz = z.ncol();
  NumericVector yaieta(n), yzeta(n), beta(nz);
  NumericMatrix ew(n, 2);

  List          fitrst(2), rst(2);
  double        mode, eta, sig2, last_diff = 10000;
  double        tmp1, tmp2;
  int           inx = 0, i, j;

  if (0 == usez) {
    nz = 0;
  }

  //initial ci
  NumericVector ci(n);
  for (i = 0; i < n; i++) {
    tmp1 = 0;
    for (j = 0; j < nz; j++) {
      tmp1 += pa[j] * z(i,j);
    }

    ci[i] = y[i] - ai[i] - tmp1;
  }

  while(inx < max_steps & last_diff > tol) {
    last_pa = clone(pa);
    last_ai = clone(ai);

    ew = cGetEw(ci, pa[nz], pa[nz+1]);
    for (i = 0; i < n; i++) {
      yaieta[i] = y[i] - ai[i] - pa[nz] * ew(i,0);
    }

    if (0 < nz) {
      beta = cSlm(yaieta, z);
    }

    for (i = 0; i < n; i++) {
      tmp1 = 0;
      for (j = 0; j < nz; j++) {
        tmp1 += beta[j] * z(i,j);
        }

      yzeta[i] = y[i] - tmp1 - pa[nz] * ew(i,0);
    }

    fitrst = cIso(yzeta, unimodal);
    ai     = fitrst[0];
    if (1 == unimodal) {
      mode = x[as<int>(fitrst[1])-1];
    } else {
      mode = pa[nz+2];
    }

    for (i = 0; i < n; i++) {
      tmp1 = 0;
      for (j = 0; j < nz; j++) {
        tmp1 += beta[j] * z(i,j);
      }
      ci[i] = y[i] - ai[i] - tmp1;
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

    //return parameter
    for (i = 0; i < nz; i++) {
      pa[i] = beta[i];
    }

    pa[nz]   = eta;
    pa[nz+1] = sig2;
    pa[nz+2] = mode;

    last_diff = 0;
    for (i = 0; i < nz+3; i++) {
      last_diff = fmax(last_diff, fabs(pa[i] - last_pa[i]));
    }

    for (i = 0; i < ai.size(); i++) {
      last_diff = fmax(last_diff, fabs(ai[i] - last_ai[i]));
    }

    inx++;
  }

  if (last_diff > tol) {
    rst = List::create(NA_REAL, NA_REAL, NA_REAL, NA_REAL);
  } else {
    rst = List::create(_["mle_pa"]   = pa,
                       _["mle_ai"]   = ai,
                       _["residual"] = ci);
  }

  return(rst);
}


//' Kernel Function
//'
//' @export
// [[Rcpp::export]]
double  get_kernel(double v) {
  double rst;

  if (abs(v) <= 1) {
    rst  = 16.0 + 35 * v - 35 * pow(v, 3);
    rst += 21 * pow(v, 5) - 5 * pow(v, 7);
    rst /= 32.0;
  } else if (v > 1) {
        rst = 1;
  } else {
        rst = 0;
  }

  return rst;
}

//' Kernel smooth function
//'
//'
//'
// [[Rcpp::export]]
NumericMatrix get_kernel_fn(NumericVector x, NumericMatrix fn, double h,
                            bool correction = true) {

  int n = x.size(), nf = fn.nrow();
  int i, j, njumps;
  double t, kt, cf;
  double a, b, min_f, range_f;

  NumericMatrix jumps(nf, 2);
  NumericMatrix rst(n, 2);

  // minimum and maximum
  min_f   = fn(0,      1);
  range_f = fn(nf - 1, 1) - fn(0, 1);

  // find jumps in fn
  njumps           = 0;
  jumps(njumps, 0) = fn(0, 0);
  jumps(njumps, 1) = 0;

  for (i = 1; i < nf; i++) {
    if (fn(i, 1) == fn(i - 1, 1))
      continue;

    njumps++;
    jumps(njumps, 0) = fn(i, 0);
    jumps(njumps, 1) = fn(i, 1) - fn(i - 1, 1);
    jumps(njumps, 1) /= range_f;
  }

  // boundary correction
  a = 0;
  b = 1; // jumps(njumps, 0)

  for (i = 0; i < n; i++) {
    cf = 0;
    for (j = 0; j <= njumps; j++) {
      t  = (x[i] - jumps(j, 0)) / h;
      kt = get_kernel(t);

      if (correction) {
        t   = (x[i] + jumps(j, 0) - 2 * a) / h;
        kt += get_kernel(t);

        t   = (2 * b - x[i] - jumps(j, 0)) / h;
        kt -= get_kernel(t);
      }

      cf += kt * jumps(j, 1);
    }

    rst(i, 0) = x[i];
    rst(i, 1) = cf * range_f + min_f;
  }

  //return
  return(rst);
}


//' Isotonic regression prediction
//'
//' @export
// [[Rcpp::export]]
NumericMatrix pred_iso(NumericVector x, NumericMatrix iso_fit,
                       double h = -1,
                       bool correction = true) {

  int           nx = x.size();
  NumericMatrix rst(nx, 2);

  int i, j;
  if (h <= 0) {
    for (i = 0; i < nx; i++) {
      j         = rank_x(x[i], iso_fit(_, 0));
      rst(i, 0) = x[i];
      rst(i, 1) = iso_fit(j, 1);

      // Rcout << i << j << x[i] << std::endl;
    }
  } else {
    rst = get_kernel_fn(x, iso_fit, h, correction);
  }

  return(rst);
}
