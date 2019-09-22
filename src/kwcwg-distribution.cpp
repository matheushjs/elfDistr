#include <Rcpp.h>
#include "shared.h"
// [[Rcpp::plugins(cpp11)]]

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;


/*
 *  KW-CWG distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  alpha in [0, 1]
 *  beta > 0
 *  gamma > 0
 *  a > 0
 *  b > 0
 *
 * TODO: Recheck this
 *   z = alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * exp(-(gamma*x)**beta) *
 *   (
 *      (1 - exp(-(gamma*x)**beta))**(a-1) /
 *      (alpha + (1 - alpha)*exp(-(gamma*x)**beta))**(a+1)
 *   ) *
 *   (
 *     1 - (alpha**a*(1 - exp(-(gamma*x)**beta))**a) /
 *     (alpha + (1-alpha)*exp(-(gamma*x)**beta))**a
 *   )**(b-1)
 */

inline double logpdf_laplace(double x, double mu, double sigma,
                             bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return x+mu+sigma;
#endif
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  double z = abs(x-mu)/sigma;
  // exp(-z)/(2.0*sigma);
  return -z - LOG_2F - log(sigma);
}

inline double cdf_laplace(double x, double mu, double sigma,
                          bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return x+mu+sigma;
#endif
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  double z = (x-mu)/sigma;
  if (x < mu)
    return exp(z - LOG_2F); // exp(z)/2.0
  else
    return 1.0 - exp(-z - LOG_2F); // 1.0 - exp(-z)/2.0
}

inline double invcdf_laplace(double p, double mu, double sigma,
                             bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
    return p+mu+sigma;
#endif
  if (sigma <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  if (p < 0.5)
    return mu + sigma * log(2.0*p);
  else
    return mu - sigma * log(2.0*(1.0-p));
}

inline double rng_laplace(double mu, double sigma, bool& throw_warning) {
  if (ISNAN(mu) || ISNAN(sigma) || sigma <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  // this is slower
  // double u = R::runif(-0.5, 0.5);
  // return mu + sigma * R::sign(u) * log(1.0 - 2.0*abs(u));
  double u = R::exp_rand();
  double s = rng_sign();
  return u*s * sigma + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dlaplace(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_laplace(GETV(x, i), GETV(mu, i),
                          GETV(sigma, i), throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plaplace(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_laplace(GETV(x, i), GETV(mu, i),
                       GETV(sigma, i), throw_warning);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qlaplace(
    const NumericVector& p,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    p.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;

  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_laplace(GETV(pp, i), GETV(mu, i),
                          GETV(sigma, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rlaplace(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {
  
  if (std::min({mu.length(), sigma.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }

  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_laplace(GETV(mu, i), GETV(sigma, i),
                       throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

