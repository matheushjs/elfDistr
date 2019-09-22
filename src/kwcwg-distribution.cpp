#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;

#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))


/*
 *  KW-CWG distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  alpha in [0, 1]
 *  beta >= 0
 *  gamma >= 0
 *  a >= 0
 *  b >= 0
 *
 * Recheck this
 *   z = alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * exp(-(gamma*x)**beta) *
 *   (
 *      (1 - exp(-(gamma*x)**beta))**(a-1) /
 *      (alpha + (1 - alpha)*exp(-(gamma*x)**beta))**(a+1)
 *   ) *
 *   (
 *     1 - (alpha**a*(1 - exp(-(gamma*x)**beta))**a) /
 *         (alpha + (1-alpha)*exp(-(gamma*x)**beta))**a
 *   )**(b-1)
 */

inline double logpdf_kwcwg(
	double x, double alpha, double beta,
	double gamma, double a, double b, bool &throw_warning)
{
#ifdef IEEE_754
	if(ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(gamma) || ISNAN(a) || ISNAN(b))
		return p+alpha+beta+gamma+a+b;
#endif

	if(alpha < 0.0 || alpha > 1.0
	   || beta < 0.0
	   || gamma < 0.0
	   || a < 0.0
	   || b < 0.0)
	{
		throw_warning = true;
		return NAN;
	}

	// Common term in the equation
	double aux1 = exp(-(gamma*x)**beta);

	// Here we will factor f(x) as being A * (B^(a-1)/C^(a-1)) * (1 - D/E)^(b-1)
	double A = pow(alpha,a) * beta * gamma * a * b * pow(gamma*x,beta-1) * aux1;
	double B = 1 - aux1;
	double C = alpha + (1-alpha)*aux1;
	double D = pow(alpha,a) * pow(1-aux1,a);
	double E = pow(alpha + (1-alpha)*aux1,a);

	// return A * (B**(a-1)/C**(a+1)) * (1 - D/E)**(b-1)

	return log(A) + (a-1)*log(B) - (a+1)*log(C) + (b-1)*log(1 - D/E);
}

inline double cdf_kwcwg(
	double x, double alpha, double beta,
	double gamma, double a, double b, bool &throw_warning)
{
#ifdef IEEE_754
	if(ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(gamma) || ISNAN(a) || ISNAN(b))
		return p+alpha+beta+gamma+a+b;
#endif

	if(alpha < 0.0 || alpha > 1.0
	   || beta < 0.0
	   || gamma < 0.0
	   || a < 0.0
	   || b < 0.0)
	{
		throw_warning = true;
		return NAN;
	}

	return 1-
		pow(1-pow(alpha,a)*
			pow(
				(1-exp(-pow(gamma*x,beta))) /
				(alpha + (1-alpha)*exp(-pow(gamma*x,beta)))
			,a)
		,b);
}

inline double invcdf_kwcwg(
	double p, double alpha, double beta,
	double gamma, double a, double b, bool &throw_warning)
{
#ifdef IEEE_754
	if(ISNAN(p) || ISNAN(alpha) || ISNAN(beta) || ISNAN(gamma) || ISNAN(a) || ISNAN(b))
		return p+alpha+beta+gamma+a+b;
#endif

	if(alpha < 0.0 || alpha > 1.0
	   || beta < 0.0
	   || gamma < 0.0
	   || a < 0.0
	   || b < 0.0
	   || !VALID_PROB(p))
	{
		throw_warning = true;
		return NAN;
	}

	// Common term
	double aux = pow(1-pow(1-p,1/b),1/a);

	return pow(log((alpha + (1-alpha)*aux)/(alpha*(1-aux))),1/beta)/gamma;
}

// Random number generation
inline double rng_kwcwg(
	double alpha, double beta, double gamma,
	double a, double b, bool &throw_warning)
{
	if(ISNAN(p) || ISNAN(alpha) || ISNAN(beta) || ISNAN(gamma) || ISNAN(a) || ISNAN(b)
	   || alpha < 0.0 || alpha > 1.0
	   || beta < 0.0
	   || gamma < 0.0
	   || a < 0.0
	   || b < 0.0
	{
		throw_warning = true;
		return NAN;
	}

	return invcdf_kwcwg(R::runif(0, 1), alpha, beta, gamma, a, b, throw_warning);
}

// [[Rcpp::export]]
NumericVector cpp_dkwcwg(
	const NumericVector& x,
	const NumericVector& alpha,
	const NumericVector& beta,
	const NumericVector& gamma,
	const NumericVector& a,
	const NumericVector& b,
	const bool& log_prob = false
){
	
	if(std::min(
		{x.length(), alpha.length(), beta.length(),
		 gamma.length(), a.length(), b.length()}) < 1)
	{
		return NumericVector(0);
	}

	int maxN = std::max({
		x.length(),
		alpha.length(),
		beta.length(),
		gamma.length(),
		a.length(),
		b.length()
	});
	NumericVector p(maxN);
	
	bool throw_warning = false;

	for(int i = 0; i < maxN; i++)
		p[i] = logpdf_kwcwg(
			GETV(x, i),
			GETV(alpha, i),
			GETV(beta, i),
			GETV(gamma, i),
			GETV(a, i),
			GETV(b, i),
			throw_warning);

	if(!log_prob)
		p = Rcpp::exp(p);
	
	if(throw_warning)
		Rcpp::warning("NaNs produced");

	return p;
}


// [[Rcpp::export]]
NumericVector cpp_pkwcwg(
	const NumericVector& x,
	const NumericVector& alpha,
	const NumericVector& beta,
	const NumericVector& gamma,
	const NumericVector& a,
	const NumericVector& b,
	const bool& lower_tail = true,
	const bool& log_prob = false
){

	if(std::min(
		{x.length(), alpha.length(), beta.length(),
		 gamma.length(), a.length(), b.length()}) < 1)
	{
		return NumericVector(0);
	}

	int maxN = std::max({
		x.length(),
		alpha.length(),
		beta.length(),
		gamma.length(),
		a.length(),
		b.length()
	});
	NumericVector p(maxN);

	bool throw_warning = false;

	for (int i = 0; i < maxN; i++)
		p[i] = cdf_kwcwg(
			GETV(x, i),
			GETV(alpha, i),
			GETV(beta, i),
			GETV(gamma, i),
			GETV(a, i),
			GETV(b, i),
			throw_warning);

	if (!lower_tail)
		p = 1.0 - p;
	
	if (log_prob)
		p = Rcpp::log(p);
	
	if (throw_warning)
		Rcpp::warning("NaNs produced");

	return p;
}


// [[Rcpp::export]]
NumericVector cpp_qkwcwg(
	const NumericVector& p,
	const NumericVector& alpha,
	const NumericVector& beta,
	const NumericVector& gamma,
	const NumericVector& a,
	const NumericVector& b,
	const bool& lower_tail = true,
	const bool& log_prob = false
){
	if(std::min(
		{x.length(), alpha.length(), beta.length(),
		 gamma.length(), a.length(), b.length()}) < 1)
	{
		return NumericVector(0);
	}

	int maxN = std::max({
		x.length(),
		alpha.length(),
		beta.length(),
		gamma.length(),
		a.length(),
		b.length()
	});
	NumericVector q(maxN);
	NumericVector pp = Rcpp::clone(p);
	
	bool throw_warning = false;

	if (log_prob)
		pp = Rcpp::exp(pp);
	
	if (!lower_tail)
		pp = 1.0 - pp;

	for (int i = 0; i < maxN; i++)
		q[i] = invcdf_kwcwg(
			GETV(pp, i),
			GETV(alpha, i),
			GETV(beta, i),
			GETV(gamma, i),
			GETV(a, i),
			GETV(b, i),
			throw_warning);
	
	if (throw_warning)
		Rcpp::warning("NaNs produced");

	return q;
}

// [[Rcpp::export]]
NumericVector cpp_rkwcwg(
	const int& n,
	const NumericVector& alpha,
	const NumericVector& beta,
	const NumericVector& gamma,
	const NumericVector& a,
	const NumericVector& b
){
	if(std::min(
		{x.length(), alpha.length(), beta.length(),
		 gamma.length(), a.length(), b.length()}) < 1)
	{
		Rcpp::warning("NAs produced");
		return NumericVector(n, NA_REAL);
	}

	NumericVector x(n);
	
	bool throw_warning = false;

	for (int i = 0; i < n; i++)
		x[i] = rng_kwcwg(
			GETV(alpha, i),
			GETV(beta, i),
			GETV(gamma, i),
			GETV(a, i),
			GETV(b, i),
			throw_warning);
	
	if (throw_warning)
		Rcpp::warning("NAs produced");

	return x;
}

