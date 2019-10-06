#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;
using Rcpp::Rcout;

#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))


/*
 *  Kw-CWG distribution
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
 *  z = alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * exp(-(gamma*x)**beta) *
 *  (
 *     (1 - exp(-(gamma*x)**beta))**(a-1) /
 *     (alpha + (1 - alpha)*exp(-(gamma*x)**beta))**(a+1)
 *  ) *
 *  (
 *    1 - (alpha**a*(1 - exp(-(gamma*x)**beta))**a) /
 *        (alpha + (1-alpha)*exp(-(gamma*x)**beta))**a
 *  )**(b-1)
 */

inline double logpdf_kwcwg(
	double x, double alpha, double beta,
	double gamma, double a, double b, bool &throw_warning)
{
#ifdef IEEE_754
	if(ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(gamma) || ISNAN(a) || ISNAN(b))
		return x+alpha+beta+gamma+a+b;
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
	double aux1 = exp(-pow(gamma*x,beta));

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
		return x+alpha+beta+gamma+a+b;
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
	if(ISNAN(alpha) || ISNAN(beta) || ISNAN(gamma) || ISNAN(a) || ISNAN(b)
	   || alpha < 0.0 || alpha > 1.0
	   || beta < 0.0
	   || gamma < 0.0
	   || a < 0.0
	   || b < 0.0)
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
	const unsigned int xL = x.length();
	const unsigned int alphaL = alpha.length();
	const unsigned int betaL = beta.length();
	const unsigned int gammaL = gamma.length();
	const unsigned int aL = a.length();
	const unsigned int bL = b.length();
	
	if(std::min({ xL, alphaL, betaL, gammaL, aL, bL }) < 1)
		return NumericVector(0);

	int maxN = std::max({ xL, alphaL, betaL, gammaL, aL, bL });
	NumericVector p = Rcpp::no_init(maxN);
	
	bool throw_warning = false;

	const int STRIDE = 8;
	int lim = maxN / STRIDE;

	#pragma omp parallel for
	for(int i = 0; i < lim*STRIDE; i += STRIDE){
		// Read parameters
		double xX[STRIDE];
		double alphaX[STRIDE];
		double betaX[STRIDE];
		double gammaX[STRIDE];
		double aX[STRIDE];
		double bX[STRIDE];

		// Common result
		double aux[STRIDE];

		// Intermediate results
		double A[STRIDE];
		double B[STRIDE];
		double C[STRIDE];
		double D[STRIDE];
		double E[STRIDE];

		#pragma omp simd
		for(int j = 0; j < STRIDE; j++){
			const int idx = i + j;
			xX[j]     = x[idx%xL];
			alphaX[j] = alpha[idx%alphaL];
			betaX[j]  = beta[idx%betaL];
			gammaX[j] = gamma[idx%gammaL];
			aX[j]     = a[idx%aL];
			bX[j]     = b[idx%bL];
		}

		#pragma omp simd
		for(int j = 0; j < STRIDE; j++){
			const int idx = i + j;
			aux[j] = exp(-pow(gammaX[j]*xX[j],betaX[j]));
		}

		#pragma omp simd
		for(int j = 0; j < STRIDE; j++){
			const int idx = i + j;
			A[j] = pow(alphaX[j],aX[j]) * betaX[j] * gammaX[j] * aX[j] * bX[j] * pow(gammaX[j]*xX[j],betaX[j]-1) * aux[j];
		}

		#pragma omp simd
		for(int j = 0; j < STRIDE; j++){
			const int idx = i + j;
			B[j] = 1 - aux[j];
		}

		#pragma omp simd
		for(int j = 0; j < STRIDE; j++){
			const int idx = i + j;
			C[j] = alphaX[j] + (1-alphaX[j])*aux[j];
		}

		#pragma omp simd
		for(int j = 0; j < STRIDE; j++){
			const int idx = i + j;
			D[j] = pow(alphaX[j],aX[j]) * pow(1-aux[j],aX[j]);
		}

		#pragma omp simd
		for(int j = 0; j < STRIDE; j++){
			const int idx = i + j;
			E[j] = pow(alphaX[j] + (1-alphaX[j])*aux[j],aX[j]);
		}

		#pragma omp simd
		for(int j = 0; j < STRIDE; j++){
			const int idx = i + j;
			p[idx] = log(A[j]) + (aX[j]-1)*log(B[j]) - (aX[j]+1)*log(C[j]) + (bX[j]-1)*log(1 - D[j]/E[j]);
		}

		for(int j = 0; j < STRIDE; j++){
			const int idx = i + j;
			#ifdef IEEE_754
			if(ISNAN(xX[j]) || ISNAN(alphaX[j]) || ISNAN(betaX[j]) || ISNAN(gammaX[j]) || ISNAN(aX[j]) || ISNAN(bX[j])){
				throw_warning = true;
				p[idx] = NAN;
			} else // XXX: ATTENTION HERE!!
			#endif

			if(alphaX[j] < 0.0 || alphaX[j] > 1.0
			   || betaX[j] < 0.0
			   || gammaX[j] < 0.0
			   || aX[j] < 0.0
			   || bX[j] < 0.0)
			{
				throw_warning = true;
				p[idx] = NAN;
			}
		}
	}

	for(int i = lim*STRIDE; i < maxN; i++)
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
		{p.length(), alpha.length(), beta.length(),
		 gamma.length(), a.length(), b.length()}) < 1)
	{
		return NumericVector(0);
	}

	int maxN = std::max({
		p.length(),
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
		{alpha.length(), beta.length(),
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

