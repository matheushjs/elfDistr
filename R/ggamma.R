
#' Generalized Gamma Probability Distribution
#'
#' Fast implementation of density, distribution function, quantile function
#' and random generation for the Generalized Gamma probability distribution.
#'
#' @param x,q	          vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param a,b,k	          Parameters of the distribution, all of which must be positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	  logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' The generalized gamma distribution proposed by Stacy (1962) has parameters
#' \eqn{a, d, p}, but here we adopt the reparametrization
#' \deqn{
#'   a = a
#' }
#' \deqn{
#'   b = p
#' }
#' \deqn{
#'   k = \frac{d}{p}
#' }
#' as is used by the R package *flexsurv*.
#'
#' Probability density function
#' \deqn{
#'    f(x) = \frac{b x^{bk-1} \exp[-(x/a)^b]}{a^{bk} \Gamma(k)}
#' }
#'
#' Cumulative density function
#' \deqn{
#'    F(x) = \frac{\gamma(k, (x/a)^p)}{\Gamma(k)}
#' }
#'
#' The above function can be written in terms of a \eqn{Gamma(\alpha, \beta)}.
#' Let \eqn{T \sim Gamma(k, 1)} and its cumulative distribution be denoted as \eqn{F_T(t)},
#' then the cumulative density function of the generalized gamma distribution can be
#' written as
#' \deqn{
#'    F(x) = F_T\bigg(\frac{x^b}{a}\bigg)
#' }
#' which allows us to write the quantile function of the generalized gamma in terms of
#' the gamma one (\eqn{Q_T(u)} is the quantile function of \eqn{T})
#' \deqn{
#'    Q(u) = (Q_T(u) \cdot a)^{1/b}
#' }
#' from which random numbers can be drawn.
#'
#' @references
#' Stacy, E. W. (1962). A generalization of the gamma distribution.
#' The Annals of mathematical statistics, 33(3), 1187-1192.
#'
#' @name G.Gamma
#' @aliases Generalized-Gamma
#' @aliases GGamma
#'
#' @keywords distribution
#' @keywords univar
#' @keywords models
#' @keywords survival
#' @concept Univariate
#' @concept Continuous
#' @concept Lifetime
#'
#' @export

dggamma = function(x, a, b, k, log=F){
	maxLength = max(length(x), length(a), length(b), length(k));
	if(length(x) != 1 && length(x) != maxLength) x = rep_len(x, maxLength)
	if(length(a) != 1 && length(a) != maxLength) a = rep_len(a, maxLength)
	if(length(b) != 1 && length(b) != maxLength) b = rep_len(b, maxLength)
	if(length(k) != 1 && length(k) != maxLength) k = rep_len(k, maxLength)

	result = log(b) - lgamma(k) + (b*k - 1)*log(x) - (b*k)*log(a) - (x/a)**b;
	if(!log) result = exp(result);
	return(result);
}

#' @rdname G.Gamma
#' @export

pggamma = function(q, a, b, k, lower.tail = TRUE, log.p = FALSE){
	cumProb = pgamma(q**b / a, shape=k, rate=1); # shape is alpha, rate is beta
	if(!lower.tail) cumProb = 1 - cumProb;
	if(log.p) cumProb = log(cumProb);
	return(cumProb);
}

#' @rdname G.Gamma
#' @export

qggamma = function(p, a, b, k, lower.tail = TRUE, log.p = FALSE){
	if(log.p) p = exp(p);
	if(!lower.tail) p = 1 - p;
	quantile = qgamma(p, shape=k, rate=1); # shape is alpha, rate is beta
	quantile = a * quantile**(1/b);
	return(quantile);
}


#' @rdname G.Gamma
#' @export

rggamma = function(n, a, b, k){
	return( qggamma(runif(0, 1, n=n), a, b, k) );
}
