#' Generalized Gamma Probability Distribution
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
