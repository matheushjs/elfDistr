

#' KW-CWG distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the KW-CWG distribution.
#'
#' @param x,q	          vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param alpha,beta,gamma,a,b	      Parameters of the distribution. 0 < alpha < 1, and the other parameters mustb e positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#'    f(x) = \alpha^a \beta \gamma a b (\gamma x)^{\beta - 1} \exp[-(\gamma x)^\beta] \cdot
#'    \frac{\{1 - \exp[-(\gamma x)^\beta]\}^{a-1}}{\{ \alpha + (1 - \alpha) \exp[-(\gamma x)^\beta] \}^{a+1}} \cdot \\
#'    \cdot \bigg\{ 1 - \frac{\alpha^a[1 - \exp[-(\gamma x)^\beta]]^a}{\{ \alpha + (1 - \alpha) \exp[-(\gamma x)^\beta] \}^a} \bigg\}
#' }
#'
#' @references
#' Afify, A.Z., Cordeiro, G.M., Butt, N.S., Ortega, E.M. and
#' Suzuki, A.K. (2017). A new lifetime model with variable shapes for
#' the hazard rate. Brazilian Journal of Probability and Statistics
#' 
#' @name KW-CWG
#' @aliases KW-CWG
#' @aliases kwcwg
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

dkwcwg <- function(x, alpha, beta, gamma, a, b, log = FALSE) {
  cpp_dkwcwg(x, alpha, beta, gamma, a, b, log[1L])
}


#' @rdname KW-CWG
#' @export

pkwcwg <- function(q, alpha, beta, gamma, a, b, lower.tail = TRUE, log.p = FALSE) {
  cpp_pkwcwg(q, alpha, beta, gamma, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname KW-CWG
#' @export

qkwcwg <- function(p, alpha, beta, gamma, a, b, lower.tail = TRUE, log.p = FALSE) {
  cpp_qkwcwg(p, alpha, beta, gamma, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname KW-CWG
#' @export

rkwcwg <- function(n, alpha, beta, gamma, a, b) {
  if (length(n) > 1) n <- length(n)
  cpp_rkwcwg(n, alpha, beta, gamma, a, b)
}

