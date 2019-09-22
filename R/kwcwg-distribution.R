

#' KW-CWG distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the KW-CWG distribution.
#'
#' @param x,q	          vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mu,sigma	      location and scale parameters. Scale must be positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' @TODO: Place the correct formulas below
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{1}{2\sigma} \exp\left(-\left|\frac{x-\mu}{\sigma}\right|\right)
#' }{
#' f(x) = 1/(2*\sigma) * exp(-|(x-\mu)/\sigma|)
#' }
#'
#' Cumulative distribution function
#' \deqn{ F(x) = \left\{\begin{array}{ll}
#' \frac{1}{2} \exp\left(\frac{x-\mu}{\sigma}\right)     & x < \mu \\
#' 1 - \frac{1}{2} \exp\left(\frac{x-\mu}{\sigma}\right) & x \geq \mu
#' \end{array}\right.
#' }{
#' F(x) = [if x < mu:] 1/2 * exp((x-\mu)/\sigma)
#'        [else:] 1 - 1/2 * exp((x-\mu)/\sigma)
#' }
#'
#' Quantile function
#' \deqn{ F^{-1}(p) = \left\{\begin{array}{ll}
#' \mu + \sigma \log(2p)     & p < 0.5 \\
#' \mu - \sigma \log(2(1-p)) & p \geq 0.5
#' \end{array}\right.
#' }{
#' F^-1(p) = [if p < 0.5:] \mu + \sigma * log(2*p)
#'           [else:] \mu - \sigma * log(2*(1-p))
#' }
#'
#' @references
#' Afify, A.Z., Cordeiro, G.M., Butt, N.S., Ortega, E.M. and
#' Suzuki, A.K. (2017). A new lifetime model with variable shapes for
#' the hazard rate. Brazilian Journal of Probability and Statistics
#'
#' @references
#' Forbes, C., Evans, M. Hastings, N., & Peacock, B. (2011).
#' Statistical Distributions. John Wiley & Sons.
#' 
#' @examples 
#' 
#' x <- rkwcwg(1e5, 5, 16)
#' hist(x, 100, freq = FALSE)
#' curve(dkwcwg(x, 5, 16), -200, 200, n = 500, col = "red", add = TRUE)
#' hist(pkwcwg(x, 5, 16))
#' plot(ecdf(x))
#' curve(pkwcwg(x, 5, 16), -200, 200, n = 500, col = "red", lwd = 2, add = TRUE)
#'
#' @name KW-CWG
#' @aliases KW-CWG
#' @aliases dkwcwg
#'
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#' @concept Lifetime
#'
#' @export

dkwcwg <- function(x, mu = 0, sigma = 1, log = FALSE) {
  cpp_dkwcwg(x, mu, sigma, log[1L])
}


#' @rdname KW-CWG
#' @export

pkwcwg <- function(q, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pkwcwg(q, mu, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname KW-CWG
#' @export

qkwcwg <- function(p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qkwcwg(p, mu, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname KW-CWG
#' @export

rkwcwg <- function(n, mu = 0, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rkwcwg(n, mu, sigma)
}

