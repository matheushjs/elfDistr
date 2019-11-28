
#' Kumaraswamy Complementary Weibull Geometric and Generalized Gamma Probability Distributions
#'
#' Density, distribution function, quantile function and random
#' generation for the Kumaraswamy Complementary Weibull Geometric
#' probability distribution (Kw-CWG) and the Generalized Gamma
#' lifetime distributions.
#' 
#' @details
#' 
#' This package follows naming convention that is consistent with base R,
#' where density (or probability mass) functions, distribution functions,
#' quantile functions and random generation functions names are followed by
#' \code{d}, \code{p}, \code{q}, and \code{r} prefixes.
#' 
#' Behaviour of the functions is consistent with base R, where for
#' not valid parameters values \code{NaN}'s are returned, while
#' for values beyond function support \code{0}'s are returned
#' (e.g. for non-integers in discrete distributions, or for
#' negative values in functions with non-negative support).
#' 
#' The most complex distributions were vectorized and coded in C++ using \pkg{Rcpp}.
#' Simpler distributions were coded directly in R, as this proved itself most efficient.
#' 
#' @docType package
#' @name elfDistr
#' 
#' @useDynLib elfDistr, .registration = TRUE 
#' @importFrom Rcpp sourceCpp
NULL
