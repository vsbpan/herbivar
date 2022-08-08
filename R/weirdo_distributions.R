#' @title Adjust Proportion Data
#' @description Adjusts proportion data to get rid of 0's and 1's, and optionally transforms the resulting adjusted proportion.
#' @param x a vector of numeric data bounded between 0 and 1 inclusive.
#' @param trans a character value indicating the transformation to be applied to the nudged data. Valid options are "identity" (default-- no transformation), "logit", "empirical_logit", "probit", and "log". See details.
#' @param nudge.method a character value indicating the method of getting rid of 0's and 1's. Valid options are "replace" (default), "smithson", "add", "subtract", and "drop". see details.
#' @param nudge.size a numeric value of the size of a small nudge (\eqn{epsilon}). Valid character values can also be supplied to choose the method for estimating \eqn{epsilon}. Valid methods are "macmillan" (default), "warton_min", and "warton_max". see details.
#' @param bounds a character value indicating bounds for which the \code{nudge.method} is applied to. Valid options are "both" (default), "zero", and "one". Only relevant for \code{nudge.method} set as "replace" or "drop".
#' @param na.action a character value indicating how to handle missing values. Valid options are "ignore", "remove", "as.is", and "fail". See details.
#'
#' @return a vector of numeric value with the same length as \code{x}
#' @details
#'
#' @references
#' Kubinec, R. 2022. Ordered Beta Regression: A Parsimonious, Well-Fitting Model for Continuous Data with Lower and Upper Bounds. Political Analysis:1–18.
#'
#' Macmillan, N. A., and C. D. Creelman. 2004. Detection Theory: A User’s Guide. Second edition. Psychology Press, New York.
#'
#' Smithson, M., and J. Verkuilen. 2006. A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. Psychological Methods 11:54–71.
#'
#' Warton, D. I., and F. K. C. Hui. 2011. The arcsine is asinine: the analysis of proportions in ecology. Ecology 92:3–10.
#' @export
#'
#' @examples
#' x <- c(0,runif(0,0,1),NA,1,1,1)
#' adjust_prop(x, nudge.method = "replace", nudge.size = 0.01)
#'
#' adjust_prop(x, nudge.method = "replace", nudge.size = "macmillan", trans = "logit")
#'
#'
adjust_prop <- function(x,
                     nudge.method = c("replace","smithson", "add", "subtract", "drop"),
                     nudge.size = c("macmillan", "warton_min", "warton_max"),
                     trans = c("identity","logit","empirical_logit","probit","log"),
                     bounds = c("both","zero","one"),
                     na.action = c("ignore","remove","as.is","fail")){

  nudge.method <- match.arg(nudge.method)
  bounds <- match.arg(bounds)
  trans <- match.arg(trans)
  na.action <- match.arg(na.action)

  bounds.numeric <- switch(bounds,
                           "both" = c(0,1),
                           "zero" = 0,
                           "one" = 1)

  if(is.character(nudge.size)){
    nudge.size <- match.arg(nudge.size)
  } else {
    nudge.size <- as.numeric(nudge.size)[1]
  }

  na.v <- is.na(x)

  if(any(x[!na.v] > 1)){
    stop("'x' contains values greater than 1")
  }
  if(any(x[!na.v] < 0)){
    stop("'x' contains negative values")
  }

  if(nudge.method == "drop"){
    x[x %in% bounds.numeric] <- NA
    if(na.action != "remove"){
      if(na.action != "ignore"){
        na.action <- "ignore"
        message("na.action defaulting to 'ignore'")
      }
    }
  }

  if(any(na.v)){
    if(na.action == "fail"){
      stop("'x' contains missing values")
    }

    if(na.action == "remove"){
      x <- x[!na.v]
    } else {
      if(na.action == "ignore"){
        x_na_rm <- x[!na.v]
        warning("'x' contains ignored missing values.")
      } else {
        warning("'x' contains missing values. Values computed as is.")
      }
      }
  } else {
    x_na_rm <- x
  }

  if(is.numeric(nudge.size)){
    epsilon <- nudge.size
    nudge.size <- "user_supplied"
  } else {
    if(nudge.method == "smithson"){
      message("epsilon defaulting to 0.5")
      nudge.size <- "smithson_default"
      epsilon <- 0.5
    } else {
      if(na.action == "ignore"){
        epsilon <- switch(nudge.size,
                          "macmillan" = 1/(2*length(x_na_rm)),
                          "warton_min" = min(x_na_rm[x_na_rm > 0]),
                          "warton_max" = min((1-x_na_rm)[(1-x_na_rm) > 0])
        )
      } else {
        epsilon <- switch(nudge.size,
                          "macmillan" = 1/(2*length(x)),
                          "warton_min" = min(x[x > 0]),
                          "warton_max" = min((1-x)[(1-x) > 0])
        )
      }
    }
  }

  if(nudge.method == "replace"){
    if(0 %in% bounds.numeric){
      x <- ifelse(x == 0,
                  epsilon,
                  x)
    }
    if(1 %in% bounds.numeric){
      x <- ifelse(x == 1,
                  1-epsilon,
                  x)
    }
  } else {
    if(nudge.method == "add"){
      x <- x + epsilon
    } else {
      if(nudge.method == "subtract"){
        x <- x - epsilon
      } else {
        if(nudge.method == "smithson"){
          if(na.action == "ignore"){
            n <- length(x_na_rm)
          } else {
            n <- length(x)
          }
          x <- (x * (n-1) + epsilon) / n
        }
      }
    }
  }

  if(any(x[!is.na(x)] > 1)){
    warning("nudged 'x' contains values greater than 1")
  }
  if(any(x[!is.na(x)] < 0)){
    warning("nudged 'x' contains negative values")
  }

  x.trans <- switch(trans,
                    "identity" = x,
                    "logit" = qlogis(x),
                    "empirical_logit" = log((x + epsilon) / (1 - x + epsilon)),
                    "probit" = qnorm(x),
                    "log" = log(x))
  x.trans <- structure(x.trans,
                       "adj.method" = list(
                         "epsilon" = epsilon,
                         "nudge.size" = nudge.size,
                         "nudge.method" = nudge.method,
                         "trans" = trans,
                         "bounds" = bounds,
                         "na.action" = na.action
                       )
                       )
  return(x.trans)
}




#' @title Continuous Bernoulli Distribution
#' @description Density, distribution function, quantile function, and random generation for a one parameter continuous Bernoulli distribution
#' @param x,q a vector of quantities defined for \eqn{0 \leq X \leq 1}
#' @param n number of observations to generate
#' @param lambda the shape parameter \eqn{0 < \lambda < 1}
#' @param log,log.p a logical value indicating whether to log transform the probabilities. Default is \code{FALSE}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P(X \leq x)} otherwise, \eqn{P(X > x)}
#' @param p a vector of probabilities
#' @return a vector of numeric values
#' @details
#' The continuous Bernoulli distribution is a one parameter distribution (\eqn{\lambda \in (0,1)}) that has support \eqn{x \in [0,1]}.
#'
#' The density function is defined as:
#' \deqn{P(X = x) = C(\lambda) \lambda^x (1 - \lambda)^{1-x},}
#'
#' where \eqn{C(\lambda)} is defined as
#' \deqn{C(\lambda) = 2}
#' if \eqn{\lambda = \frac{1}{2}}, otherwise,
#' \deqn{C(\lambda) = \frac{2 \tanh^{-1}(1-2\lambda)}{1-2\lambda}.}
#'
#' The cumulative density function is defined as
#' \deqn{P(X \leq x) = x}
#' if \eqn{\lambda = \frac{1}{2}}, otherwise,
#' \deqn{P(X \leq x) = \frac{\lambda^x (1-\lambda)^{1-x} + \lambda - 1}{2 \lambda - 1}.}
#'
#' @rdname cb
#' @export
dcb <- function(x, lambda, log = FALSE){
  if(lambda >=1 || lambda <= 0){
    stop("lambda must be between 0 and 1, not inclusive.")
  }
  if(lambda == 0.5){
    C <- 2
  } else {
    C <- 2*atanh(1-2*lambda)/(1-2*lambda)
  }

  # x has support on the boundaries!!

  p <- C * lambda^x*(1-lambda)^(1-x)
  if(log){
    p <- log(p)
  }
  return(p)
}

#' @rdname cb
#' @export
pcb <- function(q, lambda, lower.tail = TRUE, log.p = FALSE){
  if(lambda >=1 || lambda <= 0){
    stop("lambda must be between 0 and 1, not inclusive.")
  }
  if(lambda == 0.5){
    p <- q
  } else {
    p <- (lambda^q * (1-lambda)^(1-q) + lambda - 1) / (2*lambda - 1)
  }
  if(!lower.tail){
    p <- 1-p
  }

  if(log.p){
    p <- log(p)
  }
  return(p)
}

#' @rdname cb
#' @export
qcb <- function(p, lambda, lower.tail = TRUE, log.p = FALSE){
  if(lambda >=1 || lambda <= 0){
    stop("lambda must be between 0 and 1, not inclusive.")
  }
  if(log.p){
    p <- exp(p)
  }
  if(!lower.tail){
    p <- 1 - p
  }
  if(lambda == 0.5){
    q <- p
  } else {
    q <- log(p*(2*lambda - 1) + 1, base = lambda)
  }
  return(q)
}

#' @rdname cb
#' @export
rcb <- function(n, lambda){
  if(lambda >=1 || lambda <= 0){
    stop("lambda must be between 0 and 1, not inclusive.")
  }
  x <- qcb(p = runif(n,0,1), lambda = lambda, lower.tail = TRUE, log.p = FALSE)
  return(x)
}

#' @title Ordered Beta Distribution ??
#'






