#' @title Continuous Bernoulli Distribution
#' @description Density, distribution function, quantile function, and random generation for a one parameter continuous Bernoulli distribution
#' @param x,q a vector of quantities defined for \eqn{0 < X < 1}
#' @param n number of observations to generate
#' @param lambda the shape parameter \eqn{0 < \lambda < 1}
#' @param log,log.p a logical value indicating whether to log transform the probabilities. Default is \code{FALSE}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P(X \leq x)} otherwise, \eqn{P(X > x)}
#' @param p a vector of probabilities
#' @return a vector of numeric values
#' @details
#' The continuous Bernoulli distribution is a one parameter distribution (\eqn{\lambda \in (0,1)}) that has support \eqn{x \in (0,1)}.
#'
#' The density function is defined as:
#' \deqn{P(x) = C(\lambda) \lambda^x (1 - \lambda)^{1-x},}
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
#'
#' @references
#' Loaiza-Ganem, G., and J. P. Cunningham. 2019. The continuous Bernoulli: fixing a pervasive error in variational autoencoders. Page Advances in Neural Information Processing Systems. Curran Associates, Inc.
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

  p <- ifelse(x >= 1 | x <= 0,
              0,
              C * lambda^x*(1-lambda)^(1-x))

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

  p <- ifelse(
    q <= 0,
    0,
    ifelse(
      q >= 1,
      1,
      p
    )
  )

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

  q <- ifelse(p >= 1 | p <= 0,
              NaN,
              q)

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


#' @title Truncated Pareto Distribution
#' @description Density function and random number generation of a three parameter truncated Pareto distribution
#' @param x a vector of quantities
#' @param n integer value indicating the number of observations to generate
#' @param a the exponent of the Pareto distribution
#' @param b the minimum value of \code{x}
#' @param endpoint a numeric value where the Pareto distribution is truncated. Default value is 1.
#' @param log a logical value indicating whether to log transform the probabilities. Default is \code{FALSE}.
#' @details
#' The density function is defined as:
#'
#' If \eqn{x < B}
#' \deqn{P(x) = \frac{a b^a}{x^{a+1}}\frac{1}{F(x = B; a, b)},}
#' otherwise,
#' \deqn{P(x) = 0,}
#' where \eqn{B} is the end point where the Pareto distribution is truncated, \eqn{b} is the minimum value of x, \eqn{a} is the exponent of the Pareto distribution, and \eqn{F(x = B; a, b)} is the CDF of the Pareto distribution defined as
#' \deqn{F(x) = 1 - (\frac{b}{x})^a}
#'
#' @return a vector of numeric values
#' @rdname tpareto
#' @export
dtpareto<-function(x, a, b, endpoint = 1, log = FALSE){
  prob<-ifelse(
    x <= endpoint,
    extraDistr::dpareto(x, a = a, b = b) /
      extraDistr::ppareto(endpoint, a = a, b = b),
    0)
  if(log){
    prob <- log(prob)
  }
  return(prob)
}


#' @title Truncated Pareto Distribution
#' @rdname tpareto
#' @export
rtpareto<-function(n, a, b, endpoint = 1){
  x<-extraDistr::qpareto(
    p =
      runif(n = n, #Inverse sampling but shrink original support to location of truncation
            min = 0,
            max = extraDistr::ppareto(q = endpoint,
                                      a = a,
                                      b = b)),
    a = a,
    b = b
  )
  x[x>endpoint] <- NaN
  return(x)
}


#' @title Truncated Log-Normal Hurdle Distribution
#' @description Density function and random number generation of a three parameter truncated log-normal hurdle distribution. The log-normal distribution does not normally have support at zero and has an upper bound at positive infinity. Setting up the distribution as a hurdle model takes care of zeros in data. The truncation sets the upper bound to 1 or any arbitrary user defined value.
#' @param x A vector of numeric values
#' @param n integer value indicating the number of observations to generate
#' @param theta A numeric value between zero and one (inclusive) that defines the probability of non-zero values
#' @param meanlog,sdlog Mean and standard deviation of the distribution on the log scale.
#' @param endpoint A numeric value indicating the upper bound of the distribution
#' @param log a logical value indicating whether to log transform the probabilities. Default is \code{FALSE}.
#' @details
#' The density function is defined as:
#'
#' If \eqn{x = 0}
#' \deqn{f(x) = \theta}
#' If \eqn{x > 0}
#' \deqn{f(x) = (1 - \theta) \frac{g(x, \mu, \sigma)}{G(x = b; \mu ,\sigma)} }
#' where the log-normal PDF is defined as
#' \deqn{ g(x; \mu , \sigma) = \frac{1}{x \sigma \sqrt{2 \pi}} \exp{(-\frac{(\ln{x}-\mu)^2}{2\sigma^2})} }
#' and the log-normal CDF is \eqn{G(x = b; \mu ,\sigma)}.\eqn{b} is the upper bound (endpoint), \eqn{\mu} is the mean, \eqn{\sigma} is the standard deviation, and \eqn{\theta} is the proportion of zeros.
#' @rdname htlnorm
#' @return A vector of probability densities
dhtlnorm<-function(x,theta,meanlog,sdlog,endpoint =1,log=FALSE){
  p <- ifelse(
    x == 0,
    theta,
    (1-theta) * ifelse(
      x <= endpoint,
      dlnorm(x, meanlog = meanlog, sdlog = sdlog) /
        plnorm(endpoint,meanlog = meanlog,sdlog = sdlog),
      0)
  )

  if(log){
    p <- log(p)
  }
  return(p)
}

#' @title Truncated Log-Normal Hurdle Distribution
#' @rdname htlnorm
#' @export
rhtlnorm<-function(n,theta,meanlog,sdlog,endpoint=1){
  out <- ifelse(
    rbinom(n, size = 1, prob = theta) == 1,
    0,
    qlnorm(
      p =
        runif(n = n, #Inverse sampling but shrink original support to location of truncation
              min = 0,
              max = plnorm(q = endpoint,
                           meanlog = meanlog,
                           sdlog = sdlog)),
      meanlog = meanlog,
      sdlog = sdlog
    )
  )
  return(out)
}


#' @title Zero-One-Inflated Beta Distribution
#' @description Density function and random number generation of a four parameter zero-one-inflated beta distribution. The beta distribution does not have defined probabilities at 0 or 1 so this extension takes care of that.
#' @param x A vector of numeric values
#' @param n number of observations to generate
#' @param p0 The probability that \eqn{x = 0}. Must be between 0 and 1
#' @param p1 The conditional probability that \eqn{x = 1} given \eqn{x \neq 0}. Must be between 0 and 1.
#' @param alpha,beta The shape parameters of the beta distribution. Must be greater than 0.
#' @param log a logical value indicating whether to log transform the probabilities. Default is \code{FALSE}.
#' @details
#' The density function is defined as:
#'
#' If \eqn{x = 0}
#' \deqn{f(x) = p_0}
#' If \eqn{x = 1}
#' \deqn{f(x) = p_1 (1 - p_0)}
#' If \eqn{x \neq 0,1}
#' \deqn{f(x) = (1 - p_1) (1 - p_0) \frac{x^{\alpha-1}(1-x)^{\beta-1}}{B(\alpha, \beta)}}
#' @rdname zoibeta
#' @return A vector of probability densities
#' @references
#' Ospina, R., and S. L. P. Ferrari. 2008. Inflated beta distributions. Statistical Papers 51:111.
dzoibeta<-function(x,p0,p1,alpha,beta,log=FALSE){
  p <- ifelse(
    x == 0,
    p0,
    ifelse(
      x == 1,
      p1 * (1 - p0),
      (1 - p1) * (1 - p0) * dbeta(x, shape1 = alpha, shape2 = beta)
    )
  )
  if(log){
    p <- log(p)
  }
  return(p)
}

#' @title Zero-One-Inflated Beta Distribution
#' @rdname zoibeta
#' @export
rzoibeta<-function(n,p0,p1,alpha,beta){
  out<-ifelse(
    rbinom(n = n,size = 1,prob = p0) == 1,
    0,
    ifelse(
      rbinom(n = n, size = 1, prob = p1) == 1,
      1,
      rbeta(n,alpha,beta)))
  return(out)
}







