% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weirdo_distributions.R
\name{dhtlnorm}
\alias{dhtlnorm}
\alias{rhtlnorm}
\title{Truncated Log-Normal Hurdle Distribution}
\usage{
dhtlnorm(x, theta, meanlog, sdlog, endpoint = 1, log = FALSE)

rhtlnorm(n, theta, meanlog, sdlog, endpoint = 1)
}
\arguments{
\item{x}{A vector of numeric values}

\item{theta}{A numeric value between zero and one (inclusive) that defines the probability of non-zero values}

\item{meanlog, sdlog}{Mean and standard deviation of the distribution on the log scale.}

\item{endpoint}{A numeric value indicating the upper bound of the distribution}

\item{log}{a logical value indicating whether to log transform the probabilities. Default is \code{FALSE}.}

\item{n}{integer value indicating the number of observations to generat}
}
\value{
A vector of probability densities
}
\description{
Density function and random number generation of a three parameter truncated log-normal hurdle distribution. The log-normal distribution does not normally have support at zero and has an upper bound at positive infinity. Setting up the distribution as a hurdle model takes care of zeros in data. The truncation sets the upper bound to 1 or any arbitrary user defined value. The truncated probabilities are added to the end point as a Dirac delta function.
}
\details{
The density function is defined as:

If \eqn{x = 0}
\deqn{f(x) = \theta}
if \eqn{x = 1}
\deqn{f(x) = (1 - \theta) (1-G(x = b; \mu ,\sigma))}
if \eqn{1 > x > 0}
\deqn{f(x) = (1 - \theta) g(x, \mu, \sigma)}
where the log-normal PDF is defined as
\deqn{ g(x; \mu , \sigma) = \frac{1}{x \sigma \sqrt{2 \pi}} \exp{(-\frac{(\ln{x}-\mu)^2}{2\sigma^2})} }
and the log-normal CDF is \eqn{G(x = b; \mu ,\sigma)}.\eqn{b} is the upper bound (endpoint), \eqn{\mu} is the mean, \eqn{\sigma} is the standard deviation, and \eqn{\theta} is the proportion of zeros.
}
