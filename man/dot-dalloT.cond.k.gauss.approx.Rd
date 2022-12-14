% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neutral_model_ML_support.R
\name{.dalloT.cond.k.gauss.approx}
\alias{.dalloT.cond.k.gauss.approx}
\title{Calculate the kth Convolution Probability Matrix of P(phi) Via Gaussian Approximation}
\usage{
.dalloT.cond.k.gauss.approx(
  phi.T,
  k,
  min.phi = 0.005,
  max.phi = 1,
  a = 14/9,
  log = FALSE
)
}
\arguments{
\item{phi.T}{a vector of cumulative proportion herbivory}

\item{k}{a vector of integer values indicating the number of convolutions to perform}

\item{min.phi}{the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5\%).}

\item{max.phi}{the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100\%).}

\item{a}{the combined allometric scaling coefficient. Defaults to 14/9.}

\item{log}{if \code{TRUE} (default is \code{FASLE}), return probabilities on the log scale.}
}
\value{
a length(phi.T) X length(k) matrix of probabilities.
}
\description{
Internal function of \code{.dalloT.cond.k.fft.conv()} used to calculate \eqn{P(\phi^k)} via Gaussian approximation for very large values of \eqn{k}.
}
\details{
Owing to CLT,
\deqn{P(\phi_T^k) \approx \mathcal{N}(\mu = k\overline\phi,\sigma=\sqrt{k\mathbb{Var}[\phi]})}
Also see \code{?get_phi_bar()} and \code{?get_phi_var()}
}
