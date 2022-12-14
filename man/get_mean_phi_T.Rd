% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neutral_model_ML_support.R
\name{get_mean_phi_T}
\alias{get_mean_phi_T}
\title{Calculate Unlimited Mean Cumulative Proportion Herbivory From Attack Rate}
\usage{
get_mean_phi_T(lambda, min.phi = 0.005, max.phi = 1, a = 14/9)
}
\arguments{
\item{lambda}{the attack rate on the plant or leaf}

\item{min.phi}{the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5\%).}

\item{max.phi}{the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100\%).}

\item{a}{the combined allometric scaling coefficient. Defaults to 14/9.}
}
\value{
a vector of numeric values
}
\description{
Convenience function that reparameterizes attack rate \eqn{\lambda} as unlimited mean cumulative proportion herbivory \eqn{\overline{\phi_T'}}.
}
\details{
\deqn{\overline{\phi_T'} = \overline{\phi}\lambda}
}
