% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{r2_partial.brmsfit}
\alias{r2_partial.brmsfit}
\title{Partial R2 of Fitted 'brmsfit' Model Predictor}
\usage{
\method{r2_partial}{brmsfit}(
  object,
  var,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...
)
}
\arguments{
\item{object}{A 'brmsfit' object}

\item{var}{the name of the variable for which to assess the partial R2}

\item{summary}{if \code{TRUE} (default), the posterior draws are summarized}

\item{robust}{if \code{TRUE} (default is FALSE), the median instead of the mean is returned as the estimate}

\item{probs}{The lower and upper interval of posterior summary}

\item{...}{additional arguments}
}
\value{
a matrix of partial r2 draws or a vector of partial r2 summary.
}
\description{
Calculate the Bayesian partial R2 of a specified predictor variable by finding the change in r2 when the predictor variable is set to zero.
}
