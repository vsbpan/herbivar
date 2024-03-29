% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{r2_delta}
\alias{r2_delta}
\title{Delta R2 of Fitted Model Predictor}
\usage{
r2_delta(object, var, ...)
}
\arguments{
\item{object}{a supported model fit object}

\item{var}{the name of the variable for which to assess the delta R2}

\item{...}{additional parameters}
}
\value{
the delta r2
}
\description{
Generic method of finding the delta R2 of a specified predictor variable by finding the change in r2 when the predictor variable is set to zero.
}
