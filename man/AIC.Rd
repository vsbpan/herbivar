% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neutral_model_ML_support.R,
%   R/null_model_ML_support.R
\name{AIC.allo_herb_fit}
\alias{AIC.allo_herb_fit}
\alias{BIC.allo_herb_fit}
\alias{AIC.generic_null_fit}
\alias{BIC.generic_null_fit}
\title{Akaike's An Information Criterion}
\usage{
\method{AIC}{allo_herb_fit}(object, ..., k = 2)

\method{BIC}{allo_herb_fit}(object, ...)

\method{AIC}{generic_null_fit}(object, ..., k = 2)

\method{BIC}{generic_null_fit}(object, ...)
}
\arguments{
\item{object}{a fitted model object for which there exist a \code{logLik} method to extract the log-likelihood, or an object of class 'logLik'.}

\item{...}{additional arguments passed to \code{stats::AIC()} or \code{stats::BIC()} for additional fitted model objects.}

\item{k}{an atomic numeric value for the penalty. Defaults to 2.}
}
\value{
if just one object is provided, a numeric value of AIC or BIC is returned. If multiple objects are provided, a data.frame with rows corresponding to the objects and columns representing the number of parameters in the model and the AIC or BIC is returned.
}
\description{
calculate AIC or BIC for fitted model objects. see \code{stats::AIC()} for more details.
}
