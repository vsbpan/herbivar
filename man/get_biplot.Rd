% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{get_biplot}
\alias{get_biplot}
\title{Get Biplot or Triplot}
\usage{
get_biplot(
  x,
  choices = c(1, 2),
  scaling = 2,
  display = c("sites", "species", "biplot", "centroids"),
  group = NULL
)
}
\arguments{
\item{x}{an object that is supported by \code{vegan::scores()}}

\item{choices}{a vector of length 2 defining the axes to plot}

\item{scaling}{scaling argument passed to \code{vegan::scores()}}

\item{display}{a vector of characters defining what to plot}

\item{group}{a vector in the same length and order as the sites data used to fit the model that is used to color the site points}
}
\value{
a ggplot object
}
\description{
Generate biplot or triplot using ggplot2.
}
