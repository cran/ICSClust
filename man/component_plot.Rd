% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{component_plot}
\alias{component_plot}
\title{Scatterplot Matrix with densities on the diagonal}
\usage{
component_plot(
  object,
  select = TRUE,
  clusters = NULL,
  text_size_factor = 8/6.5,
  colors = NULL
)
}
\arguments{
\item{object}{a dataframe or \code{\link[ICS]{ICS-S3}} class object.}

\item{select}{a vector of indexes of variables to plot. If \code{NULL} or
\code{FALSE}, all variables are selected. If \code{TRUE} only the first three and
last three are considered.}

\item{clusters}{a vector indicating the clusters of the data to color the
plot. By default \code{NULL}.}

\item{text_size_factor}{a numeric factor for controlling the \code{axis.text}
and \code{strip.text}.}

\item{colors}{a vector of colors to use. One color for each cluster.}
}
\value{
An object of class \code{\link[GGally:ggmatrix]{"ggmatrix"}} (see
\code{\link[GGally:ggpairs]{GGally::ggpairs()}}).
}
\description{
Produces a gg-scatterplot matrix of the variables of a given dataframe
or an invariant coordinate system obtained via an ICS transformation with
densities on the diagonal for each cluster.
}
\examples{
X <- iris[,1:4]
component_plot(X)
out <- ICS(X)
component_plot(out, select = c(1,4))

}
\author{
Andreas Alfons and Aurore Archimbaud
}
