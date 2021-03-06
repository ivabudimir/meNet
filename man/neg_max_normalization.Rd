% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization_functions.R
\name{neg_max_normalization}
\alias{neg_max_normalization}
\title{Negative maximum normalization}
\usage{
neg_max_normalization(x = numeric(0))
}
\arguments{
\item{x}{Non-negative numeric vector.}
}
\value{
Normalized numeric vector.
}
\description{
Performs the negative maximum normalization on a non-negative 
numeric vector. The resulting vector is both reversed and normalized.
}
\details{
For negative maximum normalization of a vector, difference between the 
maximum and minimal value of the vector is added tothe inverse value of every 
element and the result is divided with the original maximum value. 
Used formula: `(max(x)-min(x)-x)/max(x)`.

If the negative maximum normalization is applied twice to a 
non-negative vector, obtained result is the same as if the maximum 
normalization was applied once.
}
\examples{
{
neg_max_normalization(c(1,2,3,4))
neg_max_normalization(neg_max_normalization(c(1,2,3,4)))
}

}
