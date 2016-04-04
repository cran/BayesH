\name{predict.BayesH}
\alias{predict.BayesH}
\title{Make Predictions}
\description{
    Make predictions from the results
    of bayesModel.fit function.
}
\usage{
\method{predict}{BayesH}(object, newdata, ...)

}
\arguments{
  \item{object}{An object of class \code{BayesH}.}
  \item{newdata}{An object of class matrix or data.frame with new data (explanatory variables).}
  \item{...}{Further arguments passed to or from other methods.}
}
\author{
   Renato Rodrigues Silva
}
\seealso{
  \code{bayesfit.model}.
}

\references{ Gustavo de los Campos, Hugo Naya, Daniel Gianola, Jose Crossa, Andre's Legarra,
             Eduardo Manfredi, Kent Weigel and Jose Miguel Cotes. Predicting Quantitative Traits 
             With Regression Models for Dense Molecular Markers and Pedigree, 2009. Genetics. 
             doi: 10.1534/genetics.109.101501
             
           }

