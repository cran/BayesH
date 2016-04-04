\name{summary.BayesH}
\alias{summary.BayesH}
\title{Summary}
\description{
     Provides a summary for a fitted model using bayesModel.fit function.
}
\usage{
\method{summary}{BayesH}(object, ...)

}
\arguments{
  \item{object}{An object of class \code{BayesH}.}  
  \item{...}{Further arguments passed to or from other methods.}
}
\author{
   Renato Rodrigues Silva
}


\seealso{
  \code{bayesModel.fit}.
}

\details{ 
    It is an adaptation of code written by Perez and De Los Campos (2014). The output of this function is 
    print the mean of the posterior distribution of the error variance, the variance of reponse observed values 
    (training population). When there are missing values, the correlation between the predicted values and observed values
    is printed as well.   

}

\references{ Gustavo de los Campos, Hugo Naya, Daniel Gianola, Jose Crossa, Andre's Legarra,
             Eduardo Manfredi, Kent Weigel and Jose Miguel Cotes. Predicting Quantitative Traits 
             With Regression Models for Dense Molecular Markers and Pedigree, 2009. Genetics. 
             doi: 10.1534/genetics.109.101501
             
           }
