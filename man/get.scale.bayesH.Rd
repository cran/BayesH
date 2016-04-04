\name{get.scale.bayesH}
\alias{get.scale.bayesH}

\title{Get the vector of the scale hyperparameter}

\description{Get the vector of the scale hyperparameter.}

\usage{ 
   get.scale.bayesH(nu0, R2 = 0.5, X, y, type = "bayesH")   
}

\arguments{  
  \item{nu0}{the vector with shape parameter(s) of (hyper) prior distribution for 
             variance of each the regression coefficient(s) and error variance.} 
  \item{R2}{the expected  proportion  of the variance of the response variable that 
            can be explained by the model (r squared).}
  \item{y}{the vector of response variable of the model.}
  \item{X}{the incidence matrix.} 
  \item{type}{it is a string which if were defined as ``ridge'' the function performs
              Bayesian ridge regression, otherwise, Bayes H model.}    
}

\value{ An object of class \code{vector} with values of the scale parameters.
}

\details{ The user must be provided the  vector \code{nu0} and all elements this vector 
          must be greater than zero. The strategy adopted here is  based on Perez et. al. 
          (2010).
}         

\seealso{ \code{\link{bayesModel.fit}} } 

\references{ Perez. et. al. (2010). Genomic-Enabled Prediction Based on Molecular
             Markers and Pedigree Using the Bayesian Linear Regression Package 
             in R. Plant Genome. 2010; 3(2): 106-116. doi:  10.3835/plantgenome2010.04.0005.          
           }

\author{Renato Rodrigues Silva, \email{renato.rrsilva@ufg.br} }



