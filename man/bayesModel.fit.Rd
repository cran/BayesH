\name{bayesModel.fit}
\alias{bayesModel.fit}

\title{Bayesian regression model with mixture of two scaled inverse
       chi square as hyperprior distribution. }

\description{ Performs Gibbs Sampling algorithm for fitting the
              Bayesian regression model with mixture of two scaled inverse
              chi square as hyperprior distribution for variance of each regression 
              coefficients.
}

\usage{ 
   bayesModel.fit(X, y,  nu0, s0, niter = 2000, burnin = 500, type="bayesH")   
}

\arguments{
  \item{X}{the incidence \code{matrix} of model.  }
  \item{y}{the vector of response variable of the model.}
  \item{nu0}{the degree of freedom hyperparameter(s) \code{nu0} for all mixture components.} 
  \item{s0}{the scale  hyperparameter(s) \code{s0} for all mixture components.}  
  \item{niter}{the number of iterations of Gibbs Sampling algorithm.} 
  \item{burnin}{the number of 'burn in' in a Gibbs Sampling algorithm. } 
  \item{type}{it is a string which if were defined as ``ridge'' the function performs
              Bayesian ridge regression, otherwise, Bayes H model. } 
}

\value{ The output is an object of class \code{BayesH} that contains the
        posterior distribution of intercept, posterior distribution of variance error,        
        posterior mean of regression coefficients and posterior mean of predicted values.
}

\details{ For bayesian ridge regression (type == "ridge"), the prior distribution for the error
          variance and the hyperprior distribution for variance of the regression coefficients 
          follows scaled inverse chi square with same hyperparameters \code{(nu0[1], s0[1])} 
          and \code{(nu0[2], s0[2])}, respectively.On the other hand, for hierarchical regression 
          model (type == "bayesH") is assumed that each the regression coefficient has different variance
          and each one of them follows a  mixture of scaled inverse chi square with  hyperparameters 
          (\code{nu0[1]}; \code{s0[1]}) and (\code{nu0[2]}; \code{s0[2]}), respectively. 
          In this case,  the prior distribution  for error variance also follows scaled inverse 
          chi square  with hyperparameters \code{nu0[3]} and \code{s0[3]}. NA's in the incidence 
          matrix are not allowed. All elements of vector \code{s0} must be greater than zero.
}

\seealso{ \code{\link{get.scale.bayesH}} } 


\author{Renato Rodrigues Silva, \email{renato.rrsilva@ufg.br} }

\examples{

data(example)

mod = bayesModel.fit(X = X, y = y, nu0 = c(3,30,20), 
                     s0 =  c(1e-10,0.1, 0.5), niter = 2000, 
                     burnin = 300, type = "bayesH" )
summary(mod)


}
