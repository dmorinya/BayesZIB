\encoding{utf8} 
\name{bayesZIB-package}
\alias{bayesZIB-package}
\docType{package}
\title{
Bayesian zero-inflated Bernoulli regression model
}
\description{
This package fits a Bayesian Bernoulli zero-inflated regression model handling different covariates for the zero-inflated and non zero-inflated parts.
}
\details{
\tabular{ll}{
Package: \tab bayesZIB\cr
Type: \tab Package\cr
Version: \tab 0.0.2\cr
Date: \tab 2021-5-12\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
The package implements a new Bayesian Bernoulli zero-inflated. This model is able to distinguish between two sources of zeroes (structural and non-structural) on the basis of a Bayesian framework, using \code{rstan}. All the convergence and goodness-of-fit tests from \code{rstan} are available.
}
\author{
David Moriña (Universitat de Barcelona), Pedro Puig (Universitat Autònoma de Barcelona) and Albert Navarro (Universitat Autònoma de Barcelona)

Mantainer: David Moriña Soler <dmorina@ub.edu>
}

\keyword{ package }
\seealso{
\code{\link{bayesZIB}}
}
\examples{
set.seed(1234)
x <- rbinom(20, 1, 0.4)   # Structural zeroes 
y <- rbinom(20, 1, 0.7*x) # Non-structural zeroes
fit <- bayesZIB(y~1|1, priors=list(c(0,0.5), c(0.5,1)))
print(fit$fit, pars=c("theta", "beta"))
}
