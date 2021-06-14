\encoding{utf8} 
\name{bayesZIB}
\alias{bayesZIB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Bernoulli zero-inflated regression model.
}
\description{
Fit Bernoulli zero-inflated regression models in a Bayesian framework.
}
\usage{
bayesZIB(formula, data, priors=NULL, chains=3, iter=2000, 
         adapt_delta=0.8, max_treedepth=10, verbose=FALSE, 
         cores=getOption("mc.cores", 1L))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
symbolic description of the model, see details.
}
  \item{data}{
arguments controlling formula processing via \code{model.frame}.
}
  \item{priors}{
\code{list} with two elements specifying the limits of the uniform priors for \code{w} and \code{p} respectively. It is \code{NULL} by default but should be defined if there are no covariates.
}
  \item{chains}{
a positive integer specifying the number of Markov chains. The default is 3.
}
  \item{iter}{
a positive integer specifying the number of iterations for each chain (including warmup). The default is 2000.
}
  \item{adapt_delta}{
for the No-U-Turn Sampler (NUTS), the variant of Hamiltonian Monte Carlo used used by \code{rstan}, adapt_delta is the target average proposal acceptance probability for adaptation. double, between 0 and 1, defaults to 0.8.
}
  \item{max_treedepth}{
maximum depth parameter. Positive integer, defaults to 10. When the maximum allowed tree depth is reached it indicates that NUTS is terminating prematurely to avoid excessively long execution time.
}
  \item{verbose}{
\code{TRUE} or \code{FALSE}: flag indicating whether to print intermediate output from Stan on the console, which might be helpful for model debugging.
}
  \item{cores}{
number of cores to use when executing the chains in parallel, which defaults to 1 but according to the Stan documentation it is recommended to set the \code{mc.cores} option to be as many processors as the hardware and RAM allow (up to the number of chains).
}
}

\details{
Zero-inflated models are two-component mixture models combining a point mass at zero with a proper count distribution. Thus, there are two sources of zeros: zeros may come from both the point mass and from the Bernoulli component. For modeling the unobserved state (zero vs. Bernoulli), a binary model is used that captures the probability of zero inflation. in the simplest case only with an intercept but potentially containing regressors. For this zero-inflation model, a binomial model with an appropriate link function is used.

The formula can be used to specify both components of the model: If a formula of type \code{y ~ x1 + x2} is supplied, then the same regressors are employed in both components. This is equivalent to \code{y ~ x1 + x2 | x1 + x2}. Of course, a different set of regressors could be specified for the Bernoulli and zero-inflation component, e.g., \code{y ~ x1 + x2 | z1 + z2 + z3} giving the logistic regression model \code{y ~ x1 + x2} conditional on (\code{|}) the zero-inflation model \code{y ~ z1 + z2 + z3}. A simple inflation model where all zero counts have the same probability of belonging to the zero component can by specified by the formula \code{y ~ x1 + x2 | 1}.
}

\value{
An object of class "bayesZIB", i.e., a list with components including
\item{Call}{text string with the original call to the function}
\item{x}{design matrix for the zero-inflated part}
\item{z}{design matrix for the non zero-inflated part}
\item{fit}{an object of S4 class \code{stanfit} if there are covariates or a named list with \code{iter} draws from the posterior distribution of \code{w} and \code{p}.}
}

\author{
David Moriña (Universitat de Barcelona), Pedro Puig (Universitat Autònoma de Barcelona) and Albert Navarro (Universitat Autònoma de Barcelona)

Mantainer: David Moriña Soler <dmorina@ub.edu>
}

\seealso{
\code{\link{bayesZIB-package}}
}
\examples{
set.seed(1234)
x <- rbinom(20, 1, 0.4)   # Structural zeroes 
y <- rbinom(20, 1, 0.7*x) # Non-structural zeroes
fit <- bayesZIB(y~1|1, priors=list(c(0, 0.5), c(0.5, 1)))
print(fit$fit, pars=c("theta", "beta"))
}

\keyword{ bayesZIB }
