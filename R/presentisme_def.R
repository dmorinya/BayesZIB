library(foreign)
library(rstan)
library(readr)
library(doParallel)
library(bayesZIB)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

nCores   <- detectCores()
registerDoParallel(nCores)

### Prova amb covariants
dades  <- read.dta("../Data/20180629_2.dta")
dades3 <- read.dta("../Data/MOAB_11062020.dta")
dades4 <- dades3[dades3$REGISTRO %in% dades$REGISTRO, ]
dades4$exposure <- dades$exposure
### Recodifico sp a sp_recod
dades4$sp_recod <- ifelse(dades4$sp=="No, never", 0, 1) ### En el cas d'aquestes dades, w=0.39 i p=0.7
dades4$gh <- dades$gh
dades4$sick_abs_days <- dades$sick_abs_days
dades4$age <- dades$age
dades5 <- dades4[, c("sp_recod", "gh", "salary_str", "empl_st2", "age",
                     "remplz", "basic", "unexpect", "dCJLI", "exposure", "sick_abs_days")]

### Plot prior and posteriors (Figure 2)
y <- dades5$sp_recod
### posterior distribution for w
post.w.noConst <- function(w)
{
  return((pbeta(w, sum(y)+1, length(y)-sum(y)+1)-pbeta(w/2, sum(y)+1, length(y)-sum(y)+1))/w)
}

post.w <- function(w)
{
  Cp <- integrate(post.w.noConst, 0, 0.5)$value
  return(post.w.noConst(w)/Cp)
}

### posterior distribution for p
post.p.noConst <- function(p)
{
  return((pbeta(p/2, sum(y)+1, length(y)-sum(y)+1))/p)
}

post.p <- function(p)
{
  Cp <- integrate(post.p.noConst, 0.5, 1)$value
  return(post.p.noConst(p)/Cp)
}

posterior.w.cdf <- function(x)
{
  return(integrate(post.w, 0, x)$value)
}

posterior.p.cdf <- function(x)
{
  return(integrate(post.p, 0.5, x)$value)
}

### Estimates and credible intervals
uniroot(f=function(x){posterior.w.cdf(x)-0.025}, interval=c(0.001, 0.5))$root
uniroot(f=function(x){posterior.w.cdf(x)-0.500}, interval=c(0.001, 0.5))$root
uniroot(f=function(x){posterior.w.cdf(x)-0.975}, interval=c(0.001, 0.5))$root

uniroot(f=function(x){posterior.p.cdf(x)-0.025}, interval=c(0.5, 1))$root
uniroot(f=function(x){posterior.p.cdf(x)-0.500}, interval=c(0.5, 1))$root
uniroot(f=function(x){posterior.p.cdf(x)-0.975}, interval=c(0.5, 1))$root

### Figure 2
par(mfrow=c(2, 2))
### w prior density
x <- seq(0, 0.5, 0.00005)
prior_w <- dunif(x, 0, 0.5)
plot(x, prior_w, type="l", xlab=expression(omega), ylab="Prior density")

### w posterior density
posterior_w <- post.w(x)
plot(x, posterior_w, type="l", xlab=expression(omega), ylab="Posterior density")

### p prior density
x <- seq(0.5, 1, 0.00005)
prior_p <- dunif(x, 0.5, 1)
plot(x, prior_p, type="l", xlab="p", ylab="Prior density")

### p posterior density
posterior_p <- post.p(x)
plot(x, posterior_p, type="l", xlab="p", ylab="Posterior density")
dev.off()

### Model with covariates
set.seed(1234)
fit <- bayesZIB(sp_recod~gh|remplz, data=dades5, chains = 5,
                iter=5000, adapt_delta=0.999, max_treedepth=25)

print(fit$fit, pars=c("theta", "beta"))

### Model without covariates
set.seed(1234)
fit <- bayesZIB(sp_recod~1|1, priors=list(c(0,0.5), c(0.5,1)), data=dades5)

png("diagnostics.png", width=900)
  pairs(fit$fit, pars=c("theta", "beta"))
dev.off()

### Logistic regression models
summary(mod1 <- glm(sp_recod~gh, data=dades5, family="binomial"))
confint(mod1)
summary(mod2 <- glm(sp_recod~remplz, data=dades5[dades5$exposure==1, ], family="binomial"))
confint(mod2)