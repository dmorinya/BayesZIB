bayesZIB <- function(formula, data, priors=NULL, chains=3, iter=2000, adapt_delta=0.8, max_treedepth=10, verbose=FALSE, cores=getOption("mc.cores", 1L))
{
  cl <- match.call()
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], 
                                            as.name("|"))) {
    ff <- formula
    formula[[3]][1] <- call("+")
    mf$formula <- formula
    ffc <- . ~ .
    ffz <- ~.
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    ffz[[3]] <- ff[[3]][[3]]
    ffz[[2]] <- NULL
  }
  else {
    ffz <- ffc <- ff <- formula
    ffz[[2]] <- NULL
  }
  if (inherits(try(terms(ffz), silent = TRUE), "try-error")) {
    ffz <- eval(parse(text = sprintf(paste("%s -", deparse(ffc[[2]])), 
                                     deparse(ffz))))
  }
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  mtX <- terms(ffc, data = data)
  X <- model.matrix(mtX, mf)
  mtZ <- terms(ffz, data = data)
  mtZ <- terms(update(mtZ, ~.), data = data)
  Z <- model.matrix(mtZ, mf)
  Y <- model.response(mf, "numeric")
  if (length(Y) < 1) 
    stop("empty model")
  if (all(Y > 0)) 
    stop("invalid dependent variable, minimum count is not zero")
  if (!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 
                                                     0.001))))) 
    stop("invalid dependent variable, non-integer values")
  if (any(Y < 0)) 
    stop("invalid dependent variable, negative counts")
  n <- length(Y)
  kx <- NCOL(X)
  kz <- NCOL(Z)
  if (kx==1 & kz==1 & all(X==1) & all(Z==1))
  {
    if (is.null(priors)) stop("Error: If there are no covariates, priors for w and p should be defined.")
    m=sum(Y)
    llike=function(p,w){
      theta=p*w
      m*log(theta)+(n-m)*log(1-theta)}
    thetaest=m/n
    llmax=m*log(thetaest)+(n-m)*log(1-thetaest)
    
    co=0;tot=0
    
    pr=numeric(iter);wr=numeric(iter)
    while(co<iter){
      lu=log(runif(1))
      
      # The prior is obtained from independent unif rv's 
      
      p0=runif(1,priors[[2]][1],priors[[2]][2]);w0=runif(1,priors[[1]][1],priors[[1]][2])
      
      if(lu<=llike(p0,w0)-llmax){co=co+1; pr[co]=p0; wr[co]=w0}
      tot=tot+1}
    w_mean <- round(mean(wr), 2)
    w_sd   <- round(sd(wr), 2)
    w_sem  <- round(w_sd/sqrt(length(Y)), 2)
    w_p2.5 <- round(quantile(wr, 0.025), 2)
    w_p25  <- round(quantile(wr, 0.25), 2)
    w_p50  <- round(quantile(wr, 0.5), 2)
    w_p75  <- round(quantile(wr, 0.75), 2)
    w_p97.5 <- round(quantile(wr, 0.975), 2)
    p_mean <- round(mean(pr), 2)
    p_sd   <- round(sd(pr), 2)
    p_sem  <- round(p_sd/sqrt(length(Y)), 2)
    p_p2.5 <- round(quantile(pr, 0.025), 2)
    p_p25  <- round(quantile(pr, 0.25), 2)
    p_p50  <- round(quantile(pr, 0.5), 2)
    p_p75  <- round(quantile(pr, 0.75), 2)
    p_p97.5 <- round(quantile(pr, 0.975), 2)
    means <- c(w_mean, p_mean)
    sems <- c(w_sem, p_sem)
    sds <- c(w_sd, p_sd)
    p2.5s <- c(w_p2.5, p_p2.5)
    p25s <- c(w_p25, p_p25)
    p50s <- c(w_p50, p_p50)
    p75s <- c(w_p75, p_p75)
    p97.5s <- c(w_p97.5, p_p97.5)
    tab <- cbind(means, sems, sds, p2.5s, p25s, p50s, 
                 p75s, p97.5s)
    rownames(tab) <- c("w", "p")
    colnames(tab) <- c("mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%")
    print(tab)
    output <- list()
    class(output) <- "bayesZIB"
    attr(output, "Call") <- cl
    attr(output, "x") <- X
    attr(output, "z") <- Z
    output$fit <- list(post_w=wr, post_p=pr)
    return(output)
  }else{
    mcmc_samples <- rstan::sampling(stanmodels$model, 
                                    data = list(N = length(Y), M1 = ncol(X), M2 = ncol(Z), 
                                                X1 = X, X2 = Z, y = Y, s_theta = as.array(rep(1, ncol(X))), 
                                                s = as.array(rep(1, ncol(Z)))), chains = chains, 
                                    iter = iter, control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))
    output <- list()
    class(output) <- c("stanfit", "bayesZIB")
    attr(output, "Call") <- cl
    attr(output, "x") <- X
    attr(output, "z") <- Z
    output$fit <- mcmc_samples
    return(output)
  }
}