library(bayesZIB)
library(doParallel)
library(WriteXLS)
options(mc.cores = parallel::detectCores())

nCores   <- detectCores()
registerDoParallel(nCores)

nsim <- 100

for (k in 1:nsim)
{
  n <- c(500, 1500) # Sample sizes
  theta0 <- c(-0.5, -1, -2)
  theta1 <- c(-2, -3, -4)
  theta2 <- -3
  beta0  <- c(0.5, 1, 2)
  beta1  <- c(2, 3, 4)
  beta2  <- 3
  resultat <- data.frame(expand.grid(n=n, theta0=theta0, theta1=theta1, theta2=theta2,
                                     beta0=beta0, beta1=beta1, beta2=beta2))
                                    
  genEsts <- function(i)
  {
    print(paste0("Simulation step ", i, " out of ", dim(resultat)[1]))
    n_sim <- resultat$n[i]
    theta0_sim <- resultat$theta0[i]
    theta1_sim <- resultat$theta1[i]
    theta2_sim <- resultat$theta2[i]
    beta0_sim <- resultat$beta0[i]
    beta1_sim <- resultat$beta1[i]
    beta2_sim <- resultat$beta2[i]

    ### Structural zeros (exposure)
    x1 <- rnorm(n_sim)
    x2 <- rnorm(n_sim)
    z1 <- theta0_sim + theta1_sim*x1 + theta2_sim*x2
    pr <- 1/(1+exp(-z1)) 
    exposure <- rbinom(n_sim, 1, pr)

    ### Sample zeros (observed phenomenon)
    x3  <- rnorm(n_sim)
    x4  <- rnorm(n_sim)
    z2  <- beta0_sim + beta1_sim*x3 + beta2_sim*x4
    pr2 <- 1/(1+exp(-z2)) 
    p    <- rbinom(n_sim, 1, pr2)
    pres <- ifelse(exposure==0, 0, p)
    df <- data.frame(exposure=exposure, pres=pres)

    df$x1 <- x1; df$x2 <- x2; df$x3 <- x3; df$x4 <- x4
    fit <- bayesZIB(pres~x1+x2|x3+x4, data=df, chains = 5,
                    adapt_delta=0.999, max_treedepth=15, verbose=FALSE)
    post_matrix  <- as.matrix(fit$fit, pars=c("theta", "beta"))
    theta0_p50   <- median(post_matrix[, 1])
    theta0_p2.5  <- quantile(post_matrix[, 1], 0.025)
    theta0_p97.5 <- quantile(post_matrix[, 1], 0.975)
    theta1_p50   <- median(post_matrix[, 2])
    theta1_p2.5  <- quantile(post_matrix[, 2], 0.025)
    theta1_p97.5 <- quantile(post_matrix[, 2], 0.975)
    theta2_p50   <- median(post_matrix[, 3])
    theta2_p2.5  <- quantile(post_matrix[, 3], 0.025)
    theta2_p97.5 <- quantile(post_matrix[, 3], 0.975)
    beta0_p50    <- median(post_matrix[, 4])
    beta0_p2.5   <- quantile(post_matrix[, 4], 0.025)
    beta0_p97.5  <- quantile(post_matrix[, 4], 0.975)
    beta1_p50    <- median(post_matrix[, 5])
    beta1_p2.5   <- quantile(post_matrix[, 5], 0.025)
    beta1_p97.5  <- quantile(post_matrix[, 5], 0.975)
    beta2_p50    <- median(post_matrix[, 6])
    beta2_p2.5   <- quantile(post_matrix[, 6], 0.025)
    beta2_p97.5  <- quantile(post_matrix[, 6], 0.975)
    return(c(resultat$n[i], resultat$theta0[i], resultat$theta1[i], resultat$theta2[i], 
             resultat$beta0[i], resultat$beta1[i], resultat$beta2[i], 
             theta0_p2.5, theta0_p50, theta0_p97.5, theta1_p2.5, theta1_p50, theta1_p97.5, theta2_p2.5, theta2_p50, theta2_p97.5,
             beta0_p2.5, beta0_p50, beta0_p97.5, beta1_p2.5, beta1_p50, beta1_p97.5, beta2_p2.5, beta2_p50, beta2_p97.5))
  }

  system.time(dat.fin <- foreach(j=1:dim(resultat)[1], .combine=rbind) %dopar% genEsts(j))
  colnames(dat.fin) <- c("n", "theta0", "theta1", "theta2", "beta0", "beta1", "beta2", 
                         "theta0_p2.5", "theta0_p50", "theta0_p97.5",
                         "theta1_p2.5", "theta1_p50", "theta1_p97.5",
                         "theta2_p2.5", "theta2_p50", "theta2_p97.5",
                         "beta0_p2.5", "beta0_p50", "beta0_p97.5",
                         "beta1_p2.5", "beta1_p50", "beta1_p97.5",
                         "beta2_p2.5", "beta2_p50", "beta2_p97.5")

  ### Excel exportation
  file_name <- paste0("../Results/simZIB_covs", k, ".xls")
  WriteXLS(as.data.frame(dat.fin), file_name)
}