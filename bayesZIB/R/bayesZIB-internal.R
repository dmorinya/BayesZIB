.__A__ <-
".1"
.__A__.1 <-
function (ns) 
Rcpp::loadModule(module = "stan_fit4model_mod", what = TRUE, 
    env = ns, loadNow = TRUE)
