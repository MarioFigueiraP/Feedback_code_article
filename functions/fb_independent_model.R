# This script defines the structure for implementing the independent model with feedback, with covariates, and the spatial effect.
# Therefore, as with the independent model without feedback, part of the script is dedicated to constructing the mesh required for the SPDE-FEM approach to analyze the spatial structure. 
# The necessary procedures to approximate the likelihood of the LGCP process, as described in Simpson et al. (2016), ['Going off grid: computationally efficient inference for log-Gaussian Cox processes'], are also followed."
# Since we are going to implement the feedback procedure, we will have some code related to the redefinition of the prior distributions.

fb_independent_model <- function(df.sample, df.sim, mesh_inf, prior.fb){
  # Defining the SPDE precision matrix structures, following a Matérn family covariance structrue, fixing nu = 1, d = 2, and using normal priors for a reparametrization of rho and sigma with feedback from the preferential model.
  spde <- INLA::inla.spde2.matern(mesh = mesh_inf, alpha = 2, B.tau = matrix(c(0,1,0),1,3), B.kappa = matrix(c(0,0,1),1,3), 
                                  theta.prior.mean = as.numeric(c(prior.fb$range[1], prior.fb$sigma[1])), theta.prior.prec = as.numeric(c(prior.fb$range[2]**-2, prior.fb$sigma[2]**-2)))
  A_inf <- INLA::inla.spde.make.A(mesh = mesh_inf, loc = as.matrix(df.sample[,1:2]))
  A_pred <- INLA::inla.spde.make.A(mesh = mesh_inf, loc = as.matrix(df.sim[,1:2]))
  
  # Defining the data stacks for the data and effects. The stack function, from the INLA package, is used to easily structure complex model information.
  inf.stack <- INLA::inla.stack(data = list(y = df.sample$y_sim),
                                A = list(A_inf, rep(1, times = nrow(df.sample)), df.sample$cov),
                                effects = list(
                                  list(idx.spde = 1:mesh_inf$n),
                                  list(intercept = 1),
                                  list(cov = 1)
                                ),
                                compress = FALSE,
                                remove.unused = FALSE,
                                tag = "inf.stack")
  
  # Defining the stack for the prediction of the abundance/biomass data in the study region.  
  pred.stack <- INLA::inla.stack(data = list(y = rep(NA, times = nrow(df.sim))),
                                 A = list(A_pred, 1), 
                                 effects = list(
                                   list(idx.spde = 1:mesh_inf$n),
                                   list(intercept = rep(1, times = nrow(df.sim)),
                                        cov = df.sim$cov)
                                 ),
                                 tag = "pred.stack")
  
  # Merging both stacks (inference and prediction).
  total.stack <- INLA::inla.stack(inf.stack, pred.stack)
  
  # Formula and inference, where we integrate the prior information form the preferential model to perform the feedback procedure.
  ind.formula <-  y ~ -1 +
    f(intercept, model = "linear", mean.linear = prior.fb$intercept[1], prec.linear = prior.fb$intercept[2]**-2) + 
    f(cov, model = "linear", mean.linear = prior.fb$cov[1], prec.linear = prior.fb$cov[2]**-2) + 
    f(idx.spde, model = spde)
  ind.model <- INLA::inla(data = INLA::inla.stack.data(stack = total.stack), 
                          formula = ind.formula, family = "gamma",
                          control.predictor = list(A = INLA::inla.stack.A(stack = total.stack), link = 1),
                          control.compute = list(config = TRUE, waic = FALSE, dic = FALSE), 
                          # control.family = list(hyper = list(prec=prior.fb$prec.gamma))
                          control.family = list(),
                          inla.mode = "classic"
  )
  return(list(ind.model = ind.model, total.stack = total.stack))
}
