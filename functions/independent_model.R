independent_model <- function(df.sample, df.sim, mesh_inf){
  # Building the mesh for the spatial effect
  spde <- INLA::inla.spde2.pcmatern(mesh = mesh_inf, prior.range = c(0.2,0.5), prior.sigma = c(1,0.5))
  A_inf <- INLA::inla.spde.make.A(mesh = mesh_inf, loc = as.matrix(df.sample[,1:2]))
  A_pred <- INLA::inla.spde.make.A(mesh = mesh_inf, loc = as.matrix(df.sim[,1:2]))
  
  inf.stack <- INLA::inla.stack(data = list(y = df.sample$y_sim),
                                A = list(A_inf, 1),
                                effects = list(
                                  list(idx.spde = 1:mesh_inf$n),
                                  list(
                                    intercept = rep(1, times = nrow(df.sample)),
                                    cov = df.sample$cov
                                    )
                                  ),
                                tag = "inf.stack")
  
  pred.stack <- INLA::inla.stack(data = list(y = rep(NA, times = nrow(df.sim))),
                                 A = list(A_pred, 1), 
                                 effects = list(
                                   list(idx.spde = 1:mesh_inf$n),
                                   list(intercept = rep(1, times = nrow(df.sim)),
                                        cov = df.sim$cov)
                                   ),
                                 tag = "pred.stack")
  
  total.stack <- INLA::inla.stack(inf.stack, pred.stack)
  
  ind.formula <-  y ~ -1 + intercept + cov + f(idx.spde, model = spde)
  ind.model <- INLA::inla(data = INLA::inla.stack.data(stack = total.stack), formula = ind.formula, family = "gamma",
                          control.predictor = list(A = INLA::inla.stack.A(stack = total.stack), link = 1),
                          control.compute = list(config = TRUE, waic = FALSE, dic = FALSE), 
                          control.family = list(),
                          inla.mode = "classic"
                          )
  return(list(ind.model=ind.model, total.stack = total.stack))
}
