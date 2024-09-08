# This script defines the structure for implementing the preferential model with feedback, with covariates, spatial effects, and taking into account the effect of the point pattern structure of the sampling.
# Therefore, as with the preferential model without feedback, part of the script is dedicated to constructing the mesh required for the SPDE-FEM approach to analyze the spatial structure. 
# The necessary procedures to approximate the likelihood of the LGCP process, as described in Simpson et al. (2016), ['Going off grid: computationally efficient inference for log-Gaussian Cox processes'], are also followed."
# Since we are going to implement the feedback procedure, we will have some code related to the redefinition of the prior distributions.

fb_dependent_model <- function(df.sample, df.sim, mesh_inf, cov.formula = "cos(x)**2 + sin(y)**2", extra.dependence, r.ed = 1, prior.fb){
  # Defining the vector of weights (w), which is needed as part of the integration in the non-homogeneous Poisson likelihood.
  ldomain <- unique(mesh_inf$loc[mesh_inf$segm$int$idx,1:2])
  ldomain <- rbind(ldomain, ldomain[1,])
  dmesh <- mesh.dual(mesh = mesh_inf)
  domainPolygon <- sf::st_polygon(x=list(ldomain))
  w <- sapply(1:length(dmesh), function(i) {
    if (sf::st_intersects(dmesh[[i]], domainPolygon, sparse=FALSE))
      return(sf::st_area(sf::st_intersection(dmesh[[i]], domainPolygon)))
    else return(0)
  })
  
  # Defining the SPDE precision matrix structures, following a Matérn family covariance structrue, fixing nu = 1, d = 2, and using normal priors for a reparametrization of rho and sigma with feedback from the independent model.
  spde <- INLA::inla.spde2.matern(mesh = mesh_inf, alpha = 2, B.tau = matrix(c(0,1,0),1,3), B.kappa = matrix(c(0,0,1),1,3),
                                  theta.prior.mean = as.numeric(c(prior.fb$range[1], prior.fb$sigma[1])), theta.prior.prec = as.numeric(c(prior.fb$range[2]**-2, prior.fb$sigma[2]**-2)))
  A_inf_geo <- INLA::inla.spde.make.A(mesh = mesh_inf, loc = as.matrix(df.sample[,1:2]))
  A_inf_mesh <- Diagonal(n = mesh_inf$n, x = 1)
  A_pred <- INLA::inla.spde.make.A(mesh = mesh_inf, loc = as.matrix(df.sim[,1:2]))
  
  # Defining the data stacks for the data and effects. The stack function, from the INLA package, is used to easily structure complex model information.
  inf.stack_geo <- INLA::inla.stack(data = list(y = cbind(df.sample$y_sim, NA, NA)),
                                    A = list(A_inf_geo, 1),
                                    effects = list(
                                      list(idx.spde = 1:mesh_inf$n),
                                      list(
                                        intercept = rep(1, times = nrow(df.sample)),
                                        cov = df.sample$cov
                                      )
                                    ),
                                    tag = "inf.stack_geo")
  
  # Defining the covariate values on the mesh nodes, which is required for the integration of part of the non-homogeneous Poisson likelihood.
  f_cov <- function(x,y,formula){eval(parse(text=formula))}
  cov_mesh <- f_cov(x=mesh_inf$loc[,1], y=mesh_inf$loc[,2], formula = cov.formula)
  
  # This is used for the extra spatial dependence in the point pattern of the samples.
  df.extra.dep <- rbind(data.frame(x=mesh_inf$loc[,1], y=mesh_inf$loc[,2]), df.sample[,1:2])
  cond <- if(missing(extra.dependence)){0} else{do.call(what = extra.dependence, args = list(df=df.extra.dep, r.ed=r.ed))}
  
  # The following are the other stacks needed for completely define the preferential model: 
  # the stack for the shared part between the abundance/biomass data and the point pattern; 
  inf.stack_copy <- INLA::inla.stack(data = list(y = cbind(NA, NA, rep(0, mesh_inf$n + nrow(df.sample)))),
                                     A = list(rbind(A_inf_mesh, A_inf_geo), 1, -1),
                                     effects = list(
                                       list(idx.spde = 1:mesh_inf$n),
                                       list(
                                         intercept = rep(1, mesh_inf$n + nrow(df.sample)),
                                         cov = c(cov_mesh, df.sample$cov)),
                                       list(iid.copy = 1:(mesh_inf$n + nrow(df.sample)))
                                     ),
                                     tag = "copy.stack"
  )
  # the stack for the preferential part of the model (the LGCP);
  inf.stack_dep <- INLA::inla.stack(data = list(y = cbind(NA, c(rep(0, mesh_inf$n), rep(1, nrow(df.sample))), NA), e = c(w, rep(0, nrow(df.sample)))),
                                    A = list(1,1),
                                    effects = list(
                                      list(iid.copied = 1:(mesh_inf$n + nrow(df.sample))),
                                      list(extra.dep = cond)
                                    ),
                                    tag = "inf.stack_dep")
  # the stack for the prediction of the abundance/biomass data in the study region;
  pred.stack_geo <- INLA::inla.stack(data = list(y = matrix(NA, nrow = nrow(df.sim), ncol = 3)),
                                     A = list(A_pred, 1),
                                     effects = list(
                                       list(idx.spde = 1:mesh_inf$n),
                                       list(intercept = rep(1, times = nrow(df.sim)),
                                            cov = df.sim$cov)
                                     ),
                                     tag = "pred.stack")
  
  total.stack <- INLA::inla.stack(inf.stack_geo, inf.stack_copy, inf.stack_dep, pred.stack_geo)
  # formula and inference, where we integrate the prior information form the independent model to perform the feedback procedure.
  dep.formula <-  y ~ -1 + 
    f(intercept, model = "linear", mean.linear = prior.fb$intercept[1], prec.linear = prior.fb$intercept[2]**-2) + 
    f(cov, model = "linear", mean.linear = prior.fb$cov[1], prec.linear = prior.fb$cov[2]**-2) + 
    f(idx.spde, model = spde) +
    f(iid.copy, model = "iid", hyper = list(prec = list(initial = -10, fixed = TRUE))) + extra.dep + f(iid.copied, copy = "iid.copy", fixed = FALSE)
  dep.model <- INLA::inla(data = INLA::inla.stack.data(stack = total.stack),
                          formula = dep.formula, family = c("gamma", "poisson", "gaussian"),
                          E = INLA::inla.stack.data(stack = total.stack)$e,
                          control.predictor = list(A = INLA::inla.stack.A(stack = total.stack), link = 1),
                          control.compute = list(config = TRUE, waic = FALSE, dic = FALSE),
                          # control.family = list(list(hyper = list(prec=prior.fb$prec.gamma)), list(), list(hyper = list(prec=list(initial = 10, fixed = TRUE))))
                          control.family = list(list(), list(), list(hyper = list(prec=list(initial = 10, fixed = TRUE)))),
                          inla.mode = "classic"
  )
  return(list(dep.model=dep.model, total.stack = total.stack))
}
