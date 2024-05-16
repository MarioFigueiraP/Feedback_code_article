fb_dependent_model <- function(df.sample, df.sim, mesh_inf, cov.formula = "cos(x)**2 + sin(y)**2", extra.dependence, r.ed = 1, prior.fb){
  # Building the mesh for the spatial effect
  ldomain <- unique(mesh_inf$loc[mesh_inf$segm$int$idx,1:2])
  ldomain <- rbind(ldomain, ldomain[1,])
  dmesh <- mesh.dual(mesh = mesh_inf)
  domainPolygon <- sf::st_polygon(x=list(ldomain))
  w <- sapply(1:length(dmesh), function(i) {
    if (sf::st_intersects(dmesh[[i]], domainPolygon, sparse=FALSE))
      return(sf::st_area(sf::st_intersection(dmesh[[i]], domainPolygon)))
    else return(0)
  })
  
  spde <- INLA::inla.spde2.matern(mesh = mesh_inf, alpha = 2, B.tau = matrix(c(0,1,0),1,3), B.kappa = matrix(c(0,0,1),1,3),
                                  theta.prior.mean = as.numeric(c(prior.fb$range[1], prior.fb$sigma[1])), theta.prior.prec = as.numeric(c(prior.fb$range[2]**-2, prior.fb$sigma[2]**-2)))
  A_inf_geo <- INLA::inla.spde.make.A(mesh = mesh_inf, loc = as.matrix(df.sample[,1:2]))
  A_inf_mesh <- Diagonal(n = mesh_inf$n, x = 1)
  A_pred <- INLA::inla.spde.make.A(mesh = mesh_inf, loc = as.matrix(df.sim[,1:2]))
  
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
  
  f_cov <- function(x,y,formula){eval(parse(text=formula))}
  cov_mesh <- f_cov(x=mesh_inf$loc[,1], y=mesh_inf$loc[,2], formula = cov.formula)
  
  df.extra.dep <- rbind(data.frame(x=mesh_inf$loc[,1], y=mesh_inf$loc[,2]), df.sample[,1:2])
  cond <- if(missing(extra.dependence)){0} else{do.call(what = extra.dependence, args = list(df=df.extra.dep, r.ed=r.ed))}
  
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
  
  inf.stack_dep <- INLA::inla.stack(data = list(y = cbind(NA, c(rep(0, mesh_inf$n), rep(1, nrow(df.sample))), NA), e = c(w, rep(0, nrow(df.sample)))),
                                    A = list(1,1),
                                    effects = list(
                                      list(iid.copied = 1:(mesh_inf$n + nrow(df.sample))),
                                      list(extra.dep = cond)
                                    ),
                                    tag = "inf.stack_dep")
  
  pred.stack_geo <- INLA::inla.stack(data = list(y = matrix(NA, nrow = nrow(df.sim), ncol = 3)),
                                     A = list(A_pred, 1),
                                     effects = list(
                                       list(idx.spde = 1:mesh_inf$n),
                                       list(intercept = rep(1, times = nrow(df.sim)),
                                            cov = df.sim$cov)
                                     ),
                                     tag = "pred.stack")
  
  total.stack <- INLA::inla.stack(inf.stack_geo, inf.stack_copy, inf.stack_dep, pred.stack_geo)
  
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
