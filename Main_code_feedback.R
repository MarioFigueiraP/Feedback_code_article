# Feedbback simulation analysis 

# Loading some needed packages ----

library(INLA)
library(fmesher)
library(ggplot2)
library(gridExtra)
library(ggtext)
library(dplyr)

# Loading custom function to perform the simulations, sample designs, modelization and comparison between modes ----

source("./functions/simulation_function.R")
source("./functions/sampling_functions.R")
source("./functions/independent_model.R")
source("./functions/dependent_model.R")
source("./functions/fb_independent_model.R")
source("./functions/fb_dependent_model.R")
source("./functions/several_functions.R")

cov.formula <- "1.25*x**2 + 0.25*cos(x*pi)**2 + 0.25*sin(y*pi)**2"
range <- c(0.2,0.5,1)
sigma <- c(0.5, 1, 2)
size <- c(1E2, 250, 5E2, 1E3)

A <- list(range = range, sigma = sigma, size = size)
A_DF <- combinatorialFunction(input_list = A)
A_DF <- as.data.frame(lapply(X = A_DF, FUN = rep, 1))

results <- list(ind.rmse = NULL, ind.mape = NULL, dep.rmse = NULL, dep.mape = NULL)

for(i in 1:nrow(A_DF)){
  global_seed <- i
  df_sim <- simulation_function(seed.global = global_seed, beta = c(-1.5,2), cov.formula = cov.formula, rho.spde = A_DF$range[i], sigma.spde = A_DF$sigma[i])
  ind_sample <- ind_sampling(df = df_sim, seed.sample = global_seed, size = A_DF$size[i])
  dep_sample <- dep_sampling(df = df_sim, seed.sample = global_seed, r = 1, size = A_DF$size[i], extra.dependence = extra.dependence, r.ed = 3)
  ggplot() + geom_tile(data = df_sim, mapping = aes(x = x, y = y, fill = y_sim)) + scale_fill_viridis_c(option = "turbo")
  ggplot() + geom_point(data = ind_sample, mapping = aes(x = x, y = y, colour = y_sim)) + scale_color_viridis_c(option = "turbo")
  ggplot() + geom_point(data = dep_sample, mapping = aes(x = x, y = y, colour = y_sim)) + scale_color_viridis_c(option = "turbo")
  
  mesh_inf <- fmesher::fm_mesh_2d_inla(loc.domain = matrix(c(0,0,0,1,1,1,1,0,0,0), ncol = 2, byrow = TRUE), max.edge = c(0.05, 0.10), offset = c(-0.1, -0.15))
  
  ## Independent model ----
  
  ind_model <- independent_model(df.sample = ind_sample, df.sim = df_sim, mesh_inf = mesh_inf)
  idx.pred_ind <- inla.stack.index(stack=ind_model$total.stack, tag="pred.stack")$data
  ggplot() + geom_tile(data=data.frame(x=df_sim$x, y=df_sim$y, z=ind_model$ind.model$summary.fitted.values[idx.pred_ind,"mean"]),
                       mapping=aes(x=x,y=y,fill=z)) + scale_fill_viridis_c(option="turbo")
  results$ind.rmse[i] <- rmse(y.sim = df_sim$y_sim, y.model = ind_model$ind.model$summary.fitted.values[idx.pred_ind,"mean"])
  results$ind.mape[i] <- mape(y.sim = df_sim$y_sim, y.model = ind_model$ind.model$summary.fitted.values[idx.pred_ind,"mean"])
  
  internal.range_ind <- ind_model$ind.model$internal.summary.hyperpar[2, c("mean","sd")]
  internal.sigma_ind <- ind_model$ind.model$internal.summary.hyperpar[3, c("mean","sd")]
  
  internal.gamma.prec_ind <- inla.smarginal(marginal = ind_model$ind.model$internal.marginals.hyperpar$`Intern precision-parameter for the Gamma observations`, log = TRUE)
  prec.gamma.tab_ind <- paste0("table: ",
                           paste(c(internal.gamma.prec_ind$x,
                                   internal.gamma.prec_ind$y),
                                 collapse = " ")
                           )
  
  prior_ind_to_dep <- list(range = internal.range_ind, sigma = internal.sigma_ind,
                           intercept = ind_model$ind.model$summary.fixed["intercept", c("mean", "sd")],
                           cov = ind_model$ind.model$summary.fixed["cov", c("mean", "sd")],
                           prec.gamma = list(prior = prec.gamma.tab_ind))
  
  ## Dependent model ----
  
  dep_model <- dependent_model(df.sample = dep_sample, df.sim = df_sim, mesh_inf = mesh_inf, cov.formula = cov.formula, extra.dependence = extra.dependence, r.ed = 3)
  idx.pred_dep <- inla.stack.index(stack=dep_model$total.stack, tag="pred.stack")$data
  ggplot() + geom_tile(data=data.frame(x=df_sim$x, y=df_sim$y, z=dep_model$dep.model$summary.fitted.values[idx.pred_dep,"mean"]),
                       mapping=aes(x=x,y=y,fill=z)) + scale_fill_viridis_c(option="turbo")
  results$dep.rmse[i] <- rmse(y.sim = df_sim$y_sim, y.model = dep_model$dep.model$summary.fitted.values[idx.pred_dep,"mean"])
  results$dep.mape[i] <- mape(y.sim = df_sim$y_sim, y.model = dep_model$dep.model$summary.fitted.values[idx.pred_dep,"mean"])
  
  internal.range_dep <- dep_model$dep.model$internal.summary.hyperpar[2, c("mean","sd")]
  internal.sigma_dep <- dep_model$dep.model$internal.summary.hyperpar[3, c("mean","sd")]
  
  internal.gamma.prec_dep <- inla.smarginal(marginal = dep_model$dep.model$internal.marginals.hyperpar$`Intern precision-parameter for the Gamma observations`, log = TRUE)
  prec.gamma.tab_dep <- paste0("table: ",
                               paste(c(internal.gamma.prec_dep$x,
                                       internal.gamma.prec_dep$y),
                                     collapse = " ")
  )
  
  prior_dep_to_ind <- list(range = internal.range_dep, sigma = internal.sigma_dep,
                           intercept = ind_model$ind.model$summary.fixed["intercept", c("mean", "sd")],
                           cov = ind_model$ind.model$summary.fixed["cov", c("mean", "sd")],
                           prec.gamma = list(prior = prec.gamma.tab_dep))
  
  # Independent model fed-back ----
  
  fb.ind_model <- fb_independent_model(df.sample = dep_sample, df.sim = df_sim, mesh_inf = mesh_inf, ind_model = ind_model, dep_model = dep_model, prior.fb = prior_dep_to_ind) 
  idx.pred_fb.ind <- inla.stack.index(stack=fb.ind_model$total.stack, tag="pred.stack")$data
  ggplot() + geom_tile(data=data.frame(x=df_sim$x, y=df_sim$y, z=fb.ind_model$ind.model$summary.fitted.values[idx.pred_fb.ind,"mean"]),
                       mapping=aes(x=x,y=y,fill=z)) + scale_fill_viridis_c(option="turbo")
  results$fb.ind.rmse[i] <- rmse(y.sim = df_sim$y_sim, y.model = fb.ind_model$ind.model$summary.fitted.values[idx.pred_fb.ind,"mean"])
  results$fb.ind.mape[i] <- mape(y.sim = df_sim$y_sim, y.model = fb.ind_model$ind.model$summary.fitted.values[idx.pred_fb.ind,"mean"])
  
  # Dependent model fed-back ----
  
  fb.dep_model <- fb_dependent_model(df.sample = dep_sample, df.sim = df_sim, mesh_inf = mesh_inf, cov.formula = cov.formula, extra.dependence = extra.dependence, r.ed = 3, prior.fb = prior_ind_to_dep)
  idx.pred_fb.dep <- inla.stack.index(stack=fb.dep_model$total.stack, tag="pred.stack")$data
  ggplot() + geom_tile(data=data.frame(x=df_sim$x, y=df_sim$y, z=fb.dep_model$dep.model$summary.fitted.values[idx.pred_fb.dep,"mean"]),
                       mapping=aes(x=x,y=y,fill=z)) + scale_fill_viridis_c(option="turbo")
  results$fb.dep.rmse[i] <- rmse(y.sim = df_sim$y_sim, y.model = fb.dep_model$dep.model$summary.fitted.values[idx.pred_fb.dep,"mean"])
  results$fb.dep.mape[i] <- mape(y.sim = df_sim$y_sim, y.model = fb.dep_model$dep.model$summary.fitted.values[idx.pred_fb.dep,"mean"])
  
  }

