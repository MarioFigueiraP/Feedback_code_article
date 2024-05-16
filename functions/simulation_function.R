# Simulation of the underlying abundance or biomass ecological phenomenon
simulation_function <- function(seed.global, sim_coord, beta = c(0,1), cov.formula = "cos(x)**2 + sin(y)**2", stdev.gamma = 0.1, rho.spde = 0.2, sigma.spde=1){
  if(!missing(seed.global)) set.seed(seed = seed.global)
  
  # Defining the locations for simulation
  if(!missing(sim_coord)) xy <- sim_coord
  xy <- as.matrix(expand.grid(seq(0,1,length.out=2.5E2), seq(0,1,length.out=2.5E2)))
  
  # Building the mesh for the spatial (SPDE) effect
  mesh <- fmesher::fm_mesh_2d_inla(loc.domain = matrix(c(0,0,0,1,1,1,1,0,0,0), ncol = 2, byrow = TRUE), max.edge = c(0.025,0.05), offset = c(-0.1,-0.15))
  # Simulating the spatial effect at the mesh nodes
  u_spde.nodes <- fmesher::fm_matern_sample(x = mesh, alpha = 2, rho = rho.spde, sigma = sigma.spde)
  # Creating the projection matrix from nodes to the simulation locations
  A_sim <- fmesher::fm_basis(x = mesh, loc = xy)
  # SPDE effect projected at simulation locations
  u_spde.sim <- as.vector(A_sim %*% u_spde.nodes)
  
  # Defining a covariate, e.g. bathymetry temperature, etc.
  f_cov <- function(x,y,formula){eval(parse(text=formula))}
  cov_sim <- f_cov(x=xy[,1], y=xy[,2], formula = cov.formula)
  
  # Defining the gamma generating process
  mu.gamma <- exp(beta[1] + beta[2]*cov_sim + u_spde.sim)
  var.gamma <- stdev.gamma**2
  alpha.gamma <- mu.gamma**2 / var.gamma
  beta.gamma <- mu.gamma / var.gamma
  # Simulation of the response variable (gamma driven)
  y_sim <- rgamma(n = nrow(xy), shape = alpha.gamma, rate = beta.gamma)
  
  return(data.frame(x = xy[,1], y = xy[,2], y_sim = y_sim, cov = cov_sim, u.spde_sim = u_spde.sim))
}
