# Simulation of the independent sampling, independent stratified sampling and dependent sampling

ind_sampling <- function(df, size, seed.sample){
  if(!missing(seed.sample)) set.seed(seed.sample)
  
  # Independent and completly random sample from the simulated response variable data.frame
  sample <- df[sample(x = 1:nrow(df), size = size),]
  return(sample)
}

ind_stratified_sampling <- function(seed.sample, size.per.cell, n.cells.dim, ...){
  if(!missing(seed.sample)) set.seed(seed.sample)
  
  # Building the grid stratification for the sample process
  area_length <- 1; delta_length <- 0.001; numb_rec_by_axis <- n.cells.dim # Total number of rectangles n.cells.dim*n.cells.dim (numb_rec_by_axis). The delta_length is a delta extension to avoid Na's in the over operation between points and areas
  coord_agg_grid <- expand.grid(x=seq(-delta_length/2+area_length/numb_rec_by_axis/2, 1+delta_length/2-area_length/numb_rec_by_axis/2, length.out=numb_rec_by_axis),
                                y=seq(-delta_length/2+area_length/numb_rec_by_axis/2, 1+delta_length/2-area_length/numb_rec_by_axis/2, length.out=numb_rec_by_axis))
  spcoord <- sp::SpatialPoints(coord_agg_grid)
  sppixels <- sp::SpatialPixels(spcoord)
  PolygonsPixels <- as(sppixels, Class="SpatialPolygons")
  GridPolygons <- sf::st_as_sf(PolygonsPixels)
  
  # Locations simulated from each polygon of the grid 
  xy_sample <- do.call(rbind, lapply(X = 1:nrow(GridPolygons), FUN = function(x){
    sf::st_coordinates(sf::st_sample(x = GridPolygons[x,], size = size.per.cell, type = "random"))
    })
  )
  
  # Simulation of the response variable in the above defined locations
  sample <- simulation_function(sim_coord = xy_sample, ...)
  return(sample)
}

extra.dependence <- function(df, r.ed){
  # An extra dependence in the dependent (or preferential) sampling related to the closest distance between sample points and some fixed points, e.g. the distance between ships and the closest port
  cond_Matrix <- do.call(rbind, lapply(X = list(sf::st_point(c(0,0)), sf::st_point(c(0,1)), sf::st_point(c(1,1))), FUN = function(x){
    sf::st_distance(x, sf::st_as_sf(df, coords=c(1,2)) %>% sf::st_geometry())
    })
  )
  cond <- -r.ed*apply(X = cond_Matrix, MARGIN = 2, FUN = min)
  return(cond)
}

dep_sampling <- function(df, r, seed.sample, size, extra.dependence, r.ed=1){
  if(!missing(seed.sample)) set.seed(seed.sample)

  # Dependent sampling simulation from the simulated response variable data.frame
  cond <- if(missing(extra.dependence)){0} else{do.call(what = extra.dependence, args = list(df=df, r.ed=r.ed))}
  prob <- exp(r*log(df[,3]) + cond)#/sum(exp(r*log(df[,3]) + cond)) # Defining the probability of being sample for each point
  
  idx <- sample(x = 1:nrow(df), size = size, prob = prob)
  sample <- df[idx,]
  sample$cond <- cond[idx]
  return(sample)
}
