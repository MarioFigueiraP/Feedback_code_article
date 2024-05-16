# several functions

rmse <- function(y.sim, y.model){
  return(sqrt(mean(y.sim-y.model)**2))
}

mape <- function(y.sim, y.model){
  return(mean(abs((y.sim-y.model)/y.sim))*100)
}

combinatorialFunction <- function(input_list){
  # this function returns a data.frame with the whole set of combinations of the different parameter values
  DF.old <- input_list
  DF.new <- input_list
  for(i in 2:length(input_list)){
    DF.new[[i-1]] <- rep(DF.old[[i-1]], times = length(DF.old[[i]]))
    DF.new[[i]] <- rep(DF.old[[i]], each = length(DF.old[[i-1]]))
    DF.old <- DF.new
  }
  DF <- as.data.frame(DF.new)
  return(DF)
}

mesh.dual <- function(mesh){
  library(sf)
  # Function to build a dual mesh (a Voronoi diagram) of the Constrained Refined Delaunay Triangulation 
  if (mesh$manifold=='R2'){
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0)
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1],
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2,
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      p <- p[order(atan2(yy,xx)),]
      return(rbind(p,p[1,]))
    })
    return(
      lapply(1:mesh$n, function(i)
        st_polygon(x=list(pls[[i]]))
      )
    )
  }
  else stop("It only works for R2!")
}
