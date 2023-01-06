#!/bin/sh

#  RT.sh
#  
#
#  Written by Katherine (Hudson) Gallagher, 2022
#  OVERVIEW: This function is written to calculate the residence time of particles released within the Regional Ocean Modeling System (ROMS). Residence times are calculated as the e-folding time of particle concentration within a defined region or grid. This is the time it takes for particle concentration within a defined region to drop to 1 / natural e.

#   INPUT = matrices with lat/lon positions of particles from ROMS output, number of releases within simulation, known regions to calculate residence time within as a SpatialPolygons object
#   OUTPUT = a matrix where nrow = number of regions and ncol = number of releases

#   NOTE 1: this function was written to be included in a for loop over different release numbers (n). This was done to allow the calculation to be run in parallel within R or as an array on a High Power Computing Cluster and would produce a single vector that was later concantenated with other releases.
#   NOTE 2: the units of residence time will be equal to the frequency that particle positions were saved. For example, if particle positions were saved hourly, residence time will be in hours.


calculateRT <- function(lat, lon, n, nRelease, grid){
    #lon & lat = matrices of particle positions from ROMS. *Note: latitude and longitude are in 2 seperate matrices where rows = number of particles and columns = number of particle positions saved (aka timestep).
    #n = release number
    #nRelease = total number of releases in model run
    #grid = SpatialPolygons object containing the regions or grid cells to calculate residence time within. *Note, the function is written to access SpatialPolygons but any object with accessible lat/lon coordinates will work with minor edits*
  require(sp)
  RT.e <- n.e <- vector(length = length(grid)) #create vector for residence times
  
  release <- seq(n, nrow(lat), by = nRelease) #create sequence of rows for each release
  #particle positions are saved in this matrix as row 1 = particle 1 at release time 1, row 2 = particle 1 at release time 2, etc., so creating the 'release' index allows us to subset to the desired release time
  lon.r <- lon[release,] #isolate particle lat/lons from desired release time
  lat.r <- lat[release,]

for(p in 1:length(grid)){ #for each region/grid cell

  start <- Sys.time() #save start time for benchmarking
  
  poly_point <- vector(length = nrow(lon.r)) #create index vector to determine if particle starts in grid cell
  
  particles_mat <- NULL #create matrix nrow = number of particles in release; ncol = number of timesteps-1
  #will be filled with binary (1 = yes) determining if particles are still in grid
  
  for(v in 1:length(poly_point)){ #for each particle
    lonlat <- as.data.frame(cbind(lon.r[v,], lat.r[v,])) #bind lon/lats into data frame
    w <- which(is.na(lon.r[v,])) #remove NAs from before release & if particle leaves
    lonlat <- lonlat[-w,]
    
    if(nrow(lonlat) > 0){
      poly_point <- point.in.polygon(lonlat[,1], lonlat[,2], grid@polygons[[p]]@Polygons[[1]]@coords[,1], grid@polygons[[p]]@Polygons[[1]]@coords[,2]) #determine if particle starts in grid cell
      
      if(poly_point[1] != 0){ #if the particle starts in the grid cell
          poly <- replace(poly_point, poly_point > 1, 1) #point.in.polygon produces values > 1 if point is on edge/vertex; since we just want presence/absence, replace all values greater than 1 with 1
          if(v != 1){
                  if(length(poly) < length(particles_mat)){ #pad end of poly in case particle leaves model region so that all vectors are the same length as the first one
                      poly <- c(poly, rep(NA, times = (ncol(particles_mat)-length(poly))))
                  }
            } #end if V =! 1
            
          particles_mat <- as.data.frame(rbind(poly, particles_mat)) #create matrix of 1s and 0s (and nas)
          
      } #end if particle starts in polygon
    } #end if nrow(lonlat) > 0
    #print(v)
  } #end v

    if(is.null(particles_mat) == F){
        prop.remaining <- as.vector(colSums(particles_mat, na.rm = T)/nrow(particles_mat)) #sum columns and turn into proportion of particles remaining in region
        
        i <- which(prop.remaining < exp(-1)) #determine when prop drops below 1/e
        
          if(length(i) != 0){
              RT.e[p] <- i[1] #time before prop drops below 1/3 is the RT (check?)
          } else {
              RT.e[p] <- length(prop.remaining) #if i has a length of zero, meaning prop never drops below 1/e, rt is at least the length of the time remaining in the simulation (aka length of prop.remaining)
          } #end length (i)
                  n.e[p] <- nrow(particles_mat) #count how many particles started in grid cell (helpful for troubleshooting as needed)
                  
    } else { #end if particles_mat is.null
      RT.e[p] <- NA #if no particles start in the grid cell, RT = NA
      n.e[p] <- 0
    }
    
### troubleshooting/reporting section###
  end <- Sys.time()
  print(paste('grid cell ', p, ' of ', length(grid), ' done. time elapsed: ', end-start, '. n = ', nrow(particles_mat), sep = '')) #reports time elapsed to calculate residence time in region p, also reports number of particles that started in region
  print(paste("RT = ", RT.e[p], sep = '')) #reports raw residence time
  
} #end p
  RT.n <- as.data.frame(cbind(RT.e, n.e)) #create output data frame
  return(RT.n)
} #end function
