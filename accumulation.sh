#!/bin/sh

#  accumulation.sh
#  
#
#  Written by Katherine (Hudson) Gallagher
#  OVERVIEW: This function is written to calculate the number of simulated particles present within a region/grid with particles from the Regional Ocean Modeling System (ROMS).

#   INPUT = matrices with lat/lon positions of particles from ROMS output**see note**, known regions/grid to calculate accumulation within as a SpatialPolygons object
#   OUTPUT = a matrix where nrow = number of regions and ncol = number of releases

#   NOTE 1: This function will work with raw particle positions from ROMS, however, it was designed to calculate accumulation on particle positions with the first 5 days after release removed to reduce the bias of release position in the results. This was done in a seperate code by converting the first 5 days of positions for each particle to NAs. Any ROMS particle output modified in this or similar manners should work in a similar way.




particleCounts <- function(LON5, LAT5, grid){
    #LON5 & LAT5 = matrices of particle positions from ROMS. *Note: latitude and longitude are in 2 seperate matrices where rows = number of particles and columns = number of particle positions saved (aka timestep).
    #grid = SpatialPolygons object containing the regions or grid cells to calculate residence time within. *Note, the function is written to access SpatialPolygons but any object with accessible lat/lon coordinates will work with minor edits*
    require(sp)
    
    aggs <- matrix(0, nrow = length(grid[[1]]@polygons), ncol = ncol(LAT5)) #create matrix of 0s where rows = how many regions and columns = number of positions saved (aka timesteps)
    
    for(p in 1:length(grid[[1]]@polygons)){ #for each region/grid cells
          start <- Sys.time() #for troubleshooting

        for(x in 1:nrow(LON5)){ #for each particle
              pa <- point.in.polygon(point.x = LON5[x,], point.y = LAT5[x,], pol.x = grid[[1]]@polygons[[p]]@Polygons[[1]]@coords[,1], pol.y = grid[[1]]@polygons[[p]]@Polygons[[1]]@coords[,2]) #determine if particle is in region
              pa <- replace(pa, pa > 1, 1) #point.in.polygon produces values > 1 if point is on edge/vertex; since we just want presence/absence, replace all values greater than 1 with 1
            aggs[p,] <- pa + aggs[p,] #add presence across timesteps to row
        #print(x)
        } #end x
        
        #reporting for troubleshooting/optimization
      end <- Sys.time()
      #print(paste('Region: ', o, ' Cell: ', p, sep = ''))
      print(paste(p, ' time elapsed: ', end-start, sep = ''))
    print(p)
    } #end p
    return(aggs)
} #end function
