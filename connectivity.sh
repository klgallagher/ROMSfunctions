#!/bin/sh

#  connectivity.sh
#  
#
#  Written by Katherine (Hudson) Gallagher
#  OVERVIEW: This function is written to calculate two connectivity metrics between different predefined regions using particle trajectories from the Regional Ocean Modeling System (ROMS). The two metrics are:
        #1 - the proportion of particles released in region A interact with all other regions (including region A [which should be equal to 1 as a sanity check])
        #2 - the amount of time particles released in region A spend in all other regions (including region A)

#   INPUT = matrices with lat/lon positions of particles from ROMS output, list of regions where each component of the list is a 2 column matrix with lon/lat coordinates of region boundary, list of indices of particles that start in each regions
#   OUTPUT = 2 lists of results from both metrics (see below)

connectivity <- function(starts, boxes, lon, lat, run, nRelease){
    #starts - a list with the same length as boxes that includes indexes used to subset particle position matrices
    #boxes - list of regions where each component of the list is a 2 column matrix with lon/lat coordinates of region boundary
    #lon & lat = matrices of particle positions from ROMS. *Note: latitude and longitude are in 2 seperate matrices where rows = number of particles and columns = number of particle positions saved (aka timestep).
    #run - a unique identifier associated with the model run of the particles used for saving the results (a character vector)
    #nRelease = total number of releases in model run
    
  boxList <- intList <- vector(mode = 'list', length = length(boxes)) #create empty lists for results
  
  for(x in 1:length(boxes)){ #for each region in boxes
    start <- Sys.time() #start time for optimization/troubleshooting
  
    sInd <- starts[[x]] #grab indexes of which particles start in box
    
    #for each particle release
    perMat <- vector(mode = 'list', length = length(boxes))
    intMat <- NULL
    for(y in 1:nRelease){
      #grab particles released in release y
      r <- seq(y, nrow(lat), by = nRelease)
      latR <- lat[r,]
      lonR <- lon[r,]
      
    #grab particles released within box A
      latA <- latR[sInd,]
      lonA <- lonR[sInd,]
      indNA <- which(is.na(lonA[1,])) #remove NAs from before particle release if any
      latA <- latA[,-c(indNA)]
      lonA <- lonA[,-c(indNA)]

      #create matrices for indexes
      #NOTE: 8 matrices are created here because 8 regions were examined originally. CHANGE FOR YOUR USE CASE!!!
      matA <- matB <- matC <- matD <- matE <- matF <- matG <- matH <- matrix(nrow = length(sInd), ncol = ncol(latA)) #rows = number of particles that start in region of interest, columns = number of particle positions saved (timesteps)
      matList <- list(matA, matB, matC, matD, matE, matF, matG, matH) #put matrices into a single list
      
      #loop through points that start in boxA
      for(z in 1:length(matList)){ #for each matrix
        for(z1 in 1:nrow(latA)){ #for each particle
          matList[[z]][z1,] <- point.in.polygon(point.x = lonA[z1,], point.y = latA[z1,], pol.x = boxes[[z]][,1], pol.y = boxes[[z]][,2]) #get timeseries of presence/absence in region z
          matList[[z]][z1,] <- replace(matList[[z]][z1,], matList[[z]][z1,] > 1, 1) #point.in.polygon produces values > 1 if point is on edge/vertex; since we just want presence/absence, replace all values greater than 1 with 1
        } #end z1
      } #end z
      #end result of this loop is a list of 8 matrices, each filled with 0s and 1s for each particle released in region z interacted with each of the 8 regions (matrix A is filled with presence/absence data within region A of particles released in region A; matrix B is presence/absence data within region B for particles released in region A; and so on.)
      
      #calculate percent of time particle spends in region across all regions/boxes
      matPer <- lapply(matList, FUN = function(q){rowSums(q)/ncol(latA)}) #rowSums gives the number of presences in the region of interest, and dividing by the number of timesteps in latA gives the proportion of time the particle was observed in the region of interest. matPer is a list of vectors where the length of each vector = the number of particles released in the region of interest
      perMat[[1]] <- cbind(perMat[[1]], matPer[[1]])
      perMat[[2]] <- cbind(perMat[[2]], matPer[[2]])
      perMat[[3]] <- cbind(perMat[[3]], matPer[[3]])
      perMat[[4]] <- cbind(perMat[[4]], matPer[[4]])
      perMat[[5]] <- cbind(perMat[[5]], matPer[[5]])
      perMat[[6]] <- cbind(perMat[[6]], matPer[[6]])
      perMat[[7]] <- cbind(perMat[[7]], matPer[[7]])
      perMat[[8]] <- cbind(perMat[[8]], matPer[[8]])
      #perMat[[9]] <- cbind(perMat[[9]], matPer[[9]])
      #By cbinding the vectors, this produces a matrix, where nrow = the number of particles and the ncol = number of particle release events
      
      #calculate percentage of particles released in Box A end up in other boxes (a la Kohut et al 2018)
      intVec <- lapply(matList, FUN = rowSums) #this calculates the number of times a particle was observed in region A
      intVec <- lapply(intVec, FUN = function(q){replace(q,q >= 1, 1)}) #only looking at if particle interacts with region so anything >1 replaced with 1
      intPer <- unlist(lapply(intVec, FUN = function(q){sum(q)/(length(q))})) #sum up how many particles (1s) interacted with region, then divide by the number of particles (length of vector) that started in the region
            #this produces a vector with a length = length(boxes) of the proportion of particles released in region A interacted with all regions for that release event
      intMat <- cbind(intMat, intPer)
        #intMat is a matrix where the nrows = number of regions (length(boxes)) and columns is the number of release events
    } #end y
    
    boxList[[x]] <- abind(perMat, along = 3) #convert list into 3d array with dimensions = number of particles released in region, number of releases, number of regions
    intList[[x]] <- intMat #see above
    
    #timing reports for troubleshooting/optimization
    end <- Sys.time()
    print(x)
    print(end - start)
  } #end x
  
  save(boxList, intList, file = paste('connectivityLists_run', run, '.RData', sep = ''))
} #end function
