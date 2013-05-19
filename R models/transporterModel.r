#returns the angle in rad between 2 vectors
#vectors are in the form 2x3 double with row 1 being centered a the origin
#or close - machine precision
#V1 is the primary vector
#V2 is the secondary vector
angleBetweenVec <- function(V1,V2){
  return(acos((V1[2,1]*V2[2,1] + V1[2,2]*V2[2,2] + V1[2,3]*V2[2,3])/(normVec(c(V1[2,1],V1[2,2],V1[2,3]))*normVec(c(V2[2,1],V2[2,2],V2[2,3])))))
}

#center the tmd and return
#TMD is the matrix of nx3 points
#TMDendpts is the 2x3 matrix of fitted end points - optional but saves recalculation
centerTMD <- function(TMD,TMDendpts){
  if(missing(TMDendpts)) {TMDendpts <- fitLine(TMD)}
  return(matrix(c(TMD[,1] - TMDendpts[1,1],TMD[,2] - TMDendpts[1,2],TMD[,3] - TMDendpts[1,3]),ncol = 3))
}

#function to fit a line to a 3d set of data
#input X a nx3 matrix of x,y,z coordinates
fitLine <- function(X){
  #get num of rows in X
  N <- nrow(X)
  #average each column and put in a new vector
  meanX <- apply(X, 2, mean)
  #use principle component analysis(pca) to find the best fit line
  Xpca <- princomp(X)
  #some kind of pca thing
  dirVector <- Xpca$loadings[, 1]
  
  #now fit the a line using the pca
  Xfit1 <- matrix(rep(meanX, each=N), ncol=3) + 
    Xpca$score[, 1] %*% t(dirVector)
  t <- c(min(Xpca$score[, 1])-.2, max(Xpca$score[, 1])+.2)
  return(rbind(meanX + t[1]*dirVector, meanX + t[2]*dirVector))
}

#fix for endpts being rotated oppositely
#endpts is the non rotated 2x3 matrix of endpts
#rotEndpts is the 2x3 matrix of endpts from rotateFromTo function
#returns the rotEndpts corrected: based on first row of endpts - both inputs should have the same first row
fixendpts <- function(endpts,rotEndpts){
  tol <- 0.0001
  temp <- matrix(ncol=3,nrow=2)
  if(!(normVec(rotEndpts[1,]) > (normVec(endpts[1,]) - tol) && normVec(rotEndpts[1,]) < (normVec(endpts[1,]) + tol))){
    temp[1,] <- rotEndpts[2,]
    temp[2,] <- rotEndpts[1,]
  }
  else if(!(normVec(rotEndpts[2,]) > (normVec(endpts[2,]) - tol) && normVec(rotEndpts[2,]) < (normVec(endpts[2,]) + tol))){
    temp <- rotEndpts
  }
  else{
    temp <- rotEndpts
  }
  
  return(temp)
}

#returns the norm of a vector
#V is in the form c(V_x,V_y,V_z) a 1x3 vector
#V is the vector to norm
normVec <- function(V){
  return(sqrt(V[1]^2+V[2]^2+V[3]^2))
}

#plot the 3d fit
#input TMD is the nx3 matrix of atoms with x,y,z coordinates
#input endpt is the value returned from fitLine function
#totalTMD is the matrix set of full data read in from a file
#adder is a logical which is T for first run of plotTMD and F otherwise
#     this makes the graph look nice
#datacol is the color of the data points - optional defalut:blue
#fitcol is the color of the fitline - optional defalut:red
plotTMD <- function(TMD,endpt,totalTMD,adder,datacol,fitcol){
  if(missing(adder)) {adder = T}
  if(missing(datacol)) {datacol = "blue"}
  if(missing(fitcol)) {fitcol = "red"}
  
  if(as.logical(adder)){
    plot3d(TMD[, 1], TMD[, 2], TMD[, 3], col=datacol,
           xlim = range(totalTMD[,1]),
           ylim = range(totalTMD[,2]),
           zlim = range(totalTMD[,3]),
           xlab="\305",ylab="\305",zlab="\305",add=F,type="s",size=1)}
  else {
    plot3d(TMD[, 1], TMD[, 2], TMD[, 3], col=datacol,
           xlim = range(totalTMD[,1]),
           ylim = range(totalTMD[,2]),
           zlim = range(totalTMD[,3]),
           xlab="\305",ylab="\305",zlab="\305",add=T,type="s",size=1)
  }
  
  plot3d(endpt[, 1], endpt[, 2], endpt[, 3], 
         type="l",col=fitcol, lwd=2,add=T)
}

#corrects endpoint data to have the [1,] row of data be the first entry in
#the endpt data matrix
#tmdendpts is the set of 2x3 data points representing the fit line from the function
#fitLine
#tmdNum is the tmd number, needed for referencing domainData
#totalData is the entire set of data in a data frame
#domainData is the data frame containg the tmd, start, end data
reorientTMDendpts <- function(tmdendpts,tmdNum,domainData,totalData){
  #define a tolerance for fuzzy searching of begining and end points
  tol <- abs((normVec(tmdendpts[2,]) - normVec(tmdendpts[1,]))/2)
  
  #get the start atom position of the tmd
  first <- cbind(totalData$x[domainData$start[tmdNum]],totalData$y[domainData$start[tmdNum]],
                 totalData$z[domainData$start[tmdNum]])
  
  #check to see if the norm of the first row of the tmd is "close" to the
  #norm of the first row of the tmdendpts vector
  #switch the rows of the tmdendpts if outside the tolerance
  if(!((normVec(first) + tol) > normVec(tmdendpts[1,]) && (normVec(first) - tol) < normVec(tmdendpts[1,])))
  {
    temp <- tmdendpts[1,]
    tmdendpts[1,] <- tmdendpts[2,]
    tmdendpts[2,] <- temp
  }
  
  return(tmdendpts)
}

#rotates a vector to match another
#centers the inputs then rotates the "from" vector to match the "to" vector
#from is the vector to rotate
#to is the vector the from rotates to
#from and to needs to be in the format matrix(c(x[i],y[i],z[i]),ncol=3)
#tolerance is the rotational tolerance to fit to based on the angle between the vectors - optional defalult:0.001
#returns a vector of x,y,z points with the initial endpoint conserved
rotateFromTo <- function(from,to,tolerance) {
  if(missing(tolerance)) {tolerance <- 0.001}
  
  #find the ini endpts
  fromendpts <- fitLine(from)
  toendpts <- fitLine(to)
  
  #center the tmds
  fromCenter <- centerTMD(from,fromendpts)
  toCenter <- centerTMD(to,toendpts)
  
  #fit a line to the centered tmds
  fromCenterendpts <- fitLine(fromCenter)
  toCenterendpts <- fitLine(toCenter)
  
  #rotate from to to based on the average of from and to until tolerance is met
  while(angleBetweenVec(fromCenterendpts,toCenterendpts) > tolerance){
    Q <- 0 + (toCenterendpts[2,1]/2 + fromCenterendpts[2,1]/2)*Hi + 
      (toCenterendpts[2,2]/2 + fromCenterendpts[2,2]/2)*Hj + 
      (toCenterendpts[2,3]/2 + fromCenterendpts[2,3]/2)*Hk
    Q <- Q/sqrt(Norm(Q))
    
    fromCenter <- rotate(fromCenter,Q)
    fromCenterendpts <- fitLine(fromCenter)
    
    #check to see if the enpoints were switched during fitLine operation
    #and swap back
    #if these are backward the program hangs
    if(abs(fromCenterendpts[2,1]) < abs(fromCenterendpts[1,1]) && 
         abs(fromCenterendpts[2,2]) < abs(fromCenterendpts[1,2]) &&
         abs(fromCenterendpts[2,3]) < abs(fromCenterendpts[1,3])){
      for(i in 1:3){
        temp <- fromCenterendpts[2,i]
        fromCenterendpts[2,i] <- fromCenterendpts[1,i]
        fromCenterendpts[1,i] <- temp
      }
    }
  }
  
  #shift the rotated from tmd back to the original ini position
  return(matrix(c(fromCenter[,1] + fromendpts[1,1],fromCenter[,2] + fromendpts[1,2],fromCenter[,3] + fromendpts[1,3]),ncol=3))
}

#write out a tmd to a file with the extra information
writeTMD <- function(tmd,filename,tmdNum,domains,totalAtoms){
  write.table(
    matrix(c(replicate(domains$end[tmdNum]-domains$start[tmdNum]+1,"ATOM"),
             (1:(domains$end[tmdNum]-domains$start[tmdNum]+1)),
             as.character(totalAtoms$name[c(domains$start[tmdNum]:domains$end[tmdNum])]),
             as.character(totalAtoms$AA[c(domains$start[tmdNum]:domains$end[tmdNum])]),
             totalAtoms$AANum[c(domains$start[tmdNum]:domains$end[tmdNum])],
             totalAtoms$x[tmd[,1]],
             totalAtoms$y[tmd[,2]],
             totalAtoms$z[tmd[,3]],
             replicate(domains$end[tmdNum]-domains$start[tmdNum]+1,"1.00"),
             replicate(domains$end[tmdNum]-domains$start[tmdNum]+1,"0.00"),
             as.character(totalAtoms$symbol[c(domains$start[tmdNum]:domains$end[tmdNum])])),
           ncol=11,byrow = FALSE),
    file=filename,row.names=FALSE,quote=FALSE,sep=",",col.names=FALSE)
}