require(rgl)
require(onion)

r3dDefaults$windowRect <- c(0,50, 700, 700)

##################Functions
#plot the 3d fit
#input TMD is the nx3 matrix of atoms with x,y,z coordinates
#input endpt is the value returned from fitLine function
#totalTMD is the matrix set of full data read in from a file
#adder is a logical which is T for first run of plotTMD and F otherwise
#     this makes the graph look nice
#col1 is the color of the data points - optional defalut:blue
#col2 is the color of the fitline - optional defalut:red
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
  #need a handle for when the angleBetweenVec returns the supplement of the actual angle
  #   between the vectors
  while(angleBetweenVec(fromCenterendpts,toCenterendpts) > tolerance){
    Q <- 0 + (toCenterendpts[2,1]/2 + fromCenterendpts[2,1]/2)*Hi + (toCenterendpts[2,2]/2 + fromCenterendpts[2,2]/2)*Hj + (toCenterendpts[2,3]/2 + fromCenterendpts[2,3]/2)*Hk
    Q <- Q/sqrt(Norm(Q))
    
    fromCenter <- rotate(fromCenter,Q)
    
    fromCenterendpts <- fitLine(fromCenter)
    
    #check to see if the enpoints were switched during fitLine operation
    #and swap back
    #if these are backward the program hangs
    if(abs(fromCenterendpts[2,1]) < abs(fromCenterendpts[1,1]) && 
         abs(fromCenterendpts[2,2]) < abs(fromCenterendpts[1,2]) &&
         abs(fromCenterendpts[2,3]) < abs(fromCenterendpts[1,3]){
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

#returns the angle in rad between 2 vectors
#vectors are in the form 2x3 double with row 1 being centered a the origin
#or close - machine precision
#V1 is the primary vector
#V2 is the secondary vector
angleBetweenVec <- function(V1,V2){
  return(acos((V1[2,1]*V2[2,1] + V1[2,2]*V2[2,2] + V1[2,3]*V2[2,3])/(normVec(c(V1[2,1],V1[2,2],V1[2,3]))*normVec(c(V2[2,1],V2[2,2],V2[2,3])))))
}

#returns the norm of a vector
#V is in the form c(V_x,V_y,V_z) a 1x3 vector
#V is the vector to norm
normVec <- function(V){
  return(sqrt(V[1]^2+V[2]^2+V[3]^2))
}

#center the tmd and return
#TMD is the matrix of nx3 points
#TMDendpts is the 2x3 matrix of fitted end points - optional but saves recalculation
centerTMD <- function(TMD,TMDendpts){
  if(missing(TMDendpts)) {TMDendpts <- fitline(TMD)}
  return(matrix(c(TMD[,1] - TMDendpts[1,1],TMD[,2] - TMDendpts[1,2],TMD[,3] - TMDendpts[1,3]),ncol = 3))
}

###########code

####test 1 - takes 2 rotations
#rotate 1 to look like 7
#center each one at origin

test1 <- matrix(c(lacyTMD1[,1] - TMD1endpts[1,1],lacyTMD1[,2] - TMD1endpts[1,2],lacyTMD1[,3] - TMD1endpts[1,3]),ncol = 3)
test7 <- matrix(c(lacyTMD7[,1] - TMD7endpts[1,1],lacyTMD7[,2] - TMD7endpts[1,2],lacyTMD7[,3] - TMD7endpts[1,3]),ncol = 3)
testlim <- matrix(c(test1,test7),ncol = 3)

test1endpts <-fitLine(test1)
test7endpts <-fitLine(test7)
plotTMD(test1,test1endpts,testlim,T,"blue")
plotTMD(test7,test7endpts,testlim,F,"skyblue")
rgl.light(specular="#111111",ambient="#777777",viewpoint.rel=FALSE,theta=90,phi=90)
rgl.snapshot(paste("LacY1_7Center",".png",sep=""))

testQ <- 0 + (test7endpts[2,1]/2 + test1endpts[2,1]/2)*Hi + (test7endpts[2,2]/2 + test1endpts[2,2]/2)*Hj + (test7endpts[2,3]/2 + test1endpts[2,3]/2)*Hk
testQ <- testQ/sqrt(Norm(testQ))
test1Rot <- rotate(test1,testQ)
test1Rotendpts <- fitLine(test1Rot)
plotTMD(test1Rot,test1Rotendpts,testlim,T)
plotTMD(test7,test7endpts,testlim,F,"skyblue")
rgl.light(specular="#111111",ambient="#777777",viewpoint.rel=FALSE,theta=90,phi=90)
rgl.snapshot(paste("LacY1_7CenterRot",".png",sep=""))

angleBetweenVec(test1Rotendpts,test7endpts)

#rotate again to match better
testQ <- 0 + (test7endpts[2,1]/2 + test1Rotendpts[2,1]/2)*Hi + (test7endpts[2,2]/2 + test1Rotendpts[2,2]/2)*Hj + (test7endpts[2,3]/2 + test1Rotendpts[2,3]/2)*Hk
testQ <- testQ/sqrt(Norm(testQ))
test1Rota <- rotate(test1Rot,testQ)
test1Rotaendpts <- fitLine(test1Rota)
plotTMD(test1Rota,test1Rotaendpts,testlim,T)
plotTMD(test7,test7endpts,testlim,F,"skyblue")
rgl.light(specular="#111111",ambient="#777777",viewpoint.rel=FALSE,theta=90,phi=90)
rgl.snapshot(paste("LacY1_7CenterRot",".png",sep=""))

angleBetweenVec(test1Rotaendpts,test7endpts)

###test using functional form
test1_7 <- rotateFromTo(lacyTMD1,lacyTMD7)
test1_7endpts <- fitLine(test1_7)
plotTMD(centerTMD(test1_7,test1_7endpts),fitLine(centerTMD(test1_7,test1_7endpts)),lacyPosAtoms,T)
plotTMD(test7,test7endpts,lacyPosAtoms,F,"orange","green")

############actual!!!!
test1_7 <- rotateFromTo(lacyTMD1,lacyTMD7)
test1_7endpts <- fitLine(test1_7)
plotTMD(test1_7,test1_7endpts,lacyPosAtoms,T)
plotTMD(test7,test7endpts,lacyPosAtoms,F,"orange","green")

###problem with 8 -> 2
lacyTMD8Rot <- rotateFromTo(lacyTMD8,lacyTMD2)

plotTMD(lacyTMD8,TMD8endpts,lacyPosAtoms,T,"orange","green")
plotTMD(lacyTMD2,TMD2endpts,lacyPosAtoms,F)

test8 <- matrix(c(lacyTMD8[,1] - TMD8endpts[1,1],lacyTMD8[,2] - TMD8endpts[1,2],lacyTMD8[,3] - TMD8endpts[1,3]),ncol = 3)
test2 <- matrix(c(lacyTMD2[,1] - TMD2endpts[1,1],lacyTMD2[,2] - TMD2endpts[1,2],lacyTMD2[,3] - TMD2endpts[1,3]),ncol = 3)
testlim <- matrix(c(test8,test2),ncol = 3)

test8endpts <-fitLine(test8)
test2endpts <-fitLine(test2)
plotTMD(test8,test8endpts,testlim,T,"orange","green")
plotTMD(test2,test2endpts,testlim,F)

angleBetweenVec(test8endpts,test2endpts)

testQ <- 0 + (test2endpts[2,1]/2 + test8endpts[2,1]/2)*Hi + (test2endpts[2,2]/2 + test8endpts[2,2]/2)*Hj + (test2endpts[2,3]/2 + test8endpts[2,3]/2)*Hk
testQ <- testQ/sqrt(Norm(testQ))
test8Rot <- rotate(test8,testQ)
test8Rotendpts <- fitLine(test8Rot)
plotTMD(test8Rot,test8Rotendpts,testlim,T)
plotTMD(test2,test2endpts,testlim,F)

plotTMD(test8Rot,test8Rotendpts,testlim,T,"orange","green")
plotTMD(test8,test8endpts,testlim,F)

angleBetweenVec(test8Rotendpts,test2endpts)

testQ <- 0 + (test2endpts[2,1]/2 + test8Rotendpts[2,1]/2)*Hi + (test2endpts[2,2]/2 + test8Rotendpts[2,2]/2)*Hj + (test2endpts[2,3]/2 + test8Rotendpts[2,3]/2)*Hk
testQ <- testQ/sqrt(Norm(testQ))
test8Rota <- rotate(test8Rot,testQ)
test8Rotaendpts <- fitLine(test8Rota)
plotTMD(test8Rota,test8Rotaendpts,testlim,T)
plotTMD(test2,test2endpts,testlim,F,"orange","green")

angleBetweenVec(test8Rotaendpts,test2endpts)