#fit data
outputFit <- matrix(nrow = 12,ncol = 2)
  
for(i in 1:12){
  first <- get(paste("TMD",i,"endpts",sep=""))[1,]
  last <- get(paste("TMD",i,"endpts",sep=""))[2,]
    
  outputFit[i,1] <- normVec(first)
  outputFit[i,2] <- normVec(last)
}

#from original data
outputOrig <- matrix(nrow = 12,ncol = 2)

for(i in 1:12){
  first <- cbind(lacyAtoms$x[tmdDomainNum$start[i]],lacyAtoms$y[tmdDomainNum$start[i]],lacyAtoms$z[tmdDomainNum$start[i]])
  last <- cbind(lacyAtoms$x[tmdDomainNum$end[i]],lacyAtoms$y[tmdDomainNum$end[i]],lacyAtoms$z[tmdDomainNum$end[i]])
  
  outputOrig[i,1] <- normVec(first)
  outputOrig[i,2] <- normVec(last)
}

#rotated data
outputRot <- matrix(nrow=12,ncol=2)

for(i in 1:12){
  first <- get(paste("TMD",i,"Rotendpts",sep=""))[1,]
  last <- get(paste("TMD",i,"Rotendpts",sep=""))[2,]
  
  outputRot[i,1] <- normVec(first)
  outputRot[i,2] <- normVec(last)
}