arrgh <- vector(length = 12)
numbers <- matrix(nrow = 12,ncol=3,dimnames=list(c(1:12),c("+tol","-tol","test")))
for(i in 1:12){
  placeH2 <- get(paste("TMD",i,"endpts",sep=""))[2,]
  placeH1 <- get(paste("TMD",i,"endpts",sep = ""))[1,]
  
  testtol <- abs((normVec(placeH2) - normVec(placeH1))/2)
  testFirst <- cbind(lacyAtoms$x[tmdDomainNum$start[i]],lacyAtoms$y[tmdDomainNum$start[i]],lacyAtoms$z[tmdDomainNum$start[i]])
  
  numbers[i,1] <- (normVec(testFirst) + testtol)
  numbers[i,2] <- (normVec(testFirst) - testtol)
  numbers[i,3] <- normVec(placeH1)
  
  arrgh[i] <- !((normVec(testFirst) + testtol) > normVec(placeH1) && (normVec(testFirst) - testtol) < normVec(placeH1))
}