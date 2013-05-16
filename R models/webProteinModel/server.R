library(shiny)
library(rgl)

#plot the 3d fit
#input TMD is the nx3 matrix of atoms with x,y,z coordinates
plotTMD <- function(TMD,endpt,xuLim,xdLim,yuLim,ydLim,zuLim,zdLim){
  plot3d(TMD[, 1], TMD[, 2], TMD[, 3], col="blue",
         xlim = c(xdLim,xuLim),ylim = c(ydLim,yuLim),
         zlim = c(zdLim,zuLim),
         xlab="x",ylab="y",zlab="z")
  plot3d(endpt[, 1], endpt[, 2], endpt[, 3], 
         type="l",col="red", lwd=2,add=T)
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

shinyServer(function(input, output) {
  
  #put a title above the graphs that changes when the user selects an item
  #from the dropdown menu
  formulaText <- reactive(function() {
    paste("File:", input$openFile)
  })
  # Return the formula text for printing as a caption
  output$caption <- reactiveText(function(){
    formulaText()
  })
  #plot the data
  output$mainPlot <- reactivePlot(function(){
    allAtoms <- read.csv(file.choose(),head = TRUE,sep=",")
    Ntot <- nrow(allAtoms)
    allAtomsPos <- (matrix(c(allAtoms$x, allAtoms$y, allAtoms$z),
                            nrow=Ntot,ncol=3,byrow=FALSE))
    
    atomsTMD1 <- cbind(allAtoms[c(65:338),1], allAtoms[c(65:338),2], 
                       allAtoms[c(65:338),3])
    TMD1endpts <- fitLine(atomsTMD1)
    plotTMD(atomsTMD1,TMD1endpts,min(allAtomsPos[,1]),max(allAtomsPos[,1]),
                                 min(allAtomsPos[,2]),max(allAtomsPos[,2]),
                                 min(allAtomsPos[,3]),max(allAtomsPos[,3]))
  })
})
