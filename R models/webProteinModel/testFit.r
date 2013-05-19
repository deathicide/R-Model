require(devtools)
require(rgl)
require(onion)

r3dDefaults$windowRect <- c(0,0, 1280, 756)

#rgl.snapshot(paste("pic",".png",sep=""))
#rgl.viewpoint(fov = 30,zoom = 0.9)

#load all my functions
setwd("~/Documents/School Work/Grad School/Spring 2013/Research/R models/")
source("transporterModel.r")

#read in the data from a file
lacyAtoms <- read.csv(file = "webProteinModel/LacYAtoms.csv",head = TRUE,sep=",")
tmdDomainNum <- read.csv(file = "webProteinModel/tmdDomains.csv",head = TRUE,sep=",")

Ntot <- nrow(lacyAtoms)

#put the x,y,and z position data columns into a new matrix
lacyPosAtoms <- (matrix(c(lacyAtoms$x, lacyAtoms$y, lacyAtoms$z),
                        nrow=Ntot,ncol=3,byrow=FALSE))

#now pick out a tmd
{
lacyTMD1 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[1]:tmdDomainNum$end[1]),1], 
                  lacyPosAtoms[c(tmdDomainNum$start[1]:tmdDomainNum$end[1]),2], 
                  lacyPosAtoms[c(tmdDomainNum$start[1]:tmdDomainNum$end[1]),3])
lacyTMD2 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[2]:tmdDomainNum$end[2]),1], 
                  lacyPosAtoms[c(tmdDomainNum$start[2]:tmdDomainNum$end[2]),2], 
                  lacyPosAtoms[c(tmdDomainNum$start[2]:tmdDomainNum$end[2]),3])
lacyTMD3 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[3]:tmdDomainNum$end[3]),1], 
                  lacyPosAtoms[c(tmdDomainNum$start[3]:tmdDomainNum$end[3]),2], 
                  lacyPosAtoms[c(tmdDomainNum$start[3]:tmdDomainNum$end[3]),3])
lacyTMD4 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[4]:tmdDomainNum$end[4]),1], 
                  lacyPosAtoms[c(tmdDomainNum$start[4]:tmdDomainNum$end[4]),2], 
                  lacyPosAtoms[c(tmdDomainNum$start[4]:tmdDomainNum$end[4]),3])
lacyTMD5 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[5]:tmdDomainNum$end[5]),1], 
                  lacyPosAtoms[c(tmdDomainNum$start[5]:tmdDomainNum$end[5]),2], 
                  lacyPosAtoms[c(tmdDomainNum$start[5]:tmdDomainNum$end[5]),3])
lacyTMD6 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[6]:tmdDomainNum$end[6]),1], 
                  lacyPosAtoms[c(tmdDomainNum$start[6]:tmdDomainNum$end[6]),2], 
                  lacyPosAtoms[c(tmdDomainNum$start[6]:tmdDomainNum$end[6]),3])
lacyTMD7 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[7]:tmdDomainNum$end[7]),1], 
                  lacyPosAtoms[c(tmdDomainNum$start[7]:tmdDomainNum$end[7]),2], 
                  lacyPosAtoms[c(tmdDomainNum$start[7]:tmdDomainNum$end[7]),3])
lacyTMD8 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[8]:tmdDomainNum$end[8]),1], 
                  lacyPosAtoms[c(tmdDomainNum$start[8]:tmdDomainNum$end[8]),2], 
                  lacyPosAtoms[c(tmdDomainNum$start[8]:tmdDomainNum$end[8]),3])
lacyTMD9 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[9]:tmdDomainNum$end[9]),1], 
                  lacyPosAtoms[c(tmdDomainNum$start[9]:tmdDomainNum$end[9]),2], 
                  lacyPosAtoms[c(tmdDomainNum$start[9]:tmdDomainNum$end[9]),3])
lacyTMD10 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[10]:tmdDomainNum$end[10]),1], 
                   lacyPosAtoms[c(tmdDomainNum$start[10]:tmdDomainNum$end[10]),2], 
                   lacyPosAtoms[c(tmdDomainNum$start[10]:tmdDomainNum$end[10]),3])
lacyTMD11 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[11]:tmdDomainNum$end[11]),1], 
                   lacyPosAtoms[c(tmdDomainNum$start[11]:tmdDomainNum$end[11]),2], 
                   lacyPosAtoms[c(tmdDomainNum$start[11]:tmdDomainNum$end[11]),3])
lacyTMD12 <- cbind(lacyPosAtoms[c(tmdDomainNum$start[12]:tmdDomainNum$end[12]),1], 
                   lacyPosAtoms[c(tmdDomainNum$start[12]:tmdDomainNum$end[12]),2], 
                   lacyPosAtoms[c(tmdDomainNum$start[12]:tmdDomainNum$end[12]),3])
}

#make a fit line for all tmds
{
TMD1endpts <-fitLine(lacyTMD1)
TMD2endpts <-fitLine(lacyTMD2)
TMD3endpts <-fitLine(lacyTMD3)
TMD4endpts <-fitLine(lacyTMD4)
TMD5endpts <-fitLine(lacyTMD5)
TMD6endpts <-fitLine(lacyTMD6)
TMD7endpts <-fitLine(lacyTMD7)
TMD8endpts <-fitLine(lacyTMD8)
TMD9endpts <-fitLine(lacyTMD9)
TMD10endpts <-fitLine(lacyTMD10)
TMD11endpts <-fitLine(lacyTMD11)
TMD12endpts <-fitLine(lacyTMD12)
}

#plot all tmd with fit lines
{
plotTMD(lacyTMD1,TMD1endpts,lacyPosAtoms,T,"blue")
plotTMD(lacyTMD2,TMD2endpts,lacyPosAtoms,F,"darkorange")
plotTMD(lacyTMD3,TMD3endpts,lacyPosAtoms,F,"red")
plotTMD(lacyTMD4,TMD4endpts,lacyPosAtoms,F,"yellow")
plotTMD(lacyTMD5,TMD5endpts,lacyPosAtoms,F,"green")
plotTMD(lacyTMD6,TMD6endpts,lacyPosAtoms,F,"violetred")
plotTMD(lacyTMD7,TMD7endpts,lacyPosAtoms,F,"skyblue")
plotTMD(lacyTMD8,TMD8endpts,lacyPosAtoms,F,"orange")
plotTMD(lacyTMD9,TMD9endpts,lacyPosAtoms,F,"indianred1")
plotTMD(lacyTMD10,TMD10endpts,lacyPosAtoms,F,"lightyellow")
plotTMD(lacyTMD11,TMD11endpts,lacyPosAtoms,F,"palegreen")
plotTMD(lacyTMD12,TMD12endpts,lacyPosAtoms,F,"palevioletred")

rgl.light(specular="#111111",ambient="#777777",viewpoint.rel=FALSE,theta=90,phi=90)
rgl.snapshot(paste("LacYTMD1-2",".png",sep=""))

rgl.light(specular="#111111",ambient="#777777",viewpoint.rel=FALSE,theta=90,phi=90)
rgl.snapshot(paste("LacYAll2",".png",sep=""))
}

#write out each tmd
{
writeTMD(lacyTMD1,"TMD1.txt",1,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD2,"TMD2.txt",2,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD3,"TMD3.txt",3,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD4,"TMD4.txt",4,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD5,"TMD5.txt",5,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD6,"TMD6.txt",6,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD7,"TMD7.txt",7,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD8,"TMD8.txt",8,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD9,"TMD9.txt",9,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD10,"TMD10.txt",10,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD11,"TMD11.txt",11,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD12,"TMD12.txt",12,tmdDomainNum,lacyAtoms)
}

#reorder the endpts to stay consistient
{
TMD1endpts <- reorientTMDendpts(TMD1endpts,1,tmdDomainNum,lacyAtoms)
TMD2endpts <- reorientTMDendpts(TMD2endpts,2,tmdDomainNum,lacyAtoms)
TMD3endpts <- reorientTMDendpts(TMD3endpts,3,tmdDomainNum,lacyAtoms)
TMD4endpts <- reorientTMDendpts(TMD4endpts,4,tmdDomainNum,lacyAtoms)
TMD5endpts <- reorientTMDendpts(TMD5endpts,5,tmdDomainNum,lacyAtoms)
TMD6endpts <- reorientTMDendpts(TMD6endpts,6,tmdDomainNum,lacyAtoms)
TMD7endpts <- reorientTMDendpts(TMD7endpts,7,tmdDomainNum,lacyAtoms)
TMD8endpts <- reorientTMDendpts(TMD8endpts,8,tmdDomainNum,lacyAtoms)
TMD9endpts <- reorientTMDendpts(TMD9endpts,9,tmdDomainNum,lacyAtoms)
TMD10endpts <- reorientTMDendpts(TMD10endpts,10,tmdDomainNum,lacyAtoms)
TMD11endpts <- reorientTMDendpts(TMD11endpts,11,tmdDomainNum,lacyAtoms)
TMD12endpts <- reorientTMDendpts(TMD12endpts,12,tmdDomainNum,lacyAtoms)
}

#rotate each tmd
{
lacyTMD1Rot <- rotateFromTo(lacyTMD1,lacyTMD7)
lacyTMD2Rot <- rotateFromTo(lacyTMD2,lacyTMD8)
lacyTMD3Rot <- rotateFromTo(lacyTMD3,lacyTMD9)
lacyTMD4Rot <- rotateFromTo(lacyTMD4,lacyTMD10)
lacyTMD5Rot <- rotateFromTo(lacyTMD5,lacyTMD11)
lacyTMD6Rot <- rotateFromTo(lacyTMD6,lacyTMD12)
lacyTMD7Rot <- rotateFromTo(lacyTMD7,lacyTMD1)
lacyTMD8Rot <- rotateFromTo(lacyTMD8,lacyTMD2)
lacyTMD9Rot <- rotateFromTo(lacyTMD9,lacyTMD3)
lacyTMD10Rot <- rotateFromTo(lacyTMD10,lacyTMD4)
lacyTMD11Rot <- rotateFromTo(lacyTMD11,lacyTMD5)
lacyTMD12Rot <- rotateFromTo(lacyTMD12,lacyTMD6)
}

#fit lines to the rotated data
{
TMD1Rotendpts <- fitLine(lacyTMD1Rot)
TMD2Rotendpts <- fitLine(lacyTMD2Rot)
TMD3Rotendpts <- fitLine(lacyTMD3Rot)
TMD4Rotendpts <- fitLine(lacyTMD4Rot)
TMD5Rotendpts <- fitLine(lacyTMD5Rot)
TMD6Rotendpts <- fitLine(lacyTMD6Rot)
TMD7Rotendpts <- fitLine(lacyTMD7Rot)
TMD8Rotendpts <- fitLine(lacyTMD8Rot)
TMD9Rotendpts <- fitLine(lacyTMD9Rot)
TMD10Rotendpts <- fitLine(lacyTMD10Rot)
TMD11Rotendpts <- fitLine(lacyTMD11Rot)
TMD12Rotendpts <- fitLine(lacyTMD12Rot)
}

#flip the rotated endpts if needed
{
  TMD1Rotendpts <- fixendpts(TMD1endpts,TMD1Rotendpts)
  TMD2Rotendpts <- fixendpts(TMD2endpts,TMD2Rotendpts)
  TMD3Rotendpts <- fixendpts(TMD3endpts,TMD3Rotendpts)
  TMD4Rotendpts <- fixendpts(TMD4endpts,TMD4Rotendpts)
  TMD5Rotendpts <- fixendpts(TMD5endpts,TMD5Rotendpts)
  TMD6Rotendpts <- fixendpts(TMD6endpts,TMD6Rotendpts)
  TMD7Rotendpts <- fixendpts(TMD7endpts,TMD7Rotendpts)
  TMD8Rotendpts <- fixendpts(TMD8endpts,TMD8Rotendpts)
  TMD9Rotendpts <- fixendpts(TMD9endpts,TMD9Rotendpts)
  TMD10Rotendpts <- fixendpts(TMD10endpts,TMD10Rotendpts)
  TMD11Rotendpts <- fixendpts(TMD11endpts,TMD11Rotendpts)
  TMD12Rotendpts <- fixendpts(TMD12endpts,TMD12Rotendpts)
}

#plot the rotated
{
plotTMD(lacyTMD1Rot,TMD1Rotendpts,lacyPosAtoms,T,"blue")
plotTMD(lacyTMD2Rot,TMD2Rotendpts,lacyPosAtoms,F,"darkorange")
plotTMD(lacyTMD3Rot,TMD3Rotendpts,lacyPosAtoms,F,"red")
plotTMD(lacyTMD4Rot,TMD4Rotendpts,lacyPosAtoms,F,"yellow")
plotTMD(lacyTMD5Rot,TMD5Rotendpts,lacyPosAtoms,F,"green")
plotTMD(lacyTMD6Rot,TMD6Rotendpts,lacyPosAtoms,F,"violetred")
plotTMD(lacyTMD7Rot,TMD7Rotendpts,lacyPosAtoms,F,"skyblue")
plotTMD(lacyTMD8Rot,TMD8Rotendpts,lacyPosAtoms,F,"orange")
plotTMD(lacyTMD9Rot,TMD9Rotendpts,lacyPosAtoms,F,"indianred1")
plotTMD(lacyTMD10Rot,TMD10Rotendpts,lacyPosAtoms,F,"lightyellow")
plotTMD(lacyTMD11Rot,TMD11Rotendpts,lacyPosAtoms,F,"palegreen")
plotTMD(lacyTMD12Rot,TMD12Rotendpts,lacyPosAtoms,F,"palevioletred")
}
rgl.light(specular="#111111",ambient="#777777",viewpoint.rel=FALSE,theta=90,phi=90)
rgl.snapshot(paste("LacYRotAll",".png",sep=""))

#write out the rotated tmd
{
writeTMD(lacyTMD1Rot,"TMD1Rot.txt",1,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD2Rot,"TMD2Rot.txt",2,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD3Rot,"TMD3Rot.txt",3,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD4Rot,"TMD4Rot.txt",4,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD5Rot,"TMD5Rot.txt",5,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD6Rot,"TMD6Rot.txt",6,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD7Rot,"TMD7Rot.txt",7,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD8Rot,"TMD8Rot.txt",8,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD9Rot,"TMD9Rot.txt",9,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD10Rot,"TMD10Rot.txt",10,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD11Rot,"TMD11Rot.txt",11,tmdDomainNum,lacyAtoms)
writeTMD(lacyTMD12Rot,"TMD12Rot.txt",12,tmdDomainNum,lacyAtoms)
}