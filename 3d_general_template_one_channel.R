#----------- 3D DOSE RESPONSE PROGRAM ----------- 

#  THE PURPOSE OF THIS PROGRAM IS TO TAKE A .CSV FILE AND CREATE A 3D WATERFALL PLOT

#-- The following "install" should be ran if it is the first time
#-- you run this program because most likely you will not have these
#-- libraries installed


# install.packages("rgl")
# install.packages("RColorBrewer")




#----------- LOG10 ----------- 
#- titration concentrations from the qHTS
#- This should be your titration points
conc <- c(
  -9.231,
  -9.11,
  -8.633,
  -8.532,
  -8.155,
  -8.131,
  -7.833,
  -7.678,
  -7.449,
  -7.201,
  -7.134,
  -6.867,
  -6.724,
  -6.435,
  -6.247,
  -6.235,
  -5.77,
  -5.736,
  -5.603,
  -5.293,
  -5.037,
  -4.971,
  -4.816,
  -4.338
  
)

titrations <- length(conc)

lowerBound <- min(conc)    #- Modify as needed according to your titration range
upperBound <- max(conc)     #- Modify as needed according to titration range


#- Color of points in graph

col_1 <- c("royalblue3")   # poitns in titrations of importance
col_2 <- c("gray")   # points / noise - non-actives
col_3 <- c("darkgreen")   # curve color

## all R colors can be found at http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf


#----------- S T O P !!!!!!
# ====================== ONLY CHANGE THINGS BEYOND THIS POINTS IF YOU KNOW WHAT YOU ARE DOING ========== 
#----------- OTHERWISE YOU MAY BREAK THE PROGRAM

library(rgl)
library(RColorBrewer)


# Choose the file
ifile <- file.choose()

#- set the starting value for S0, baseline at 0 or 
#- basline at 100 (for loss of signal)
baseline <- 0

#- HOW TRANSPARENT DO YOU WANT THE CURVES? 0-1

alpha_1 <- 1


useFilter <- TRUE
responseUpperFilter <- 20.0

pal <- as.list(brewer.pal(8, "Set1"))
names(pal) <- c("NHH", "NHL", "NLH", "NLL", "PHH", "PHL", "PLH", "PLL")



cdata <- read.csv(ifile, header=TRUE, na.string="null")

dataCols <- c("Protocol.Name", "Sample.Data.Type")
for(i in 1:titrations) {
  dataCols <- c(dataCols, paste("Data",(i-1), sep=""));
}

cdata <- cdata[,dataCols]
names(cdata)<-c(dataCols)

mainMatrix <- data.frame(protocolname=character(), readout=character(), x=double(), y=double(), z=double())

prevCmpd<-1
currCmpd<-2

for (i in 1:nrow(cdata)) {
  
  for(k in 1:titrations) {
    if( !is.na(cdata[i,(paste("Data",(k-1),sep=""))])) {
      if(useFilter == FALSE) {
        mainMatrix <- rbind(mainMatrix, data.frame(protocolname=cdata$Protocol.Name[i], readout=cdata$Sample.Data.Type[i], x=conc[k],  z=i, y=(cdata[i,(paste("Data",(k-1),sep=""))])))
      } else if ( as.numeric((cdata[i,(paste("Data",(k-1),sep=""))])) < responseUpperFilter) {
        mainMatrix <- rbind(mainMatrix, data.frame(protocolname=cdata$Protocol.Name[i], readout=cdata$Sample.Data.Type[i], x=conc[k],  z=i, y=(cdata[i,(paste("Data",(k-1),sep=""))])))
      }
    } 
  }
}



#----------- CORRECTING DATA TYPE -----------

for(i in mainMatrix)
{
  mainMatrix$protocolname <- as.character(mainMatrix$protocolname)
  mainMatrix$readout <- as.character(mainMatrix$readout)
  mainMatrix$x <- as.double(mainMatrix$x)
  mainMatrix$y <- as.double(mainMatrix$y)
  mainMatrix$z <- as.double(mainMatrix$z)
}


#----------- OBTAINING THE POINTS TO BE PLOTTED ----------- 


if(useFilter == FALSE) {
  waterfall_POINTS_data <- mainMatrix
  
} else {
  waterfall_POINTS_data <- mainMatrix
  
}

waterfall_POINTS_data_1 <- data.frame(x=double(), y=double(), z=double())
waterfall_POINTS_data_2 <- data.frame(x=double(), y=double(), z=double())

## adding size to the table for graphing points later

for(i in 1:length(waterfall_POINTS_data$readout)){
  if( waterfall_POINTS_data$readout[i] == 'active') {
    
    waterfall_POINTS_data_1 <- rbind(waterfall_POINTS_data_1, data.frame(x=waterfall_POINTS_data$x[i],
                                                                         y=waterfall_POINTS_data$y[i],
                                                                         z=waterfall_POINTS_data$z[i]))
  } else {
    waterfall_POINTS_data_2 <- rbind(waterfall_POINTS_data_2, data.frame(x=waterfall_POINTS_data$x[i],
                                                                         y=waterfall_POINTS_data$y[i],
                                                                         z=waterfall_POINTS_data$z[i]))
    
    
  }
  
}


#----------- RUN THE FIT RoUTINE ----------- 


cdata <- read.csv(ifile, header=TRUE, na.string="null")

cdata <- cdata[,c("Fit_Output","Protocol.Name", "CC.v2", "Sample.Data.Type", "Log.AC50..M.", "Hill.Coef", "Inf.Activity", "Zero.Activity")]
names(cdata)<-c("Fit_Output", "protocolname","klass","readout", "LAC50", "HILL", "INF", "ZERO")

#cdata$LAC50 <- as.numeric(as.character(cdata$LAC50))
#cdata$HILL <- as.numeric(as.character(cdata$HILL))


cdata$LAC50[is.na(cdata$LAC50)] <- log10(10)
cdata$HILL[ is.na(cdata$HILL)] <- 1

interleave <- function(x) {
  unlist(lapply(1:(length(x)-1), function(i) c(x[i], x[i+1])))
}

f <- function(params, concs, interleave=TRUE) {
  xx <- seq(min(concs), max(concs), length=100)
  yy <- with(params, ZERO + (INF-ZERO)/(1 + 10^( (LAC50-xx)*HILL) ))
  return(data.frame(x=xx, y=yy))
}


#f <- function(params, concs, interleave=TRUE) {
#  xx <- seq(min(concs)-0.05, max(concs)+0.05, length=100)
#  yy <- with(params, ZERO + (INF-ZERO)/(1 + 10^( (LAC50-xx)*HILL) ))
#  return(data.frame(x=xx, y=yy))
#}


#f0 <- function(params, concs, interleave=TRUE) {
#  xx <- seq(min(concs)*1.1, max(concs)*0.9, length=100)
#print(paste("*********",params$LAC50,"*******",params$HILL))

#  yy <- with(params, ZERO + (INF-ZERO)/(1 + 10^((LAC50)-xx)*HILL) )-100)
#  return(data.frame(x=xx, y=yy))
#}


mainMatrix <- data.frame(protocolname=character(),ccv2=character(),datatype=character(),cpd=integer(), x=double(),y=double(),z=double())

#mainMatrix <- data.frame(x=double(),y=double(),z=double())




prevCmpd<-1
currCmpd<-2

rowIndex = 0;
for (i in 1:nrow(cdata)) {
  
  if(cdata[i,"Fit_Output"]==1) {
    rowIndex = rowIndex+1
    if(baseline == 0) {
      d1 <- data.frame(f(cdata[i,], c(lowerBound, upperBound)),z=i)
    } else {
      d1 <- data.frame(f(cdata[i,], c(lowerBound, upperBound)),z=i)
    }
    
    #add multiple rows
    mainMatrix <- rbind(mainMatrix, data.frame(protocolname=cdata$protocolname[i], ccv2=cdata$klass[i], readout=cdata$readout[i], cpd=d1[,3], x=d1[,1], z=i, y=d1[,2]))
    
    #add the final trailing row
    mainMatrix<-rbind(mainMatrix, data.frame(protocolname=cdata$protocolname[i], ccv2='', readout=cdata$readout[i], cpd=prevCmpd, x='', z='', y=1))
    
  }
}



#----------- CORRECTING DATA TYPE -----------

mainMatrix$protocolname <- as.character(mainMatrix$protocolname)
mainMatrix$readout <- as.character(mainMatrix$readout)
mainMatrix$x <- as.double(mainMatrix$x)
mainMatrix$y <- as.double(mainMatrix$y)
mainMatrix$z <- as.double(mainMatrix$z)
mainMatrix$ccv2 <- as.character(mainMatrix$ccv2)
mainMatrix$cpd <- as.integer(mainMatrix$cpd)



#----------- OBTAINING THE LINES TO BE PLOTTED -----------

waterfall_LINES_data <- mainMatrix




# ====================== 3D GRAPHING ============================




#----------- CHANGING POP-UP WINDOW PARAMETERS -----------

# Changing parameters defaults from open3d() using par3d() for what I consider
# an initial good view of graph and window size
# The open3() window can be manually enlarged and the graph can be rotated 
# using user's mouse

## CHANGE AS NEEDED
# graph view coordinates
newMatrix = t(matrix(c(-0.6416085, 0.01792431, -0.7668231,  0, -0.1346225, 0.98157728,  0.1355843,  0,
                       0.7551258, 0.19022357, -0.6273754,  0, 0.0000000, 0.00000000,  0.0000000,  1), 
                     ncol = 4, nrow =4))

# Matrix was transposed since f(matrix) always reads y1,y2,y3,y4 
# col first then rows 

## cHANGE AS NEEDED
# Pop-up window size
newWindowRect = c(138, 161, 886, 760)



#----------- SCALING THE 3D WINDOW ----------- 

## CAREFULLY CHANGE AS NEEDED
#scale=c(concentration, %, #samples),  # 
open3d(userMatrix = newMatrix, windowRect = newWindowRect)

#----------- PLOTTING POINTS IN 3D GRAPH ----------- 

## DO NOT CHANGE
for(i in waterfall_POINTS_data_1$x)
{
  points3d(x=waterfall_POINTS_data_1$x[i], y = waterfall_POINTS_data_1$y[i],
           z=waterfall_POINTS_data_1$z[i], col = col_1,
           size = 3)
  break
}

for(i in waterfall_POINTS_data_2$x)
{
  points3d(x=waterfall_POINTS_data_2$x[i],y=waterfall_POINTS_data_2$y[i],
           z=waterfall_POINTS_data_2$z[i], col = col_2,
           size = 1)
  break
}

#----------- SELECTING VALUES FOR LINES TO BE PLOTTED ----------- 

# DO NOT CHANGE
wfl <- waterfall_LINES_data[,c("x","y","z")]

#----------- PLOTTING LINES IN 3D GRAPH ----------- 

## CHANGE AS NEEDED
for(i in wfl$x)
{ 
  lines3d(wfl[i], col = col_3, alpha = alpha_1)
  break
}


#----------- CREATE BOX AROUND EDGES FOR THE GRAPH ----------- 

axes3d(expand = 1.03, box = FALSE, xunit = 'pretty', yunit = "pretty", zunit = 'pretty')

# Adding Grid lines

grid3d(c("x","y","z+"))

aspect3d(1,1,3)




#----------- EXPORTING THE IMAGE FILES ----------- 

## the following code is put as a comment to allow the user to position the interactive 
## graph however they please
## NAME THE FILE HOWEVER YOU LIKE WITHIN THE "" MARKS BELOW

# rgl.snapshot(filename = ".png")
# rgl.postscript(".svg", fmt="svg")


