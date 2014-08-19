# Dependencies: VMS_ETL_helper_fns.R
#helper wrapper for stats for use in apply/sapply/lapply 

plystats <- function(df) {
	win <- 3
	stats(df$Latitude,df$Longitude,df$Date_Time,df$onland,win)
}

stats <- function(x,y,tm,onland,win) {
	n = length(x)
	# Can be rewritten/compiled in C later for speed purposes.. 
	# iterate over number of rows in dataset
	speeds <- rep(0,n); 
	angles <- rep(0,n); 
	dists <- rep(0,n); 
	
	#for sinuosity
	slidingC <- rep(0,n); p <- rep(0,n); b <- rep(0,n); 
	lagDist <- rep(0,n); 
	sumDist <- 0; 

	vect1 <- rbind(0,0)
	for (i in (win+1):(n-1)) {
		if (onland[i] == 0) {
		vect2 <- c(x[i+1] - x[i], y[i+1] - y[i])
		dists[i] <- gcd.hf(x[i+1],y[i+1],x[i],y[i])
		#Assume time is given as POSIXct
		timediff <- difftime(tm[i],tm[i-1],units="hours")
		speeds[i] <- dists[i] / as.numeric(timediff)
		

		angles[i] <- atan(vect1[2]/vect2[1]) - atan(vect1[2] / vect1[1] )	
		if (sqrt(sum(vect2^2)) == 0) angles[i] = 0 

		slidingC[i] <- mean(cos(angles[(i-win):i]))
		p[i] <- mean(dists[(i-win):i]) 
		b[i] <- sd(dists[(i-win):i ])

		expDist <- dist(rbind(c(x[i],y[i]),c(x[i-win],y[i-win])))
		lagDist[i] <- expDist / sum(dists[(i-win):i])
		if (expDist == 0) lagDist[i] = 0

		vect1 <- vect2

		}
	}

	sinuosity <- rep(0,n)
	for (i in (win+1):n) {
		if (onland[i] == 0) {
		ratio <- (1 + slidingC[i]) / (1 - slidingC[i])
		sinuosity[i] = 2 / sqrt(p[i] * ratio + b[i]*b[i])
		}
	}
	#retList <- list("speeds" = speeds, "angles" = angles, "sinuosity" = sinuosity, "lagDist" = lagDist, "meanCos" = slidingC)
	retList <- data.frame("speeds" = speeds, "angles" = angles, "sinuosity" = sinuosity, "lagDist" = lagDist, "meanCos" = slidingC)

	return(retList)
}

#######################################
# Calculate the statistics on the interpolated data 
# Report resolution in minutes, convert to hours
stats_interp <- function(x,y,resolution,win) {
	n = length(x)
	# Can be rewritten/compiled in C later for speed purposes.. 
	# iterate over number of rows in dataset
	speeds <- rep(0,n); 
	angles <- rep(0,n); 
	dists <- rep(0,n); 
	
	#for sinuosity
	slidingC <- rep(0,n); p <- rep(0,n); b <- rep(0,n); 
	lagDist <- rep(0,n); 
	sumDist <- 0; 

	timediff <- resolution / 60


	vect1 <- rbind(0,0)
	for (i in (win+1):(n-1)) {

		vect2 <- c(x[i+1] - x[i], y[i+1] - y[i])
		dists[i] <- gcd.hf(x[i+1],y[i+1],x[i],y[i])

		####################
		# this is the section different from stats
		#Assume time is given as POSIXct
		# Calculate difference in time in hours
		speeds[i] <- dists[i] / timediff
		###############

		angles[i] <- atan(vect1[2]/vect2[1]) - atan(vect1[2] / vect1[1] )	
		if (sqrt(sum(vect2^2)) == 0) angles[i] = 0 

		slidingC[i] <- mean(cos(angles[(i-win):i]))
		p[i] <- mean(dists[(i-win):i]) 
		b[i] <- sd(dists[(i-win):i ])

		expDist <- dist(rbind(c(x[i],y[i]),c(x[i-win],y[i-win])))
		lagDist[i] <- expDist / sum(dists[(i-win):i])
		if (expDist == 0) lagDist[i] = 0

		vect1 <- vect2

		
	}

	sinuosity <- rep(0,n)
	for (i in (win+1):n) {
		ratio <- (1 + slidingC[i]) / (1 - slidingC[i])
		sinuosity[i] = 2 / sqrt(p[i] * ratio + b[i]*b[i])
		
	}
	#retList <- list("speeds" = speeds, "angles" = angles, "sinuosity" = sinuosity, "lagDist" = lagDist, "meanCos" = slidingC)
	retList <- data.frame("speeds" = speeds, "angles" = angles, "sinuosity" = sinuosity, "lagDist" = lagDist, "meanCos" = slidingC)

	return(retList)
}

#From Gurarie's dissertation and documentation
stats_BCPA_regular <- function(Z,Times,ws) {
Times <- Times / 3600 # seconds into hours 
dZ <- diff(Z) 
# orientation
Phi <- Arg(dZ) 
# turning angles
Theta <- diff(Phi) 
# step lengths
S <- Mod(dZ)
# time intervals
dT <- diff(Times)
# magnitude of linear velocity between points
V <- S / dT 
V <- V[-1]
VC <- V*cos(Theta) 
VS <- V*sin(Theta) 

T <- (Times[-1] + Times[-length(Times)]) / 2
T <- T[-1]

vc.sweep <- WindowSweep(VC, T, windowsize = ws, windowstep = 1) 
vc.output <- PartitionParameters(vc.sweep, T, windowsize=ws,windowstep = 1) 

PlotBCPA(T, VC, vc.sweep, vc.output, threshold=10)
 PathPlot(T, Z, vc.sweep, vc.output, threshold = 10)

 return( list(vc.sweep,vc.output, Phi, Theta, S, V, VC, VS ) )

}

stats_BCPA_VMS <- function(Z,Times,ws) {
dZ <- diff(Z) 
# orientation
Phi <- Arg(dZ) 
# turning angles
Theta <- diff(Phi) 
# step lengths
S <- Mod(dZ)
# time intervals
dT <- diff(Times)
dT <- as.numeric(dT) / 60 # in minutes
# magnitude of linear velocity between points
V <- S / dT
V <- V[-1]
VC <- V*cos(Theta) 
VS <- V*sin(Theta) 

T <- (as.numeric(Times[-1]) + as.numeric(Times[-length(Times)])) / 2
T <- T[-1]

vc.sweep <- WindowSweep(VC, T, windowsize = ws, windowstep = 1) 
vc.output <- PartitionParameters(vc.sweep, T, windowsize=ws,windowstep = 1) 

PlotBCPA(T, VC, vc.sweep, vc.output, threshold=10)
 PathPlot(T, Z, vc.sweep, vc.output, threshold = 10)

 return( list("VCSweep" = vc.sweep,"VCOutput" = vc.output, "Phi" = Phi, "Theta" = Theta, "S" = S, "V" = V, 
 	"VC" = VC, "VS" = VS ) )

}

stats_BCPA_regular_notrain <- function(Z,Times,ws) {
Times <- Times / 3600 # seconds into hours 
dZ <- diff(Z) 
# orientation
Phi <- Arg(dZ) 
# turning angles
Theta <- diff(Phi) 
# step lengths
S <- Mod(dZ)
# time intervals
dT <- diff(Times)
# magnitude of linear velocity between points
V <- S / dT 
V <- V[-1]
VC <- V*cos(Theta) 
VS <- V*sin(Theta) 

T <- (Times[-1] + Times[-length(Times)]) / 2
T <- T[-1]


 return( list(Phi, Theta, S, V, VC, VS ) )

}