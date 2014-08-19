# Interpolation functions for interpolating between VMS times
# Dependencies: VMS_ETL_helper_fns.R
# Available interpolations: cHsInterpolate, from Hintzen 2010
# CRmInterpolate, from 2011

# I think M0, M1, speed are all vectors 
# Following Hinzen 2010 
cHsInterpolate <- function(heading,M0,M1,speed) {

H0 <- speed * sin(heading)
H1 <- speed * cos(heading)

flon <- function(t) {
	tcubed <- t^3
	tsqd <- t*t
	
	(2*tcubed - 3*tsqd + 1)*M0[1] + (tcubed - 2*tsqd + t) * H0[1] + (-2*tcubed + 3*tsqd)*M1[1] + (tcubed - 2*tsqd) * H1[1]
}
flat <- function(t) {
	tcubed <- t^3
	tsqd <- t*t
	 	
	(2*tcubed - 3*tsqd + 1)*M0[2] + (tcubed - 2*tsqd + t) * H0[2] + (-2*tcubed + 3*tsqd)*M1[2] + (tcubed - 2*tsqd) * H1[2]

}
}



# 
# Following Russo, Parisi, Cataudella, 2010; modified Catmull-Rom that doesn't require parameter fitting
# x and y are vectors of spatial coordinates at sampled times
# headings is vector of vessel heading, as reported by VMS
# resolution: minutes btn observations 
# times is the POSIXct date-time vector from the VMS data to provide a properly regular interpolation

# Returns a dataframe of interpolated tracks 
CRmInterpolate <- function(x,y,headings,speeds,resolution,times) {

#################################
# Convert speed as knots into degrees latitude/longitude via nautical miles
# One knot is one arcminute of latitude per hour 
speeds <- speeds / 60 # decimal degrees latitude 

#private helper fnctions for calculation 

#Xright = X_{i+1}
#Xleft = X_{i-1}
#X = X_i
# Calculate the estimated heading for a data point and its neighbors 
Hcr <- function(Xleft, X, Xright) {
	H_i_CR <- (Xright - Xleft) / 2

	if (norm_vec(H_i_CR) == 0) return(cbind(0,0))
	# if using Euclidean distances
	# Hest <- ( dist(rbind(Xleft,X)) +  dist(rbind(Xright,X)))  * H_i_CR  / norm_vec(H_i_CR) 
	# using the gcd_hf helper function for lat/long great circle calculations
	# Corner case for stationary movement? Interpolated points should all be stationary as well?. 
	##################### GREAT CIRCLE DISTANCE: KM
	# dist1 <- gcd.hf(Xleft[[1]],Xleft[[2]],X[[1]],X[[2]])
	# dist2 <- gcd.hf(X[[1]],X[[2]],Xright[[1]],Xright[[2]])

	##################### EUCLIDEAN DISTANCE APPROXIMATION
	dist1 <- as.numeric( sqrt( (X[1]-Xleft[1])*(X[1]-Xleft[1]) + (X[2]-Xleft[2])*(X[2]-Xleft[2]) ) )
	dist2 <- as.numeric( sqrt( (Xright[1]-X[1])*(Xright[1]-X[1]) + (Xright[2]-X[2])*(Xright[2]-X[2]) ) )

	Hest <- (dist1 + dist2) * H_i_CR / norm_vec(H_i_CR)

	return((Hest))
}
###########################################################################################################
# These are the calculations for the cubic Hermite spline; t is a time on the interpolated
# interval [0,1]; M0 and M1 are location vectors of first and second point, H0 and H1 are 
# heading estimates, here obtained from the modified Catmull Rom algorithm 
flon <- function(t, M0, M1, H0, H1) {
	tcubed <- t*t*t #don't recompute every time
	tsqd <- t*t
	(2*tcubed - 3*tsqd + 1)*M0[1] + (tcubed - 2*tsqd + t) * H0[1] + (-2*tcubed + 3*tsqd)*M1[1] + (tcubed - 2*tsqd) * H1[1]
}
flat <- function(t, M0, M1, H0, H1) {
	tcubed <- t*t*t #don't recompute every time
	tsqd <- t*t 	

	(2*tcubed - 3*tsqd + 1)*M0[2] + (tcubed - 2*tsqd + t) * H0[2] + (-2*tcubed + 3*tsqd)*M1[2] + (tcubed - 2*tsqd) * H1[2]
}

# need to calculate H_est from 1 to N-1, return in a list
# somehow transform x and y lists into lists of 2-vectors 
vect_xy <- data.frame(x,y)
n <- length(x) #number of data points
H_est = as.data.frame(matrix( ncol=2, nrow=(n - 1) ))
H_p = as.data.frame(matrix( ncol=2, nrow=(n - 1) ))

X_left <- vect_xy[1,]
X_curr <- vect_xy[2,]
# Calculate estimated headings/tangents at each time point
for (i in 2:( n - 1 )) {
	X_right <- vect_xy[i + 1,]
	H_est[i,] <- Hcr(X_left, X_curr, X_right) 
	X_left <- X_curr
	X_curr <- X_right 

	H_p[i,] <- cbind(speeds[i]*cos(deg2rad(headings[i])), speeds[i]*sin(deg2rad(headings[i])))
}
	H_p <- H_p[2:(n-1),]
	H_est <- H_est[2:(n-1),] # drop the first data point for data structure alignment
### QUESTION: What does it mean for a vessel heading to have a magnitude? 


# Question: do i need to convert all of these into unit vectors? 
# headings_unit <- sapply(headings, FUN = heading2unitvec)
# # need to transpose headings_unit 
# headings_unit <- t(headings_unit) #headings_unit is a n by 2 data frame, first column x, second column y 

H_drift_i <- H_est - H_p
H_drift_med <- cbind(median(H_drift_i[,1]),median(H_drift_i[,2]))

H_crm <- H_drift_med + H_est 


#####################################################################
# return interpolated tracks up to specified resolution 
# call flon, flat for linspace(0,1,resolution) at each time 
# using time of first observation
# need to set which time interval of observation you need
# This method doesn't include the OG data by default? 


# write / initialize a new dataframe for results 
# convert interval to seconds
 resolution <- resolution * 60 
# Use it from the perspective of first heading estimated 
duration <- seconds(as.duration(as.interval( times[2], tail(times,2)[1] )))
n_obs <- as.numeric(floor( duration / resolution ))
# TIME UNITS ???????

interp_result = as.data.frame(matrix(ncol=3, nrow=(n_obs) ) )
# index through data.frame
j <- 1
k <- 1 #index rows of interpolation result

#NOTE: design better-aligned data structures next time!!!!!
# Interpolate at interval "interval" from 0 to end of time 

## Some throat-clearing initializations 
subindex <- 0
leftover <- 0 #"leftover" time for the boundary cases between intervals 
# Need to truncate vect_xy, original observations, to account for CRm undefined at first and last value
vect_xy <- vect_xy[2:(n-1),]

######## iterate through from first heading calculations to one-before-last 
TIME <- 0 
for (i in 1:(n - 3) ) {
	oneStep <- as.numeric(seconds(as.duration(as.interval(times[i], times[i+1]))))
	subindex <- leftover
	if (is.na(subindex)) cat("subindex NA")
	if (is.na(oneStep)) {
		cat ("oneStep NA\n")
		cat(time[i],"\n")
		cat(time[i+1],"\n")
	}
	while (subindex < oneStep) {
		interp_result[k,1] <- flon( (subindex / oneStep), vect_xy[i,], vect_xy[i + 1,], H_crm[i,], H_crm[i + 1,] )
		interp_result[k,2] <- flat( (subindex / oneStep), vect_xy[i,], vect_xy[i + 1,], H_crm[i,], H_crm[i + 1,] )
		interp_result[k,3] <- TIME
		k <- k + 1
		j <- j + 1
		subindex <- subindex + resolution 
		TIME <- TIME + resolution
	} 

	leftover <- subindex - oneStep
} 
return(interp_result)


#function end
}





