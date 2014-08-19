# File of private helper functions used in various other calculations for
# cleaning VMS data 

# gcd.hf: great-circle distance calculation
# printDist: diagnostic method for checking distances between VMS data points

###########################################################
# Various useful one-liners

heading2unitvec <- function(x) cbind(cos(deg2rad(x)),sin(deg2rad(x)))
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
deg2rad <- function(deg) return(deg*pi/180)
norm_vec <- function(x) sqrt(sum(x^2)) 


##################################
# SPLITTING THE DATA FRAME INTO SEPARATE TRIPS 
# Return the dataframe of VMS speeds not consecutively the same for n times


split_VMS <- function(df,n) {
# Use this on cumulative sum of speed
# FIRST: clean the NA values 
df$Avg_Speed[which(is.na(df$Avg_Speed))] <- 0 
#
df$Cum_Speed <- cumsum(df$Avg_Speed) 
df$trips <- 0 
notconsec_VMS1_Ind <- not_consec(df$Cum_Speed, n)
# find breakpoints
notconsec_breaks <- notconsec_VMS1_Ind$Ind[which(diff(notconsec_VMS1_Ind$Ind) != 1)]

# need to split on these repeats  
# doing it the brute force way 

# if difference between diffed vector is > 1 
break_ind <- diff(notconsec_VMS1_Ind$Ind) 
for (i in 1:(length(notconsec_VMS1_Ind$Ind)-1)) {
# write i to the interval vect$ind[breaks[i]:breaks[i+1]] 
# else write 0 
  if (break_ind[i] != 1) {
    df$trips[notconsec_VMS1_Ind$Ind[i]:notconsec_VMS1_Ind$Ind[i+1]] <- i 
  }
  else {
    df$trips[notconsec_VMS1_Ind$Ind[i]] <- 0 
  }
}

# split on this matrix
df_split <- split(df, df$trips) 
return(df_split)

}


# Great Circle Distance Calc in R, using Haversine for small-scale computational stability. 
# http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
# Calculates the geodesic distance between two points specified by #DEGREE# latitude/longitude using the
# Haversine formula (hf)
# Accepts measurements in degrees, converts them to radians
gcd.hf <- function(long1, lat1, long2, lat2) {
  deg2rad <- function(deg) return(deg*pi/180)
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1) 
  long2 <- deg2rad(long2) 
  lat2 <- deg2rad(lat2) 

  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

# Helper function to remove values repeated consecutively N times. 
# Returns a vector of the indices to keep 
not_consec <- function(X, N)
{
 .rle <- rle(sort(X))
 res <- .rle$values[.rle$lengths >= N]
 res <- res[res > 0]
 inds <- X %in% res
 X[inds] <- NA 
 list(X = X, Ind = which(inds)) 
}



#Check the distances jumped by each 
#Helper function that doesn't work for bounday cases (first and last element) 
printDist <- function(ind) {
	dist <- gcd.hf(VMS_1$Longitude[ind-1],VMS_1$Latitude[ind-1],VMS_1$Longitude[ind],VMS_1$Latitude[ind])
	print(dist)
	diff <- VMS_1$Date_Time[ind] - VMS_1$Date_Time[ind-1]
	if (diff < 1) print(diff)
}


#############
# Iterate through a list of data frames and print out where 
# missing values are 
which_isna_iterate <- function(df_list)
for (i in 1:length(df_list)) {
  df <- df_list[[i]]
  # check by column 
  for (j in 1:dim(df)[2]) {
    derp <- which(is.na(df[,j])) # check each column for missing values 
    if (length(derp) > 0) {
    cat("Missing at frame #", i, ", column: ", j, ": ", derp,"\n")  # print for debugging purposes 
    }
  }
}


#############
# Iterate through a list of data frames and remove missing values
which_isna_remove <- function(df_list) {
for (i in 1:length(df_list)) {
  df <- df_list[[i]]
  # check by column 
  for (j in 1:dim(df)[2]) {
    derp <- which(is.na(df[,j])) # check each column for missing values 
    if (length(derp) > 0) {
    cat("Missing at frame #", i, ", column: ", j, ": ", derp,"\n")  # print for debugging purposes 
    df_list[[i]] <- df[-derp,] # Remove the offending row 

    }
  }
}
return(df_list)
}


##########################
# Helper function to write the interpolation for a list of data frames that need interpolation
# Resolution: # minutes between observations
interpolate_list <- function(df_list, resolution) {
interpolated_list <- list() 
for (i in 1:length(df_list)) {
  df <- df_list[[i]]
  interpolated_list[[i]] <- CRmInterpolate(df$Longitude,df$Latitude,df$Avg_Direction,df$Avg_Speed,10,df$Date_Time)
  #print(i) #print for diagnostic purposes
}
return(interpolated_list)

}