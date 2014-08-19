#Various tools for a first exploratory analysis of the VMS data 



VMS <- read.csv('jittered_VMS.csv',sep=',') 
IDs <- length(unique(VMS[,2]))

#segment the data by vessel ID and -> get what is a season, essentially 
#for each i in IDs 
#season$i = data.frame([VMS[which(VMS[,2] == i)],])
################################################################3
# Legacy cleaning 
VMS_1 <- VMS[VMS$vessel_id==1,]
VMS_2 <- VMS[VMS$vessel_id==2,] 
VMS_3 <- VMS[VMS$vessel_id==3,]

# Split on vessel id

#order by date, once segmented by vessel id
VMS_1 <- VMS_1[order(VMS_1$Date_Time),]
VMS_2 <- VMS_2[order(VMS_2$Date_Time),]
VMS_3 <- VMS_3[order(VMS_3$Date_Time),]

VMS_1_stats <- stats(VMS_1$Latitude,VMS_1$Longitude,VMS_1$Date_Time,VMS_1$onland,4)
VMS_2_stats <- stats(VMS_2$Latitude,VMS_2$Longitude,VMS_2$Date_Time,VMS_2$onland,4)
VMS_3_stats <- stats(VMS_3$Latitude,VMS_3$Longitude,VMS_3$Date_Time,VMS_3$onland,4)

#Isolate the indices of the times when the onland indicator is on
onland <- which(VMS_1$onland == 1) 
# Choose the indices where the fisher is not continuously on land
jumpsonland <- onland[which(diff(onland) != 1)]
lapply(jumpsonland,printDist) 



# Extract viable speeds 
# plot(VMS_1$Avg_Speed[(VMS_1$Avg_Speed > 0.4) & (VMS_1$onland == 0) & (VMS_1$Avg_Speed < 50)])
# plot(VMS_1$Avg_Speed[(VMS_1$Avg_Speed > 0.4) & (VMS_1$onland == 0) & (VMS_1$Avg_Speed < 50)])

#using the lubridate package to extract hour from POSIX
VMS_1_offland <- VMS_1[which(VMS_1$onland == 0),]
VMS_1_offland$hour <- hour(VMS_1_offland$Date_time) 
awake <- subset(VMS_1_offland, VMS_1_offland$hour > 5)


#more cleaning
awake <- awake[(!which(is.na(awake)),]
awake <- awake[(awake$Avg_Speed < 50),]
 

speeds <- as.numeric(VMS_1$Avg_Speed) 
plot(cumsum(as.numeric(speeds[1:2800])),xlab="Measurement #",ylab="Cumulative Speed from VMS data")
plot(cumsum(VMS_1$Avg_Speed[2000:7500]),xlab="index",ylab="Cumulative Sum of speed",main="Vessel 3, overview of cumulative speed", col=as.factor(VMS_1$onland))

VMS_1_split <- split_VMS(VMS_1,split) 



### check for missing values
####################
####### NOTE: should properly be wrapped in a function with rm or not as a flag 
for (i in 2:length(VMS_1_split)) {
	df <- VMS_1_split[[i]]
	# check by column 
	for (j in 1:dim(df)[2]) {
		derp <- which(is.na(df[,j])) # check each column for missing values 
		if (length(derp) > 0) {
		cat("Missing at iteration", i, ", column: ", j, ": ", derp,"\n")  # print for debugging purposes 
		VMS_1_split[[i]] <- df[-derp,] # Remove the offending row 
		#################################
		# NOTE: It looks like the missing values are mostly headings, which could potentially be calculated from position data
		# TODO - more proper missing data 
		}
	}
}

# iterate through again and check that everything was removed 
which_isna_iterate <- function(df_list)
for (i in 2:length(VMS_1_split)) {
	df <- VMS_1_split[[i]]
	# check by column 
	for (j in 1:dim(df)[2]) {
		derp <- which(is.na(df[,j])) # check each column for missing values 
		if (length(derp) > 0) {
		cat("Missing at iteration", i, ", column: ", j, ": ", derp,"\n")  # print for debugging purposes 
		}
	}
}


# After splitting: check that each trip has enough data points, arbitrarily set to 20 
VMS_1_indices <- sapply(VMS_1_split, function(x) any(dim(x)[1] > 20) )
VMS_2_indices <- sapply(VMS_2_split, function(x) any(dim(x)[1] > 20) )
VMS_3_indices <- sapply(VMS_3_split, function(x) any(dim(x)[1] > 20) )


########################
# Compute the cubic hermite spline interpolation for each trip and write to a list of 
# such data frames 
# CRmInterpolate is modified catmull-rom from Russo, Parisi, Cataudella, 2011
# 
# resolution: # minutes between observations
resolution <- 10
interpolated_list <- list() 
for (i in 2:length(VMS_1_split)) {
	df <- VMS_1_split[[i]]
	interpolated_list[[i]] <- CRmInterpolate(df$Longitude,df$Latitude,df$Avg_Direction,df$Avg_Speed,10,df$Date_Time)
	#print(i) #print for diagnostic purposes
}

stats_list <- list() 
for (i in 2:length(interpolated_list)) {
	df <- interpolated_list[[i]]
	# if (df != NULL) {
	stats_list[[i]] <- stats_interp(df[,1],df[,2],10,4)
	# }
	#print(i) #print for diagnostic purposes
}



# first element is data frame composed of the consecutively 0 speed states 
# need to index with double square brackets 
stats_list <- lapply(interpolated_list, plystats) 



#################
#################
################# MISCELLANEOUS OLD TESTING CODE
output <- clustCombi(awake$Avg_Speed) 


night <- (notconsec_VMS1$hour > 21 | notconsec_VMS1$hour < 5)
notconsec_day_VMS1 <- subset(notconsec_VMS1,notconsec_VMS1$night == FALSE)

plot(cumsum(notconsec_day_VMS1$Avg_Speed[1:500]), col=as.factor(notconsec_day_VMS1$night)) 



test_bcpa <- CRmInterpolate(VMS_1_split[[2]]$Longitude, VMS_1_split[[2]]$Latitude, VMS_1_split[[2]]$Avg_Direction, VMS_1_split[[2]]$Avg_Speed,10,VMS_1_split[[2]]$Date_Time) 
stats_BCPA(test_bcpa[,3],test_bcpa[,4]) 
stats_BCPA_VMS(complex(re=VMS_1_split_lt[[35]]$Longitude,im=VMS_1_split_lt[[35]]$Latitude),VMS_1_split_lt[[35]]$Date_Time,30) 


 # # BCPA TESTING
HSMM_interp_raw <- CRmInterpolate(raw_data_HSMM_interp$Longitude,raw_data_HSMM_interp$Latitude,raw_data_HSMM_interp$Avg_Direction,raw_data_HSMM_interp$Avg_Speed,10,raw_data_HSMM_interp$Date_Time)
HSMM_interp_stats <- stats_BCPA_regular_notrain(complex(re=HSMM_interp_raw[,1],im=HSMM_interp_raw[,2]),HSMM_interp_raw[,3],50) 

sum_stats
HSMM <- hsmm.viterbi(sum_stats)