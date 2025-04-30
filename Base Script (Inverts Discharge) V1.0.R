#install.packages("installr")
#library(installr)
#updateR()

                                                      ##VERSION 1.0##

                                               ##Invert Discharge Analysis##

#After installing packages on your PC once, you do not need to run the packages again, but you do need to run
#the libraries every time you reopen RStudio.
#install.packages('vegan')
library(vegan)
#install.packages("dplyr")
library(dplyr)                            #####################################
                                          ##Start here after loading packages##
                                          #####################################

#Copy the file address for the exported Database spreadsheets here. Make sure they are saved as .csv files.
Q<-read.csv("C:/Users/ccheri/Documents/R/DB Exports/Discharge/BUFF_Invert_ExportDischargeDetailData.csv")  #This is the raw discharge data.

#Before we get started, let's change any -999 values to NA's in our dataframe.
Q[Q == "-999"] <- as.numeric("NA")

#First we need to add some necessary columns to calculate total discharge.We also need to convert depth from cm to m.

#We'll create a column that recognizes the previous bank distance readings so we can calculate incremental widths
#between DischargeNo (velocity) readings. The previous value for every first bank distance reading that is 0 will be NA.
Q <- Q %>% group_by(Season, LocationID) %>% 
  mutate("prev_value" = ifelse(DistanceFromBank_m!=0 & DischargeNo==1, DistanceFromBank_m, lag(DistanceFromBank_m)))

#Here we are creating a dataframe with the incremental width values by subtracting each bank distance reading by
#the previous one (again the first will usually be "0" because nothing *should* come before the first reading).
width_calc <- Q %>%
  filter(!is.na(Season)) %>%
  group_by(Season, LocationID, DischargeNo) %>%      
  summarize(Increment_width = ifelse(DistanceFromBank_m!=0 & DischargeNo==1, prev_value-0, (DistanceFromBank_m-prev_value)))

#Now we'll bind the incremental width values to our main dataframe "Q".
Q <- cbind(Q, Increment_width = width_calc$Increment_width)

#Next we'll convert the depth readings from cm to m and then rename that column.
Q$Depth_cm <- Q$Depth_cm * .01
names(Q)[names(Q)=="Depth_cm"] <- "Depth_m"

#Then we'll add a cross-sectional area column with area calculations for each cross-sectional point where a reading
#was taken.
Q <- cbind(Q, Area = Q$Increment_width * Q$Depth_m)

#Before we calculate incremental discharge we'll force all velocity readings to be positive numbers.
Q$Velocity_ms = abs(Q$Velocity_ms)

#Now we'll calculate incremental discharge values where each reading was taken and bind those to the dataframe as well.
Q <- cbind(Q, Increment_Q = Q$Area * Q$Velocity_ms)

#The following code will produce a discharge stats summary for each year and site in the data.
Q_summary <- Q %>%
  group_by(LocationID, Season) %>%      
  summarize(Discharge = sum(Increment_Q, na.rm = TRUE))

#Let's see what we got!
Q_summary

#Lastly we'll save the summary as a .csv file so that we can load it into a new R script that we'll use to make
#various plots.

write.csv(Q_summary, "C:/Users/ccheri/Documents/R/Summaries/Park Discharge Summaries/BUFF_Inverts_Qsummary.csv") #You will need to change 
                                                                                                                 #the name and address accordingly.
