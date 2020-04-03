
NetCDFExtract = function(name){
library(dplyr)
library(ncdf4)



Met = nc_open(paste0(name,"_Met.nc"))
Flux = nc_open(paste0(name,"_Flux.nc"))

time = ncvar_get(Met,"time")
Tair = ncvar_get(Met,"Tair")
Precip = ncvar_get(Met,"Precip")
VPD = ncvar_get(Met,"VPD")
NEE = ncvar_get(Flux,"NEE")
GPP = ncvar_get(Flux,"GPP_DT")

ncatt_get(Met, "time")
time_from = substr(ncatt_get(Met, "time")$units, 15, 33)
time <- as.POSIXct(time, origin=time_from)

nc_close(Met)
nc_close(Flux)

# Define unit conversions
SEC_TO_30MIN = 1800
KEL_TO_CEL = -273.15

# Create data frame with the units converted
# Precip is now in mm
# Tair is now in Celsius
df <- data.frame(time, NEE, GPP, Tair+KEL_TO_CEL,Precip*SEC_TO_30MIN,VPD)
colnames(df)<- c("time","NEE","GPP","Tair","Precip","VPD")


df_day <- df %>%
      mutate(day=as.Date(time, format="%Y-%m-%d")) %>%
      group_by(day) %>%               # group by the day column
      summarise(NEE=sum(NEE),GPP=sum(GPP),Tair=mean(Tair),Precip=sum(Precip),VPD=mean(VPD))

df_month <- df_day %>%
      mutate(month=format(as.Date(df_day$day),"%Y-%m")) %>%
      group_by(month) %>%               # group by the day column
      summarise(NEE=sum(NEE),GPP=sum(GPP),Tair=mean(Tair),Precip=sum(Precip),VPD=mean(VPD))


# Create the NEE input matrix

df_day$Year = format(df_day$day,"%Y")
NEE = aggregate(df_day['NEE'],by=df_day['Year'],sum)
NEE$YearID = 1:nrow(NEE)

# Create the GPPinput matrix

df_day$Year = format(df_day$day,"%Y")
GPP = aggregate(df_day['GPP'],by=df_day['Year'],sum)
GPP$YearID = 1:nrow(GPP)

# Create our Precip input matrix
Precip = data.frame(NEE$Year,matrix(df_month$Precip,ncol = 12, byrow = TRUE))
colnames(Precip) <- c("Year","ppt1","ppt2","ppt3","ppt4","ppt5","ppt6","ppt7","ppt8","ppt9","ppt10","ppt11","ppt12")

# Create our temp input matrix
Tair = data.frame(NEE$Year,matrix(df_month$Tair,ncol = 12, byrow = TRUE))
colnames(Tair) <- c("Year","Tair1","Tair2","Tair3","Tair4","Tair5","Tair6","Tair7","Tair8","Tair9","Tair10","Tair11","Tair12")


Flux = list("Precip"=Precip,"NEE"=NEE,"GPP"=GPP,"Tair"=Tair)
}

