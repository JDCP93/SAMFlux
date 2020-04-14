# Source required packages
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Clear workspace
rm(list=ls())

# Source the function to extract the FLUXNET data
source("NCExtraction.R")

# Define conversion parameter
InchTomm = 25.4

#*******************************************************************************
# Load the data
#*******************************************************************************

# Above ground biomass in g/m2/yr
BDK_Met = read.csv("BDK_MET_Proc.csv")
BDK_NPP = read.csv("BDK_NPP_Proc.csv")

TWM_Met = read.csv("TWM_MET_Proc.csv")
TWM_NPP = read.csv("TWM_NPP_Proc.csv")

KRS_Met = read.csv("KRS_MET_Proc.csv")
KRS_NPP = read.csv("KRS_NPP_Proc.csv")

DZH_Met = read.csv("DZH_MET_Proc.csv")
DZH_NPP = read.csv("DZH_NPP_Proc.csv")

# GPP in g/m2/yr
Var = NetCDFExtract("US-Var_2001-2014_FLUXNET2015")
Wkg = NetCDFExtract("US-Wkg_2005-2014_FLUXNET2015")
Ha1 = NetCDFExtract("US-Ha1_1992-2012_FLUXNET2015")

# NPP in g/m2/yr
FC_Met = read.table("data/dataset3.csv",header=TRUE,stringsAsFactors = FALSE)
FC_NPP = read.table("data/dataset2.csv",header=TRUE,stringsAsFactors = FALSE)
FC_NPP = FC_NPP[,1:3]

#*******************************************************************************
# Calculations
#*******************************************************************************


# Calculate yearly MAP and mean GPP/NPP/Biomass

BDK_MAP = mean(rowSums(BDK_Met[,2:13]))
BDK_MAGPP = mean(BDK_NPP[,2],na.rm=TRUE)

TWM_MAP = mean(rowSums(TWM_Met[,2:13]))
TWM_MAGPP = mean(TWM_NPP[,2],na.rm=TRUE)

KRS_MAP = mean(rowSums(KRS_Met[,2:13]))
KRS_MAGPP = mean(KRS_NPP[,2],na.rm=TRUE)

DZH_MAP = mean(rowSums(DZH_Met[,2:13]))
DZH_MAGPP = mean(DZH_NPP[,2],na.rm=TRUE)

# Calculated VAR is very similar to that given by ORNL website
Var_MAP = mean(rowSums(Var$Precip[,2:13]))
Var_MAGPP = mean(Var$GPP[,2],na.rm=TRUE)

# NOTE THAT ORNL WEBSITE GIVES MAP OF WKG AS 407 NOT 292 AS CALCULATED!!
Wkg_MAP = mean(rowSums(Wkg$Precip[,2:13]))
Wkg_MAGPP = mean(Wkg$GPP[,2],na.rm=TRUE)

# Note MAP was in inches so convert
FC_MAP = mean(rowSums(FC_Met[,2:13]))*InchTomm
FC_MAGPP = mean(FC_NPP[,2],na.rm=TRUE)


# Put into dataframe
MAP = c(BDK_MAP,TWM_MAP,KRS_MAP,DZH_MAP,Var_MAP,Wkg_MAP,FC_MAP)
MAGPP = c(BDK_MAGPP,TWM_MAGPP,KRS_MAGPP,DZH_MAGPP,Var_MAGPP,Wkg_MAGPP,FC_MAGPP)
MAT = c(12.8,18.7,5.2,5,15.8,15.6,8.8)
Site = c("BDK","TWM","KRS","DZH","VAR","WKG","FCO")
Method = c("AGB","AGB","AGB","AGB","GPP","GPP","ANPP")
Lat = c(35.68,-24.90,51.67,49.33,38.41,31.74,40.82)
Lon = c(62,28.35,36.5,46.78,-120.95,-109.94,-104.75)

df = data.frame(Site,MAP,MAGPP,MAT,Method,Lat,Lon)

#*******************************************************************************
# Plot the data
#*******************************************************************************


# Plot GPP vs Precip
plot1 <- ggplot(df, aes(x=MAP,y=MAT,fill=Method)) + geom_point(aes(color=MAGPP),size=log(MAGPP)/2) +
  scale_color_gradient(low="blue",high="red") +
  geom_label_repel(label=Site,hjust=0, vjust=0) +
  xlab("Mean Annual Precipitation (mm)") +
  ylab("Mean Annual Temperature (degC)") +
  coord_cartesian(xlim=c(0,700),ylim=c(0,20),expand=FALSE) +
  ggtitle("Precipitation vs Productivity for Grassland Datasets")
        
# Plot location of sites
theme_set(theme_bw())
world <- ne_countries(scale = "medium", returnclass = "sf")

plot2 <- ggplot(data = world) +
  geom_sf() +
  geom_point(data=df,aes(x=Lon,y=Lat),color='red') +
  geom_label_repel(data=df,aes(x=Lon,y=Lat,label=Site)) +
  ggtitle("Location of Sites")

# Put them together
grid.arrange(plot2,plot1)

