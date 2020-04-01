source("SAMFunctionTair.R")
source("NCExtraction.R")
library(ggplot2)
library(gridExtra)

Flux = NetCDFExtract("US-Var_2001-2014_FLUXNET2015")

timeblock = SpringBlock(3)
Nlag = timeblock$Nlag
block = timeblock$block
# For the timeblock created, run the model
SAMTair(Flux$NEE[-(1:3),],Flux$Tair,Nlag,block,prior=FALSE)

load("SAM_Tair_posterior_3_6_20200401_125725.Rdata")
NPPobs = data.frame(Year=1:nrow(Flux$NEE)+2001,
                    NPP_obs = Flux$NEE[,2])
NPPmod = data.frame(Year=1:length(SAM_Tair_posterior_3_6_20200401_125725$NPPmod$mean)+2004,
                    NPP_mod = SAM_Tair_posterior_3_6_20200401_125725$NPPmod$mean,
                    NPP_modmin = SAM_Tair_posterior_3_6_20200401_125725$NPPmod$min,
                    NPP_modmax = SAM_Tair_posterior_3_6_20200401_125725$NPPmod$max)
plot1 <- ggplot(NPPobs) +
  geom_line(data=NPPobs,aes(Year,NPP_obs),color='steelblue',size=2) +
  geom_point(data=NPPobs,aes(Year,NPP_obs),color='steelblue',size=2,na.rm=TRUE) +
  geom_ribbon(data=NPPmod, aes(x=Year, ymin=NPP_modmin, ymax=NPP_modmax), fill="grey70", alpha=0.4) +
  geom_line(data=NPPmod,aes(Year,NPP_mod)) +
  theme_bw() +
  theme(axis.line=element_line(colour = "black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank()) 

# Assemble the dataframe for spring/non-spring
plotWeights = data.frame("Weights"=as.vector(t(SAM_Tair_posterior_3_6_20200401_125725$monthlyWeights$mean)),
                         "YearIntoPast"=rep(0:2,each=12),
                         "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                       "Jun", "May","Apr","Mar","Feb","Jan"),3))
# Assign the factors
plotWeights$Month = factor(plotWeights$Month, levels = c("Jan","Feb","Mar","Apr",
                                                         "May","Jun","Jul","Aug",
                                                         "Sep","Oct","Nov","Dec"))
# Create plot ready for grid.assign
plot2 <- ggplot(plotWeights,aes(x=YearIntoPast,y=Weights,fill=Month)) +
  geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
  facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y") +
  ggtitle(paste0("VAR 3 year lag, spring/non-spring TEMP (sum) (R2 = ",
                 signif(SAM_Tair_posterior_3_6_20200401_125725$R2,2), ", DIC = ",
                 signif(SAM_Tair_posterior_3_6_20200401_125725$DIC,4), ", MAE = ",
                 signif(SAM_Tair_posterior_3_6_20200401_125725$MAE,4),")"))
grid.arrange(plot2, plot1)




# 
# Load Ogle data
Tair = read.table("data/OgleTemp.csv",header=TRUE,stringsAsFactors = FALSE,sep=",")
ppt = read.table("data/dataset3.csv",header=TRUE,stringsAsFactors = FALSE)
# Load NPP data
NPP = read.table("data/dataset2.csv",header=TRUE,stringsAsFactors = FALSE)
NPP = NPP[,1:3]

SAMTair(NPP,ppt*25.4,Tair,Nlag,block,prior=FALSE)


load("SAM_Tair_posterior_3_6_20200401_133609.Rdata")
NPPobs = data.frame(Year=1:nrow(NPP)+1938,
                    NPP_obs = NPP[,2])
NPPmod = data.frame(Year=1:length(SAM_Tair_posterior_3_6_20200401_133609$NPPmod$mean)+1938,
                    NPP_mod = SAM_Tair_posterior_3_6_20200401_133609$NPPmod$mean,
                    NPP_modmin = SAM_Tair_posterior_3_6_20200401_133609$NPPmod$min,
                    NPP_modmax = SAM_Tair_posterior_3_6_20200401_133609$NPPmod$max)
plot3 <- ggplot(NPPobs) +
  geom_line(data=NPPobs,aes(Year,NPP_obs),color='steelblue',size=2) +
  geom_point(data=NPPobs,aes(Year,NPP_obs),color='steelblue',size=2,na.rm=TRUE) +
  geom_ribbon(data=NPPmod, aes(x=Year, ymin=NPP_modmin, ymax=NPP_modmax), fill="grey70", alpha=0.4) +
  geom_line(data=NPPmod,aes(Year,NPP_mod)) +
  theme_bw() +
  theme(axis.line=element_line(colour = "black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank()) 

# Assemble the dataframe for spring/non-spring
plotWeights = data.frame("Weights"=as.vector(t(SAM_Tair_posterior_3_6_20200401_133609$monthlyWeights$mean)),
                         "YearIntoPast"=rep(0:2,each=12),
                         "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                       "Jun", "May","Apr","Mar","Feb","Jan"),3))
# Assign the factors
plotWeights$Month = factor(plotWeights$Month, levels = c("Jan","Feb","Mar","Apr",
                                                         "May","Jun","Jul","Aug",
                                                         "Sep","Oct","Nov","Dec"))
# Create plot ready for grid.assign
plot4 <- ggplot(plotWeights,aes(x=YearIntoPast,y=Weights,fill=Month)) +
  geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
  facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y") +
  ggtitle(paste0("Ogle 3 year lag, spring/non-spring TEMP (sum) (R2 = ",
                 signif(SAM_Tair_posterior_3_6_20200401_133609$R2,2), ", DIC = ",
                 signif(SAM_Tair_posterior_3_6_20200401_133609$DIC,4), ", MAE = ",
                 signif(SAM_Tair_posterior_3_6_20200401_133609$MAE,4),")"))
grid.arrange(plot4, plot3)




load("SAM_Tair_posterior_3_6_20200401_144851.Rdata")
NPPobs = data.frame(Year=1:nrow(NPP)+1938,
                    NPP_obs = NPP[,2])
NPPmod = data.frame(Year=1:length(SAM_Tair_posterior_3_6_20200401_144851$NPPmod$mean)+1938,
                    NPP_mod = SAM_Tair_posterior_3_6_20200401_144851$NPPmod$mean,
                    NPP_modmin = SAM_Tair_posterior_3_6_20200401_144851$NPPmod$min,
                    NPP_modmax = SAM_Tair_posterior_3_6_20200401_144851$NPPmod$max)
plot5 <- ggplot(NPPobs) +
  geom_line(data=NPPobs,aes(Year,NPP_obs),color='steelblue',size=2) +
  geom_point(data=NPPobs,aes(Year,NPP_obs),color='steelblue',size=2,na.rm=TRUE) +
  geom_ribbon(data=NPPmod, aes(x=Year, ymin=NPP_modmin, ymax=NPP_modmax), fill="grey70", alpha=0.4) +
  geom_line(data=NPPmod,aes(Year,NPP_mod)) +
  theme_bw() +
  theme(axis.line=element_line(colour = "black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank()) 

# Assemble the dataframe for spring/non-spring
plotWeights = data.frame("Weights"=as.vector(t(SAM_Tair_posterior_3_6_20200401_144851$monthlyWeights$mean)),
                         "YearIntoPast"=rep(0:2,each=12),
                         "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                       "Jun", "May","Apr","Mar","Feb","Jan"),3))
# Assign the factors
plotWeights$Month = factor(plotWeights$Month, levels = c("Jan","Feb","Mar","Apr",
                                                         "May","Jun","Jul","Aug",
                                                         "Sep","Oct","Nov","Dec"))
# Create plot ready for grid.assign
plot6 <- ggplot(plotWeights,aes(x=YearIntoPast,y=Weights,fill=Month)) +
  geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
  facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y") +
  ggtitle(paste0("Ogle 3 year lag, spring/non-spring TEMP (sum) (R2 = ",
                 signif(SAM_Tair_posterior_3_6_20200401_144851$R2,2), ", DIC = ",
                 signif(SAM_Tair_posterior_3_6_20200401_144851$DIC,4), ", MAE = ",
                 signif(SAM_Tair_posterior_3_6_20200401_144851$MAE,4),")"))
grid.arrange(plot6, plot5, plot4, plot3,nrow=4)




# The below used Precip values in inches... tut tut, we want mm.
# However, I fully expect absolutely no difference in result

load("SAM_Tair_posterior_3_6_20200401_152932.Rdata")
NPPobs = data.frame(Year=1:nrow(NPP)+1938,
                    NPP_obs = NPP[,2])
NPPmod = data.frame(Year=1:length(SAM_Tair_posterior_3_6_20200401_152932$NPPmod$mean)+1938,
                    NPP_mod = SAM_Tair_posterior_3_6_20200401_152932$NPPmod$mean,
                    NPP_modmin = SAM_Tair_posterior_3_6_20200401_152932$NPPmod$min,
                    NPP_modmax = SAM_Tair_posterior_3_6_20200401_152932$NPPmod$max)
plot7 <- ggplot(NPPobs) +
  geom_line(data=NPPobs,aes(Year,NPP_obs),color='steelblue',size=2) +
  geom_point(data=NPPobs,aes(Year,NPP_obs),color='steelblue',size=2,na.rm=TRUE) +
  geom_ribbon(data=NPPmod, aes(x=Year, ymin=NPP_modmin, ymax=NPP_modmax), fill="grey70", alpha=0.4) +
  geom_line(data=NPPmod,aes(Year,NPP_mod)) +
  theme_bw() +
  theme(axis.line=element_line(colour = "black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank()) +
ggtitle(paste0("NPP = ",
               signif(SAM_Tair_posterior_3_6_20200401_152932$alphas$mean[1],3), " + ",
               signif(SAM_Tair_posterior_3_6_20200401_152932$alphas$mean[2],3), " * TAir + ",
               signif(SAM_Tair_posterior_3_6_20200401_152932$alphas$mean[3],3)," * PPT"))

# Assemble the dataframe for spring/non-spring
plotWeights_T = data.frame("Temp_Weights"=as.vector(t(SAM_Tair_posterior_3_6_20200401_152932$monthlyWeights_T$mean)),
                         "YearIntoPast"=rep(0:2,each=12),
                         "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                       "Jun", "May","Apr","Mar","Feb","Jan"),3))
# Assign the factors
plotWeights_T$Month = factor(plotWeights_T$Month, levels = c("Jan","Feb","Mar","Apr",
                                                         "May","Jun","Jul","Aug",
                                                         "Sep","Oct","Nov","Dec"))
# Create plot ready for grid.assign
plot8 <- ggplot(plotWeights_T,aes(x=YearIntoPast,y=Temp_Weights,fill=Month)) +
  geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
  facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y")  +
ggtitle(paste0("R2 = ",
               signif(SAM_Tair_posterior_3_6_20200401_152932$R2,2), ", DIC = ",
               signif(SAM_Tair_posterior_3_6_20200401_152932$DIC,4), ", MAE = ",
               signif(SAM_Tair_posterior_3_6_20200401_152932$MAE,4)))


# Assemble the dataframe for spring/non-spring
plotWeights_P = data.frame("PPT_Weights"=as.vector(t(SAM_Tair_posterior_3_6_20200401_152932$monthlyWeights_P$mean)),
                         "YearIntoPast"=rep(0:2,each=12),
                         "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                       "Jun", "May","Apr","Mar","Feb","Jan"),3))
# Assign the factors
plotWeights_P$Month = factor(plotWeights_P$Month, levels = c("Jan","Feb","Mar","Apr",
                                                         "May","Jun","Jul","Aug",
                                                         "Sep","Oct","Nov","Dec"))
# Create plot ready for grid.assign
plot9 <- ggplot(plotWeights_P,aes(x=YearIntoPast,y=PPT_Weights,fill=Month)) +
  geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
  facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y") +
  guides(color = "none") +
  ggtitle(paste0("Data: Ogle Nlag: 3 , Block: spring/non-spring, Input: Temp + Precip"))
grid.arrange(plot9, plot8, plot7,nrow=3)


# The model below uses precipitation in mm!


load("SAM_Tair_posterior_3_6_20200401_155434.Rdata")
NPPobs = data.frame(Year=1:nrow(NPP)+1938,
                    NPP_obs = NPP[,2])
NPPmod = data.frame(Year=1:length(SAM_Tair_posterior_3_6_20200401_155434$NPPmod$mean)+1938,
                    NPP_mod = SAM_Tair_posterior_3_6_20200401_155434$NPPmod$mean,
                    NPP_modmin = SAM_Tair_posterior_3_6_20200401_155434$NPPmod$min,
                    NPP_modmax = SAM_Tair_posterior_3_6_20200401_155434$NPPmod$max)
plot7 <- ggplot(NPPobs) +
  geom_line(data=NPPobs,aes(Year,NPP_obs),color='steelblue',size=2) +
  geom_point(data=NPPobs,aes(Year,NPP_obs),color='steelblue',size=2,na.rm=TRUE) +
  geom_ribbon(data=NPPmod, aes(x=Year, ymin=NPP_modmin, ymax=NPP_modmax), fill="grey70", alpha=0.4) +
  geom_line(data=NPPmod,aes(Year,NPP_mod)) +
  theme_bw() +
  theme(axis.line=element_line(colour = "black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank()) +
  ggtitle(paste0("NPP = ",
                 signif(SAM_Tair_posterior_3_6_20200401_155434$alphas$mean[1],3), " + ",
                 signif(SAM_Tair_posterior_3_6_20200401_155434$alphas$mean[2],3), " * TAir + ",
                 signif(SAM_Tair_posterior_3_6_20200401_155434$alphas$mean[3],3)," * PPT"))

# Assemble the dataframe for spring/non-spring
plotWeights_T = data.frame("Temp_Weights"=as.vector(t(SAM_Tair_posterior_3_6_20200401_155434$monthlyWeights_T$mean)),
                           "YearIntoPast"=rep(0:2,each=12),
                           "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                         "Jun", "May","Apr","Mar","Feb","Jan"),3))
# Assign the factors
plotWeights_T$Month = factor(plotWeights_T$Month, levels = c("Jan","Feb","Mar","Apr",
                                                             "May","Jun","Jul","Aug",
                                                             "Sep","Oct","Nov","Dec"))
# Create plot ready for grid.assign
plot8 <- ggplot(plotWeights_T,aes(x=YearIntoPast,y=Temp_Weights,fill=Month)) +
  geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
  theme(legend.position = c(0.5,-0.28)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y")  +
  ggtitle(paste0("R2 = ",
                 signif(SAM_Tair_posterior_3_6_20200401_155434$R2,2), ", DIC = ",
                 signif(SAM_Tair_posterior_3_6_20200401_155434$DIC,4), ", MAE = ",
                 signif(SAM_Tair_posterior_3_6_20200401_155434$MAE,4)))


# Assemble the dataframe for spring/non-spring
plotWeights_P = data.frame("PPT_Weights"=as.vector(t(SAM_Tair_posterior_3_6_20200401_155434$monthlyWeights_P$mean)),
                           "YearIntoPast"=rep(0:2,each=12),
                           "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                         "Jun", "May","Apr","Mar","Feb","Jan"),3))
# Assign the factors
plotWeights_P$Month = factor(plotWeights_P$Month, levels = c("Jan","Feb","Mar","Apr",
                                                             "May","Jun","Jul","Aug",
                                                             "Sep","Oct","Nov","Dec"))
# Create plot ready for grid.assign
plot9 <- ggplot(plotWeights_P,aes(x=YearIntoPast,y=PPT_Weights,fill=Month)) +
  geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
  facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y") +
  theme(legend.position = "none") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  ggtitle(paste0("Data: Ogle, Nlag: 3, Block: spring/non-spring, Input: Temp + Precip"))
grid.arrange(plot9, plot8, plot7,nrow=3)


# HMMM alphas are wildly dependent on absolute values of Tair and PPT...
