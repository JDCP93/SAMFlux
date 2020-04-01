source("SAMFunction.R")
source("NCExtraction.R")
library(ggplot2)
library(gridExtra)

Flux = NetCDFExtract("US-Var_2001-2014_FLUXNET2015")



# timeblock = SpringBlock(3)
# Nlag = timeblock$Nlag
# block = timeblock$block
# # For the timeblock created, run the model
# SAM(Flux$NEE[-(1:3),],Flux$Precip,Nlag,block,prior=FALSE)

# VAR

load("SAM_posterior_3_6_20200331_093002.Rdata")
NPPobs = data.frame(Year=1:nrow(Flux$NEE)+2001,
                    NPP_obs = Flux$NEE[,2])
NPPmod = data.frame(Year=1:length(SAM_posterior_3_6_20200331_093002$NPPmod$mean)+2004,
                    NPP_mod = SAM_posterior_3_6_20200331_093002$NPPmod$mean,
                    NPP_modmin = SAM_posterior_3_6_20200331_093002$NPPmod$min,
                    NPP_modmax = SAM_posterior_3_6_20200331_093002$NPPmod$max)
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
plotWeights = data.frame("Weights"=as.vector(t(SAM_posterior_3_6_20200331_093002$monthlyWeights$mean)),
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
  ggtitle(paste0("VAR 3 year lag, spring/non-spring (R2 = ",
                 signif(SAM_posterior_3_6_20200331_093002$R2,2), ", DIC = ",
                 signif(SAM_posterior_3_6_20200331_093002$DIC,4), ", MAE = ",
                 signif(SAM_posterior_3_6_20200331_093002$MAE,4),")"))
grid.arrange(plot2, plot1)

#******************************************************************************

# block = matrix(1:36,nrow=3,byrow=TRUE)
# Nlag = 3
# # For the timeblock created, run the model
# SAM(Flux$NEE[-(1:3),],Flux$Precip,Nlag,block,prior=FALSE)
# 


load("SAM_posterior_3_36_20200331_094006.Rdata")
NPPobs = data.frame(Year=1:nrow(Flux$NEE)+2001,
                    NPP_obs = Flux$NEE[,2])
NPPmod = data.frame(Year=1:length(SAM_posterior_3_36_20200331_094006$NPPmod$mean)+2004,
                    NPP_mod = SAM_posterior_3_36_20200331_094006$NPPmod$mean,
                    NPP_modmin = SAM_posterior_3_36_20200331_094006$NPPmod$min,
                    NPP_modmax = SAM_posterior_3_36_20200331_094006$NPPmod$max)
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
plotWeights = data.frame("Weights"=as.vector(t(SAM_posterior_3_36_20200331_094006$monthlyWeights$mean)),
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
  ggtitle(paste0(" VAR 3 year lag, 36 blocks (R2 = ",
                 signif(SAM_posterior_3_36_20200331_094006$R2,2), ", DIC = ",
                 signif(SAM_posterior_3_36_20200331_094006$DIC,4), ", MAE = ",
                 signif(SAM_posterior_3_36_20200331_094006$MAE,4),")"))
grid.arrange(plot4, plot3)

#******************************************************************************


grid.arrange(plot2,plot1,plot4,plot3,nrow=4)


#******************************************************************************

# VAR


 # block = matrix(1:24,nrow=2,byrow=TRUE)
 # Nlag = 2
 # # For the timeblock created, run the model
 # SAM(Flux$NEE[-(1:2),],Flux$Precip,Nlag,block,prior=FALSE)

 
 load("SAM_posterior_2_24_20200331_100726.Rdata")
 NPPobs = data.frame(Year=1:nrow(Flux$NEE)+2001,
                     NPP_obs = Flux$NEE[,2])
 NPPmod = data.frame(Year=1:length(SAM_posterior_2_24_20200331_100726$NPPmod$mean)+2003,
                     NPP_mod = SAM_posterior_2_24_20200331_100726$NPPmod$mean,
                     NPP_modmin = SAM_posterior_2_24_20200331_100726$NPPmod$min,
                     NPP_modmax = SAM_posterior_2_24_20200331_100726$NPPmod$max)
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
 plotWeights = data.frame("Weights"=as.vector(t(SAM_posterior_2_24_20200331_100726$monthlyWeights$mean)),
                          "YearIntoPast"=rep(0:1,each=12),
                          "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                        "Jun", "May","Apr","Mar","Feb","Jan"),2))
 # Assign the factors
 plotWeights$Month = factor(plotWeights$Month, levels = c("Jan","Feb","Mar","Apr",
                                                          "May","Jun","Jul","Aug",
                                                          "Sep","Oct","Nov","Dec"))
 # Create plot ready for grid.assign
 plot6 <- ggplot(plotWeights,aes(x=YearIntoPast,y=Weights,fill=Month)) +
   geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
   facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y") +
   ggtitle(paste0(" VAR 2 year lag, 24 blocks (R2 = ",
                  signif(SAM_posterior_2_24_20200331_100726$R2,2), ", DIC = ",
                  signif(SAM_posterior_2_24_20200331_100726$DIC,4), ", MAE = ",
                  signif(SAM_posterior_2_24_20200331_100726$MAE,4),")"))
 grid.arrange(plot6, plot5)

 
 #******************************************************************************
  # VAR 
 
 
  # block = matrix(rep(c(1,2,3,4,5,6,7,8),each=3),nrow=2,byrow=TRUE)
  # Nlag = 2
  # # For the timeblock created, run the model
  # SAM(Flux$NEE[-(1:2),],Flux$Precip,Nlag,block,prior=FALSE)
  # 
 
 load("SAM_posterior_2_8_20200331_102101.Rdata")
 NPPobs = data.frame(Year=1:nrow(Flux$NEE)+2001,
                     NPP_obs = Flux$NEE[,2])
 NPPmod = data.frame(Year=1:length(SAM_posterior_2_8_20200331_102101$NPPmod$mean)+2003,
                     NPP_mod = SAM_posterior_2_8_20200331_102101$NPPmod$mean,
                     NPP_modmin = SAM_posterior_2_8_20200331_102101$NPPmod$min,
                     NPP_modmax = SAM_posterior_2_8_20200331_102101$NPPmod$max)
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
         panel.background=element_blank()) 
 
 # Assemble the dataframe for spring/non-spring
 plotWeights = data.frame("Weights"=as.vector(t(SAM_posterior_2_8_20200331_102101$monthlyWeights$mean)),
                          "YearIntoPast"=rep(0:1,each=12),
                          "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                        "Jun", "May","Apr","Mar","Feb","Jan"),2))
 # Assign the factors
 plotWeights$Month = factor(plotWeights$Month, levels = c("Jan","Feb","Mar","Apr",
                                                          "May","Jun","Jul","Aug",
                                                          "Sep","Oct","Nov","Dec"))
 # Create plot ready for grid.assign
 plot8 <- ggplot(plotWeights,aes(x=YearIntoPast,y=Weights,fill=Month)) +
   geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
   facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y") +
   ggtitle(paste0(" VAR 2 year lag, 8 blocks (R2 = ",
                  signif(SAM_posterior_2_8_20200331_102101$R2,2), ", DIC = ",
                  signif(SAM_posterior_2_8_20200331_102101$DIC,4), ", MAE = ",
                  signif(SAM_posterior_2_8_20200331_102101$MAE,4),")"))
 grid.arrange(plot8, plot7)
 
 
 grid.arrange(plot8,plot7,plot6,plot5,nrow=4)
 #******************************************************************************
 
 # Ogle with the original model
 
 load("~/PHD/SAMCode-master/SAM_posterior_3_6_20200331_104434.Rdata")
 NPPobs = data.frame(Year=1:nrow(NPP)+1938,
                     NPP_obs = NPP[,2])
 NPPmod = data.frame(Year=1:length(SAM_posterior_3_6_20200331_104434$NPPmod$mean)+1938,
                     NPP_mod = SAM_posterior_3_6_20200331_104434$NPPmod$mean,
                     NPP_modmin = SAM_posterior_3_6_20200331_104434$NPPmod$min,
                     NPP_modmax = SAM_posterior_3_6_20200331_104434$NPPmod$max)
 plot9 <- ggplot(NPPobs) +
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
 plotWeights = data.frame("Weights"=as.vector(t(SAM_posterior_3_6_20200331_104434$monthlyWeights$mean)),
                          "YearIntoPast"=rep(0:2,each=12),
                          "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                        "Jun", "May","Apr","Mar","Feb","Jan"),3))
 # Assign the factors
 plotWeights$Month = factor(plotWeights$Month, levels = c("Jan","Feb","Mar","Apr",
                                                          "May","Jun","Jul","Aug",
                                                          "Sep","Oct","Nov","Dec"))
 # Create plot ready for grid.assign
 plot10 <- ggplot(plotWeights,aes(x=YearIntoPast,y=Weights,fill=Month)) +
   geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
   facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y") +
   ggtitle(paste0(" Ogle 3 year lag, 6 blocks OG Model (R2 = ",
                  signif(SAM_posterior_3_6_20200331_104434$R2,2), ", DIC = ",
                  signif(SAM_posterior_3_6_20200331_104434$DIC,4), ", MAE = ",
                  signif(SAM_posterior_3_6_20200331_104434$MAE,4),")"))
 grid.arrange(plot10, plot9)
 
 #****************************************************************************
 
 # This is for the HA1 dataset
 
 # timeblock = SpringBlock(3)
 # Nlag = timeblock$Nlag
 # block = timeblock$block
 # # # For the timeblock created, run the model
 # SAM(NEE[-(1:3),],Precip,Nlag,block,prior=FALSE)
 
 
 
 load("SAM_posterior_3_6_20200331_152929.Rdata")
 NPPobs = data.frame(Year=1:nrow(NEE)+1992,
                     NPP_obs = NEE[,2])
 NPPmod = data.frame(Year=1:length(SAM_posterior_3_6_20200331_152929$NPPmod$mean)+1995,
                     NPP_mod = SAM_posterior_3_6_20200331_152929$NPPmod$mean,
                     NPP_modmin = SAM_posterior_3_6_20200331_152929$NPPmod$min,
                     NPP_modmax = SAM_posterior_3_6_20200331_152929$NPPmod$max)
 plot11 <- ggplot(NPPobs) +
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
 plotWeights = data.frame("Weights"=as.vector(t(SAM_posterior_3_6_20200331_152929$monthlyWeights$mean)),
                          "YearIntoPast"=rep(0:2,each=12),
                          "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                        "Jun", "May","Apr","Mar","Feb","Jan"),3))
 # Assign the factors
 plotWeights$Month = factor(plotWeights$Month, levels = c("Jan","Feb","Mar","Apr",
                                                          "May","Jun","Jul","Aug",
                                                          "Sep","Oct","Nov","Dec"))
 # Create plot ready for grid.assign
 plot12 <- ggplot(plotWeights,aes(x=YearIntoPast,y=Weights,fill=Month)) +
   geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
   facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y") +
   ggtitle(paste0(" VAR 2 year lag, 8 blocks (R2 = ",
                  signif(SAM_posterior_3_6_20200331_152929$R2,2), ", DIC = ",
                  signif(SAM_posterior_3_6_20200331_152929$DIC,4), ", MAE = ",
                  signif(SAM_posterior_3_6_20200331_152929$MAE,4),")"))
 grid.arrange(plot12, plot11)
 
 
 
 
 Nlag = 3
 block = matrix(1:36,nrow=3,byrow=TRUE)
 # # # For the timeblock created, run the model
 # SAM(NEE[-(1:3),],Precip,Nlag,block,prior=FALSE)
 
 
 
 load("SAM_posterior_3_36_20200331_154314.Rdata")
 NPPobs = data.frame(Year=1:nrow(NEE)+1992,
                     NPP_obs = NEE[,2])
 NPPmod = data.frame(Year=1:length(SAM_posterior_3_36_20200331_154314$NPPmod$mean)+1995,
                     NPP_mod = SAM_posterior_3_36_20200331_154314$NPPmod$mean,
                     NPP_modmin = SAM_posterior_3_36_20200331_154314$NPPmod$min,
                     NPP_modmax = SAM_posterior_3_36_20200331_154314$NPPmod$max)
 plot13 <- ggplot(NPPobs) +
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
 plotWeights = data.frame("Weights"=as.vector(t(SAM_posterior_3_36_20200331_154314$monthlyWeights$mean)),
                          "YearIntoPast"=rep(0:2,each=12),
                          "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                        "Jun", "May","Apr","Mar","Feb","Jan"),3))
 # Assign the factors
 plotWeights$Month = factor(plotWeights$Month, levels = c("Jan","Feb","Mar","Apr",
                                                          "May","Jun","Jul","Aug",
                                                          "Sep","Oct","Nov","Dec"))
 # Create plot ready for grid.assign
 plot14 <- ggplot(plotWeights,aes(x=YearIntoPast,y=Weights,fill=Month)) +
   geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
   facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y") +
   ggtitle(paste0(" VAR 2 year lag, 8 blocks (R2 = ",
                  signif(SAM_posterior_3_36_20200331_154314$R2,2), ", DIC = ",
                  signif(SAM_posterior_3_36_20200331_154314$DIC,4), ", MAE = ",
                  signif(SAM_posterior_3_36_20200331_154314$MAE,4),")"))
 grid.arrange(plot14, plot13)
 