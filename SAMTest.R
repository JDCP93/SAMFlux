source("SAMFunction.R")

timeblock = SpringBlock(3)
Nlag = timeblock$Nlag
block = timeblock$block
Precip = Flux$Precip
NEE = Flux$NEE[-(1:Nlag),]
# For the timeblock created, run the model
SAM(Precip,NEE,Nlag,block,prior=FALSE)



load("SAM_posterior_3_6_20200325_165952.Rdata")
library(ggplot2)
library(gridExtra) # Source gridExtra for better plots
NPPobs = data.frame(Year=1:nrow(Flux$NEE)+2000,
                    NPP_obs = Flux$NEE[,2])
NPPmod = data.frame(Year=1:length(SAM_posterior_3_6_20200325_165952$NPPmod$mean)+2003,
                    NPP_mod = SAM_posterior_3_6_20200325_165952$NPPmod$mean,
                    NPP_modmin = SAM_posterior_3_6_20200325_165952$NPPmod$min,
                    NPP_modmax = SAM_posterior_3_6_20200325_165952$NPPmod$max)
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
plotWeights = data.frame("Weights"=as.vector(t(SAM_posterior_3_6_20200325_165952$monthlyWeights$mean)),
                         "YearIntoPast"=rep(0:2,each=12),
                         "Month"=rep(c("Dec","Nov","Oct","Sep", "Aug","Jul",
                                       "Jun", "May","Apr","Mar","Feb","Jan"),3))
# Assign the factors
plotWeights$Month = factor(plotWeights$Month, levels = c("Jan","Feb","Mar","Apr",
                                                         "May","Jun","Jul","Aug",
                                                         "Sep","Oct","Nov","Dec"))
# Create plot ready for grid.assign
plot1 <- ggplot(plotWeights,aes(x=YearIntoPast,y=Weights,fill=Month)) +
      geom_bar(stat="identity", position=position_dodge(), linetype = "solid",size=0.5,color="black") +
      facet_grid(.~YearIntoPast,scales = "free_x",switch = "x", space = "free_y") +
      ggtitle(paste0("3 year lag, spring/non-spring (R2 = ",
                     signif(SAM_posterior_3_6_20200325_165952$R2,2), ", DIC = ",
                     signif(SAM_posterior_3_6_20200325_165952$DIC,4), ", MAE = ",
                     signif(SAM_posterior_3_6_20200325_165952$MAE,4),")"))
grid.arrange(plot1, plot3)

