

# Diaz et al. 2024 "Diverse and ecologically unique Mesophotic Coral Ecosystems in the central Indian Ocean"

# a iNext step by step script is also available in Daragmeh & El-Khaled 2021, "iNext4steps workflow for biodiversity assessment and comparison"


########## SPECIES RICHNESS CALCULATIONS EGMONT ATOLL ##############
install.packages("iNEXT")
install.packages("devtools")

library(devtools)
library(httr)
library(iNEXT)
library(ggplot2)
library(ggpubr)
library(grid)

########################################################################################

#set working directory

setwd("C:/Users/cdiaz1/Desktop")


#import data - species by site abundance matrix (species as rows, images as columns)
#IDR (Ile Des Rats), MA (Manta alley) are the two study sites
IDR_15_20m<-read.delim(file = "IDR-15-20.txt", header = TRUE, row.names = 1)
MA_15_20m<-read.delim(file = "MA-15-20.txt", header = TRUE, row.names = 1)
IDR_30_40m<-read.delim(file = "IDR-30-40.txt", header = TRUE, row.names = 1)
MA_30_40m<-read.delim(file = "MA-30-40.txt", header = TRUE, row.names = 1)
IDR_60_70m<-read.delim(file = "IDR-60-70.txt", header = TRUE, row.names = 1)
MA_60_70m<-read.delim(file = "MA-60-70.txt", header = TRUE, row.names = 1)
IDR_80_90m<-read.delim(file = "IDR-80-90.txt", header = TRUE, row.names = 1)
MA_80_90m<-read.delim(file = "MA-80-90.txt", header = TRUE, row.names = 1)
IDR_110_120m<-read.delim(file = "IDR-110-120.txt", header = TRUE, row.names = 1)
MA_110_120m<-read.delim(file = "MA-110-120.txt", header = TRUE, row.names = 1)
IDR_150_160m<-read.delim(file = "IDR-150-160.txt", header = TRUE, row.names = 1)
MA_150_160m<-read.delim(file = "MA-150-160.txt", header = TRUE, row.names = 1)

#transform data to iNext input format - sum of abundances for each OTU (Operational Taxonomic Unit)
IDR_15_20m_count<-(as.abucount(IDR_15_20m))
MA_15_20m_count<-(as.abucount(MA_15_20m))
IDR_30_40m_count<-(as.abucount(IDR_30_40m))
MA_30_40m_count<-(as.abucount(MA_30_40m))
IDR_60_70m_count<-(as.abucount(IDR_60_70m))
MA_60_70m_count<-(as.abucount(MA_60_70m))
IDR_80_90m_count<-(as.abucount(IDR_80_90m))
MA_80_90m_count<-(as.abucount(MA_80_90m))
IDR_110_120m_count<-(as.abucount(IDR_110_120m))
MA_110_120m_count<-(as.abucount(MA_110_120m))
IDR_150_160m_count<-(as.abucount(IDR_150_160m))
MA_150_160m_count<-(as.abucount(MA_150_160m))

# create a list of your matrices (renamed so the output looks nice)
IDR_MA_all =  list(IDR_15_20m = IDR_15_20m_count, MA_15_20m = MA_15_20m_count,
                   IDR_30_40m = IDR_30_40m_count, MA_30_40m = MA_30_40m_count,
                   IDR_60_70m = IDR_60_70m_count, MA_60_70m = MA_60_70m_count, 
                   IDR_80_90m = IDR_80_90m_count, MA_80_90m = MA_80_90m_count, 
                   IDR_110_120m = IDR_110_120m_count, MA_110_120m = MA_110_120m_count, 
                   IDR_150_160m = IDR_150_160m_count, MA_150_160m = MA_150_160m_count)

IDR_all=  list(IDR_15_20m = IDR_15_20m_count,
               IDR_30_40m = IDR_30_40m_count,
               IDR_60_70m = IDR_60_70m_count, 
               IDR_80_90m = IDR_80_90m_count, 
               IDR_110_120m = IDR_110_120m_count, 
               IDR_150_160m = IDR_150_160m_count)

MA_all=  list(MA_15_20m = MA_15_20m_count,
              MA_30_40m = MA_30_40m_count,
              MA_60_70m = MA_60_70m_count, 
              MA_80_90m = MA_80_90m_count, 
              MA_110_120m = MA_110_120m_count, 
              MA_150_160m = MA_150_160m_count)

# have a look at the raw data - Check q value (0 = richness, 1 = shannon diversity, 2 = simpson diversity)
out.IDR_15_20 <- iNEXT(IDR_15_20m_count, q = 0, endpoint = NULL, datatype="abundance")
out.MA_15_20 <- iNEXT(MA_15_20m_count, q = 0, endpoint = NULL, datatype="abundance")
out.IDR_30_40 <- iNEXT(IDR_30_40m_count, q = 0, endpoint = NULL, datatype="abundance")
out.MA_30_40 <- iNEXT(MA_30_40m_count, q = 0, endpoint = NULL, datatype="abundance")
out.IDR_60_70 <- iNEXT(IDR_60_70m_count, q = 0, endpoint = NULL, datatype="abundance")
out.MA_60_70 <- iNEXT(MA_60_70m_count, q = 0, endpoint = NULL, datatype="abundance")
out.IDR_80_90 <- iNEXT(IDR_80_90m_count, q = 0, endpoint = NULL, datatype="abundance")
out.MA_80_90 <- iNEXT(MA_80_90m_count, q = 0, endpoint = NULL, datatype="abundance")
out.IDR_110_120 <- iNEXT(IDR_110_120m_count, q = 0, endpoint = NULL, datatype="abundance")
out.MA_110_120 <- iNEXT(MA_110_120m_count, q = 0, endpoint = NULL, datatype="abundance")
out.IDR_150_160 <- iNEXT(IDR_150_160m_count, q = 0, endpoint = NULL, datatype="abundance")
out.MA_150_160 <- iNEXT(MA_150_160m_count, q = 0, endpoint = NULL, datatype="abundance")

out.IDR_MA_all_rich<- iNEXT(IDR_MA_all, q = 0, endpoint = NULL, datatype="abundance")
out.IDR_MA_all_shan<- iNEXT(IDR_MA_all, q = 1, endpoint = NULL, datatype="abundance")
out.IDR_MA_all_simp<- iNEXT(IDR_MA_all, q = 2, endpoint = NULL, datatype="abundance")

out.IDR_all_rich<- iNEXT(IDR_all, q = 0, endpoint = NULL, datatype="abundance")
out.IDR_all_shan<- iNEXT(IDR_all, q = 1, endpoint = NULL, datatype="abundance")
out.IDR_all_simp<- iNEXT(IDR_all, q = 2, endpoint = NULL, datatype="abundance")

out.MA_all_rich<- iNEXT(MA_all, q = 0, endpoint = NULL, datatype="abundance")
out.MA_all_shan<- iNEXT(MA_all, q = 1, endpoint = NULL, datatype="abundance")
out.MA_all_simp<- iNEXT(MA_all, q = 2, endpoint = NULL, datatype="abundance")

# to check if confidence intervals overlap
out.IDR_MA_15 <- iNEXT(IDR_MA_15, q = 0, endpoint = NULL, datatype="abundance")
out.IDR_MA_30 <- iNEXT(IDR_MA_30, q = 0, endpoint = NULL, datatype="abundance")
out.IDR_MA_80 <- iNEXT(IDR_MA_80, q = 0, endpoint = NULL, datatype="abundance")          
out.IDR_MA_110150<- iNEXT(IDR_MA_110150, q = 0, endpoint = NULL, datatype="abundance")
out.IDR_MA_1530 <- iNEXT(IDR_MA_1530, q = 0, endpoint = NULL, datatype="abundance")               
out.IDR_MA_153060 <- iNEXT(IDR_MA_153060, q = 0, endpoint = NULL, datatype="abundance") 
out.IDR_MA_6080 <- iNEXT(IDR_MA_6080, q = 0, endpoint = NULL, datatype="abundance") 

out.IDR_15_20
out.MA_15_20
out.IDR_30_40
out.MA_30_40
out.IDR_60_70
out.MA_60_70
out.IDR_80_90
out.MA_80_90
out.IDR_110_120
out.MA_110_120
out.IDR_150_160
out.MA_150_160

out.IDR_MA_all_rich
out.IDR_MA_all_shan
out.IDR_MA_all_simp

out.IDR_MA_15
out.IDR_MA_30
out.IDR_MA_80
out.IDR_MA_110150
out.IDR_MA_1530

###### make bar chart with the hill numbers ##### 
# species richness
q0 <- data.frame (site =rep(c("Ile des Rats", "Manta Alley"),each=6),
                  depth = rep(c("15-20m", "30-40m", "60-70m", "80-90m", "110-120m", "150-160m"),2),               
                  sprich =c(170,159,150,143,93,52,170,140,146,130,89,44))
q0

q0plot <- ggplot(data=q0, aes(x=depth, y=sprich, fill=site)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_discrete(limits=c("15-20m", "30-40m", "60-70m", "80-90m", "110-120m", "150-160m"))

q0plot

# shannon 
q1 <- data.frame (site =rep(c("Ile des Rats", "Manta Alley", "Sandes Seamount"),each=6),
                  depth = rep(c("15-20m", "30-40m", "60-70m", "80-90m", "110-120m", "150-160m"),3),               
                  shannon =c()
q1
                  
q1plot <- ggplot(data=q1, aes(x=depth, y=shannon, fill=site)) +
                    geom_bar(stat="identity", position=position_dodge()) +
                    scale_x_discrete(limits=c("15-20m", "30-40m", "60-70m", "80-90m", "110-120m", "150-160m"))
q1plot
                  
# simpson
q2 <- data.frame (site =rep(c("Ile des Rats", "Manta Alley", "Sandes Seamount"),each=6),
                  depth = rep(c("15-20m", "30-40m", "60-70m", "80-90m", "110-120m", "150-160m"),3),                                       simpson =c()))
q2

q2plot <- ggplot(data=q2, aes(x=depth, y=simpson, fill=site)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_discrete(limits=c("15-20m", "30-40m", "60-70m", "80-90m", "110-120m", "150-160m"))

q2plot


###### iNext plots #####
# type 1 = sample-size-based rarefaction/extrapolation curve; 
# type 2 = sample completeness curve 
# type 3 = coverage-based rarefaction/extrapolation curve.
# change the type number to make a different plot 


###########  per site (the two sites on a same plot is too much)  #############

## plot of IDR - sp richness - type 1
IDR_rich_type1 <- ggiNEXT(out.IDR_all_rich, type=1) + 
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom",
        legend.box = "vertical", 
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_colour_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B' # colours are not in the shading order because R takes it as alphabetical
  ),
  labels = c('IDR 110 - 120m','IDR 15 - 20m', 'IDR 150 - 160m', 'IDR 30 - 40m', 'IDR 60 - 70m', 'IDR 80 - 90m'
  )) +
  labs (x="Number of individuals", y= "Morphospecies diversity") +
  scale_fill_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B'
  )) +
  xlim(0, 10000) +
  ylim(0, 250)

IDR_rich_type1

## plot of MA - sp richness -type 1
MA_rich_type1 <- ggiNEXT(out.MA_all_rich, type=1) + 
  theme_bw(base_size = 18) + 
  theme(legend.position = "none", 
        legend.box = "vertical", 
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank())+ 
  scale_colour_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
 ),
  labels = c('MA 110 - 120m', 'MA 15 - 20m', 'MA 150 - 160m', 'MA 30 - 40m', 'MA 60 - 70m', 'MA 80 - 90m'
  )) +
  labs (x="Number of individuals", y= "Morphospecies diversity") +
  scale_fill_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
  )) +
  xlim(0, 10000) +
  ylim(0, 250)

MA_rich_type1

## plot of IDR - shannon - type 1
IDR_shan_type1 <- ggiNEXT(out.IDR_all_shan, type=1) + 
  theme_bw(base_size = 18) + 
  theme(legend.position = "none", 
        legend.box = "vertical",
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_colour_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B'
  ),
  labels = c('IDR 110 - 120m', 'IDR 15 - 20m', 'IDR 150 - 160m', 'IDR 30 - 40m', 'IDR 60 - 70m', 'IDR 80 - 90m'
  )) +
  labs (x="Number of individuals", y= "Morphospecies diversity") +
  scale_fill_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B'
  ))+
  xlim(0, 10000) +
  ylim(0, 50)

IDR_shan_type1

## plot of MA - shannon - type 1
MA_shan_type1 <- ggiNEXT(out.MA_all_shan, type=1) + 
  theme_bw(base_size = 18) + 
  theme(legend.position = "none", 
        legend.box = "vertical",
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank())+
  scale_colour_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
  ),
  labels = c('MA 110 - 120m', 'MA 15 - 20m', 'MA 150 - 160m', 'MA 30 - 40m', 'MA 60 - 70m', 'MA 80 - 90m'
  )) +
  labs (x="Number of individuals", y= "Morphospecies diversity") +
  scale_fill_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
  )) +
  xlim(0, 10000) +
  ylim(0, 50)

MA_shan_type1

## plot of IDR  - simpson - type 1
IDR_sim_type1 <- ggiNEXT(out.IDR_all_simp, type=1) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "none", 
        legend.box = "vertical", 
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_colour_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B'
  ),
  labels = c('IDR 110 - 120m', 'IDR 15 - 20m', 'IDR 150 - 160m', 'IDR 30 - 40m', 'IDR 60 - 70m', 'IDR 80 - 90m'
  )) +
  labs (x="Number of individuals", y= "Morphospecies diversity") +
  scale_fill_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B'                             
  )) +
  xlim(0, 10000) +
  ylim(0, 30)

IDR_sim_type1

## plot of MA - simpson - type 1
MA_sim_type1 <- ggiNEXT(out.MA_all_simp, type=1) + 
  theme_bw(base_size = 18) + 
  theme(legend.position = "none", 
        legend.box = "vertical", 
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank())+ 
  scale_colour_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
  ),
  labels = c('MA 110 - 120m', 'MA 15 - 20m', 'MA 150 - 160m', 'MA 30 - 40m', 'MA 60 - 70m', 'MA 80 - 90m'
  )) +
  labs (x="Number of individuals", y= "Morphospecies diversity") +
  scale_fill_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
  )) +
  xlim(0, 10000) +
  ylim(0, 30)

MA_sim_type1

### combine plots with ggarrange
plot_type_1 <- ggarrange(IDR_rich_type1, MA_rich_type1, IDR_shan_type1, MA_shan_type1, IDR_sim_type1, MA_sim_type1,ncol = 2, nrow =3, IDR_rich_type1 + rremove ("ylab") + rremove ("xlab"), MA_rich_type1 + rremove ("ylab") + rremove ("xlab"), IDR_shan_type1 + rremove ("ylab") + rremove ("xlab"), MA_shan_type1 + rremove ("ylab") + rremove ("xlab"), IDR_sim_type1 + rremove ("ylab") + rremove ("xlab"), MA_sim_type1 + rremove ("ylab") + rremove ("xlab"), align = "hv")

plot_type_1


#### plots for type 3: sample coverage and Cmax added, which is 96.4% in our case ####

##plot of IDR - sp richness - type 3
IDR_rich_type3 <- ggiNEXT(out.IDR_all_rich, type=3) + 
  theme_bw(base_size = 18) + 
  geom_vline(xintercept = 0.969, linetype="dashed",color = "black", size=.3)+
  theme(legend.position = "none", 
        legend.box = "vertical", 
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_colour_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B'
  ),
  labels = c('IDR 110 - 120m', 'IDR 15 - 20m', 'IDR 150 - 160m', 'IDR 30 - 40m', 'IDR 60 - 70m', 'IDR 80 - 90m'
  )) +
  labs (y="Morphospecies diversity", x= "Sample coverage") +
  scale_fill_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B'
  )) +
  ylim(0, 250)

IDR_rich_type3

## plot of MA - sp richness -type 3 
MA_rich_type3 <- ggiNEXT(out.MA_all_rich, type=3) + 
  theme_bw(base_size = 18) + 
  geom_vline(xintercept = 0.969, linetype="dashed",color = "black", size=.3)+
  theme(legend.position = "none", 
        legend.box = "vertical", 
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank())+ 
  scale_colour_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
  ),
  labels = c('MA 110 - 120m', 'MA 15 - 20m', 'MA 150 - 160m', 'MA 30 - 40m', 'MA 60 - 70m', 'MA 80 - 90m'
  )) +
  labs (y="Morphospecies diversity", x= "Sample coverage") +
  scale_fill_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
  )) +
  ylim(0, 250)

MA_rich_type3

## plot of IDR - shannon -type 3
IDR_shan_type3 <- ggiNEXT(out.IDR_all_shan, type=3) + 
  theme_bw(base_size = 18) + 
  geom_vline(xintercept = 0.969, linetype="dashed",color = "black", size=.3)+
  theme(legend.position = "none", 
        legend.box = "vertical",
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_colour_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B'
  ),
  labels = c('IDR 110 - 120m', 'IDR 15 - 20m', 'IDR 150 - 160m', 'IDR 30 - 40m', 'IDR 60 - 70m', 'IDR 80 - 90m'
  )) +
  labs (y="Morphospecies diversity", x= "Sample coverage") +
  scale_fill_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B'
  ))+
  ylim(0, 50)

IDR_shan_type3

## plot of MA - shannon - type 3
MA_shan_type3 <- ggiNEXT(out.MA_all_shan, type=3) + 
  theme_bw(base_size = 18) + 
  geom_vline(xintercept = 0.969, linetype="dashed",color = "black", size=.3)+
  theme(legend.position = "none", 
        legend.box = "vertical",
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank())+ 
  scale_colour_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
  ),
  labels = c('MA 110 - 120m', 'MA 15 - 20m', 'MA 150 - 160m', 'MA 30 - 40m', 'MA 60 - 70m', 'MA 80 - 90m'
  )) +
  labs (y="Morphospecies diversity", x= "Sample coverage") +
  scale_fill_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
  )) +
  ylim(0, 50)

MA_shan_type3

## plot of IDR - simpson - type 2
IDR_sim_type3 <- ggiNEXT(out.IDR_all_simp, type=3) + 
  theme_bw(base_size = 18) + 
  geom_vline(xintercept = 0.969, linetype="dashed",color = "black", size=.3)+
  theme(legend.position = "none", 
        legend.box = "vertical", 
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_colour_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B'
  ),
  labels = c('IDR 110 - 120m', 'IDR 15 - 20m', 'IDR 150 - 160m', 'IDR 30 - 40m', 'IDR 60 - 70m', 'IDR 80 - 90m'
  )) +
  labs (y="Morphospecies diversity", x= "Sample coverage") +
  scale_fill_manual(values = c('#006969','#A9DDDD', '#004343','#76BDBD','#429F9F','#008B8B'                            
  )) +
  ylim(0, 30)

IDR_sim_type3

## plot of MA - simpson - type 3
MA_sim_type3 <- ggiNEXT(out.MA_all_simp, type=3) +
  theme_bw(base_size = 18) +
  geom_vline(xintercept = 0.969, linetype="dashed",color = "black", size=.3)+
  theme(legend.position = "none", 
        legend.box = "vertical", 
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank())+ 
  scale_colour_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
  ),
  labels = c('MA 110 - 120m', 'MA 15 - 20m', 'MA 150 - 160m', 'MA 30 - 40m', 'MA 60 - 70m', 'MA 80 - 90m'
  )) +
  labs (y="Morphospecies diversity", x= "Sample coverage") +
  scale_fill_manual(values = c('#AF4F00','#FFCB9F','#703300','#FFAE6A','#FC9947','#E86900'
  )) +
  ylim(0, 30)

MA_sim_type3

### combine plots with ggarrange
plot_type_3<- ggarrange(IDR_rich_type3, MA_rich_type3, IDR_shan_type3, MA_shan_type3, IDR_sim_type3, MA_sim_type3,ncol = 2, nrow =3, IDR_rich_type3 + rremove ("ylab") + rremove ("xlab"), MA_rich_type3 + rremove ("ylab") + rremove ("xlab"), IDR_shan_type3 + rremove ("ylab") + rremove ("xlab"), MA_shan_type3 + rremove ("ylab") + rremove ("xlab"), IDR_sim_type3 + rremove ("ylab") + rremove ("xlab"), MA_sim_type3 + rremove ("ylab") + rremove ("xlab"), align = "hv")

plot_type_3



########## Estimate D #########
# for coverage base, m= nb of individuals needed; order 0 (sp rich), 1 (shannon), 2(simpson); sc = level; qD = result of sp rich, shannon or simpson; LCL= lower confidence level; UCL = upper confidence level.

estimateD(IDR_MA_15, "abundance", base = "coverage", level = 0.99, conf = 0.95)
estimateD(IDR_MA_30, "abundance", base = "coverage", level = 0.99, conf = 0.95)
estimateD(IDR_MA_6080, "abundance", base = "coverage", level = 0.99, conf = 0.95)
estimateD(IDR_MA_110150, "abundance", base = "coverage", level = 0.99, conf = 0.95)

##### transform the tables into incidence frequency to know the sample completeness #####
# based on the script from Daragmeh and El-Khaled 2021

IDR_15_20m_incfreq<-read.delim(file = "IDR-15-20-freq.txt", header = TRUE, row.names = 1)
IDR_30_40m_incfreq<-read.delim(file = "IDR-30-40-freq.txt", header = TRUE, row.names = 1)

IDR_trial_list =  list(IDR_15_20m_incfreq ,
                       IDR_30_40m_incfreq)
names(IDR_trial_list)<-c("IDR15","IDR30")

IDR_MA_all_list =  list(IDR_15_20m , MA_15_20m,
                        IDR_30_40m , MA_30_40m ,
                        IDR_60_70m , MA_60_70m , 
                        IDR_80_90m , MA_80_90m , 
                        IDR_110_120m , MA_110_120m , 
                        IDR_150_160m , MA_150_160m )


# Transform the list objects into an incidence frequency format as required for iNEXT's online tools
# This is done using iNEXT's as.incfreq function

IDR_15_20m_incfreq<-read.delim(file = "IDR-15-20-freq.txt", header = TRUE, row.names = 1)
IDR_30_40m_incfreq<-read.delim(file = "IDR-30-40-freq.txt", header = TRUE, row.names = 1)
IDR_60_70m_incfreq<-read.delim(file = "IDR-60-70-freq.txt", header = TRUE, row.names = 1)
IDR_80_90m_incfreq<-read.delim(file = "IDR-80-90-freq.txt", header = TRUE, row.names = 1)
IDR_110_120m_incfreq<-read.delim(file = "IDR-110-120-freq.txt", header = TRUE, row.names = 1)
IDR_150_160m_incfreq<-read.delim(file = "IDR-150-160-freq.txt", header = TRUE, row.names = 1)
MA_15_20m_incfreq<-read.delim(file = "MA-15-20-freq.txt", header = TRUE, row.names = 1)
MA_30_40m_incfreq<-read.delim(file = "MA-30-40-freq.txt", header = TRUE, row.names = 1)
MA_60_70m_incfreq<-read.delim(file = "MA-60-70-freq.txt", header = TRUE, row.names = 1)
MA_80_90m_incfreq<-read.delim(file = "MA-80-90-freq.txt", header = TRUE, row.names = 1)
MA_110_120m_incfreq<-read.delim(file = "MA-110-120-freq.txt", header = TRUE, row.names = 1)
MA_150_160m_incfreq<-read.delim(file = "MA-150-160-freq.txt", header = TRUE, row.names = 1)

IDR_MA_all_list =  list(IDR_15_20m_incfreq , MA_15_20m_incfreq,
                        IDR_30_40m_incfreq , MA_30_40m_incfreq ,
                        IDR_60_70m_incfreq , MA_60_70m_incfreq , 
                        IDR_80_90m_incfreq , MA_80_90m_incfreq , 
                        IDR_110_120m_incfreq , MA_110_120m_incfreq , 
                        IDR_150_160m_incfreq , MA_150_160m_incfreq )
names(IDR_MA_all_list)<-c("IDR15","MA15","IDR30","MA30","IDR60","MA60","IDR80","MA80","IDR110","MA110","IDR150","MA150")

IDR_MA_all_incfreq <- lapply(IDR_MA_all_list, as.incfreq)

# Save this object as .rds file for further iNEXT analysis
saveRDS(IDR_MA_all_incfreq,"IDR_MA_all_incfreq.rds")

# Write a tab-delimited txt-file which  will be formated outside of RStudio using excel
write.table(IDR_MA_all_incfreq,"IDR_MA_all_inc.txt",sep="\t",col.names = NA)


### Sample-sized based rarefaction and extrapolation curves ###
# Loading the .rds-file saved  
IDR_MA_all_incfreq<-readRDS("IDR_MA_all_incfreq.rds")

# Compute the curves with the iNEXT function and 500 bootstraps
IDR_MA_all_inext<-iNEXT(IDR_MA_all_incfreq,datatype = "incidence_freq",q=c(0,1,2),nboot=500)

# Write the results of sample-size based rarefaction and extrapolation to a file
write.table(IDR_MA_all_inext$iNextEst,"IDR_MA_all_sample_size_data.txt",sep="\t",col.names = NA)

