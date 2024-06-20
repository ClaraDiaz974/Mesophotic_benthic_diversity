

#Diaz et al. 2024 "Diverse and ecologically unique Mesophotic Coral Ecosystems in the central Indian Ocean"


######## Heatmap plots to show beta diversity for MCEs in Chagos #######

install.packages("Rtools")
install.packages("betapart")
library(betapart)
library(vegan)
library(ggplot2)
library(ggpubr)

#set working directory
setwd("C:/Users/cdiaz1/Desktop")

####### Beta diversity with depth for each study site #######

#Read in presence-absence data - morphospecies in columns, depths in rows - depth should be in order of increasing depth
#morphospecies catalogue is available on zenodo (Diaz et al.2023)
#IDR (Ile Des Rats), MA (Manta Alley) are the 2 study sites
data_MA <- read.csv("MA_Beta_P-A.csv", row.names = 1)
data_IDR <- read.csv("IDR_Beta_P-A.csv", row.names = 1)

#calculate beta diversity for P-A data
beta.pair(data_MA, index.family = "jaccard")
beta.pair(data_IDR, index.family = "jaccard")
#copy results from two rows above into excel and create new tables for heatmap - jtu = turnover matrix, jne = nestedness matrix, jac = dissimilarity matrix

#after making a file with all the beta results manually, here how to plot; based on Perez-Rosales et al., 2022 script on GitHub
betaresults <- read.csv(file = "C:/Users/cdiaz1/Desktop/Beta_results_IDR_MA.csv", header = T, row.names = 1)

betaresults$depth1 <- as.factor (betaresults$depth1)
betaresults$depth2 <- as.factor (betaresults$depth2)

betaresults$depth1 <- factor(betaresults$depth1, levels = c("15-20","30-40","60-70","80-90","110-120"))
betaresults$depth2 <- factor(betaresults$depth2, levels = c("30-40","60-70","80-90","110-120", "150-160"))

# IDR dissimilarity plot
IDRdiss <- ggplot(betaresults, aes(depth2, depth1, fill=IDR_dissimilarity))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.5, limit=c(0,1), breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), space="Lab", name="Jaccard's\ndissimilarity")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")+
  coord_fixed()+
  labs(x="",y="Depth (m)", title = "Total dissimilarity") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))

IDRdiss

# IDR turnover 
IDRturn <- ggplot(betaresults, aes(depth2, depth1, fill=IDR_turnover))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.5, limit=c(0,1),
                       breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),                     
                       space="Lab", name="Jaccard's\ndissimilarity")+ 
  # the\n is to have the title in two lines instead of one
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+ #to center the title, left by default
  theme(legend.position = "none")+  
  coord_fixed()+
  labs(x="",y="", title = "Turnover") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))
IDRturn

# IDR nestedness
IDRnest <- ggplot(betaresults, aes(depth2, depth1, fill=IDR_nestedness))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.5, limit=c(0,1), breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), space="Lab", name=" Jaccard's\ndissimilarity")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed()+
  labs(x="",y="", title = "Nestedness") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))

IDRnest

# MA dissimilarity
MAdiss <- ggplot(betaresults, aes(depth2, depth1, fill=MA_dissimilarity))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.5, limit=c(0,1), breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), space="Lab", name="Jaccard's\ndissimilarity")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  theme(legend.position = "none")+
  coord_fixed()+
  labs(x="",y="Depth (m)", title = "") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))

MAdiss

# MA turnover 
MAturn <- ggplot(betaresults, aes(depth2, depth1, fill=MA_turnover))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.5, limit=c(0,1),
                       breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),                     
                       space="Lab", name="Jaccard's\ndissimilarity")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  theme(legend.position = "none")+
  coord_fixed()+
  labs(x="",y="", title = "") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))
MAturn

# MA nestedness
MAnest <- ggplot(betaresults, aes(depth2, depth1, fill=MA_nestedness))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.5, limit=c(0,1), breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), space="Lab", name="Jaccard's\ndissimilarity")+ 
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="",y="", title = "") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
       strip.text = element_text(size = 10, colour="black"))

MAnest

#### Heat map figure
# Plot beta diversity graphs together: IDR and MA

figure <- ggarrange(IDRdiss ,  IDRturn , IDRnest , MAdiss , MAturn , MAnest,
                    nrow = 2, ncol= 3,
                    common.legend = TRUE, legend = "right")
figure #some aesthetic tweaks were made on illustrator or ppt afterwards. 

ggsave("C:/Users/cdiaz1/Desktop/Figure_heatmap_combined.jpg", figure,  width = 12, height = 8)


###### Beta diversity between study sites, for a same depth ######

setwd("C:/Users/cdiaz1/Desktop")

#Read in presence-absence data 
data_15m <- read.csv("BothSites_15m__Beta_P-A.csv", row.names = 1)
data_30m <- read.csv("BothSites_30m__Beta_P-A.csv", row.names = 1)
data_60m <- read.csv("BothSites_60m__Beta_P-A.csv", row.names = 1)
data_80m <- read.csv("BothSites_80m__Beta_P-A.csv", row.names = 1)
data_110m <- read.csv("BothSites_110m__Beta_P-A.csv", row.names = 1)
data_150m <- read.csv("BothSites_150m__Beta_P-A.csv", row.names = 1)

#calculate beta diversity for P-A data
beta.pair(data_15m, index.family = "jaccard")
beta.pair(data_30m, index.family = "jaccard")
beta.pair(data_60m, index.family = "jaccard")
beta.pair(data_80m, index.family = "jaccard")
beta.pair(data_110m, index.family = "jaccard")
beta.pair(data_150m, index.family = "jaccard")

#copy results from two rows above into excel and create new tables for heatmap - jtu = turnover matrix, jne = nestedness matrix, jac = dissimilarity matrix

Beta_results <-read.csv("Beta_div_Sites_forStackAreaPlot.csv", row.names = 1)
Beta_results

Figure <- ggplot(Beta_results, aes(x=factor(Depth), y=Diss_value)) +
  geom_bar(aes(fill=factor(Diss_type)), stat = "identity", width = 0.7) +
  scale_x_discrete(name ="Depth bands (m)", labels= c("15-20", "30-40", "60-70", "80-90", "110-120", "150-160")) +
  scale_y_continuous(name ="Jaccard's dissimilarity between study sites", limits=c(0,0.7), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)) +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "bottom", legend.text = element_text(size=10)) +
  guides(fill = guide_legend(ncol = 2, title = NULL)) +
  scale_fill_manual(values = c("#F05945","#5EAAA8"))
Figure

ggsave("D:/OneDrive/OneDrive - University of Plymouth/PhD Plymouth/Chapters/CH3 - Sp richness & Beta diversity/New catalogue/Beta diversity/Fig_BetaDiv_Sites_new_title.png", Figure, width = 5.3, height = 4.5)

##### end of the script #####

