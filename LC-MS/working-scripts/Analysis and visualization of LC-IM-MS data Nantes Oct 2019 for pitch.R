### Analysis and visualization of LC-IM-MS data, Nantes Oct 2019

# general information about the data --------------------------------------

# Intensities manually extracted from MassLynx (see power point presentation
# with results), binned by m/z, retention time and drift time
# peak at m/z 851.27, DT 6.5 and RT 25.4-25.9 has weird shape -> RT bin 
# enlarged on 03. May 2021. In addition, DT estimate changed for no enzyme, 
# 30 min: at RT 15.5 DT of 6.7 changed to 6.5, at RT 23.6 DT of 6.5 changed
# to 6.3

# packages ----------------------------------------------------------------

require(ggpubr)
require(dplyr)
require(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(tidyverse)
library(wesanderson)
library(reshape2)


fonts <- list(
  sans = "Helvetica",
  mono = "Consolas",
  `Times New Roman` = "DejaVu Serif"
)
font_family<-"sans"
font_size<-8
font_size_titles<-10

'%!in%' <- function(x,y)!('%in%'(x,y))

# import data -------------------------------------------------------------

all<-read.table("C:/Users/admin/ownCloud/Deep-water kelp/Nantes2019/Results/LC-IM-MS_mod03May2021.txt",h=T,dec=",")


lam_stds<-all[all$Treatment=="Standard",]
lam_stds$RTandDT<-paste("RT=",lam_stds$RT_min," min DT=",lam_stds$DT_ms," ms",sep="")


wo_stds<-all[all$Treatment!="Standard",]
wo_stds$RTandDT<-paste("RT=",wo_stds$RT_min," min DT=",wo_stds$DT_ms," ms",sep="")



#filtering based on abundances in time series
wo_stds_filtered<-wo_stds%>%
  filter(!RTandDT%in%c("RT=12.6 min DT=6 ms", "RT=24.2 min DT=6.6 ms", "RT=18.7 min DT=6.8 ms"))

wo_stds_filtered<-wo_stds_filtered[order(wo_stds_filtered$Intensity_IM,decreasing=T),]
wo_stds_filtered<-wo_stds_filtered[order(wo_stds_filtered$DP),]


palette_9extractions<-c("#74A089","#FDDDA0","#C93312","#0B775E", "#FAD510" ,"#F2300F", "#899DA4", "#FAEFD1", "#DC863B")

Fig_samples_filtered<-ggplot()+
  geom_point(data=wo_stds_filtered,aes(x=Position,y=Intensity_IM,fill=RTandDT),pch=21, size=3)+
  scale_y_log10()+
  facet_wrap(facets=wo_stds_filtered$DP,nrow=1)+
  scale_fill_manual(values=c(palette_9extractions[1:5],palette_9extractions[1:7],palette_9extractions[1:7]))+
  guides(ncol=1)+
  theme(axis.text.x= element_text(colour="#000000", size = font_size, family = font_family),
        axis.text.y = element_text(colour ="#000000", size = font_size, family = font_family),
        axis.title.y = element_text(colour="#000000", size = font_size, family = font_family),
        axis.title.x = element_text(colour="#000000", size = font_size, family = font_family),
        legend.title= element_text(colour="#000000",face="bold", size = font_size, family = font_family),
        legend.text = element_text(colour="#000000", size = font_size, family = font_family),
        #        legend.position="none",
        legend.box.background = element_rect(colour="black"),
        panel.background = element_rect(colour="#000000", fill=NA),
        plot.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour="#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())


# check standards ---------------------------------------------------------

stds<-all[which(all$Treatment %in% c("Standard","none")),]
stds$RTandDT<-paste("RT=",stds$RT_min," min DT=",stds$DT_ms," ms",sep="")

ggplot(stds[stds$DP==5,])+
  geom_bar(aes(x=RTandDT, y=Intensity_MZ, fill=Treatment),stat="identity",position=position_dodge2())

overSN<-stds[stds$Intensity_IM>=7000,]
ggplot(overSN[overSN$DP==3,])+
  geom_bar(aes(x=RTandDT, y=Intensity_MZ, fill=Treatment),stat="identity",position=position_dodge2())



# compare std vs sample features ------------------------------------------

comparison<-all
comparison$RTandDT<-paste("RT=",comparison$RT_min," min DT=",comparison$DT_ms," ms",sep="")

comparison_end<-comparison[which(comparison)]
ggplot(comparison[comparison$DP==3,])+
  geom_bar(aes(x=RTandDT, y=Intensity_MZ, fill=Treatment),stat="identity",position=position_dodge2())




# margi's approach --------------------------------------------------------


# highest features in digested samples
high_digests<-wo_stds %>% 
  group_by(Incubation_time_min, Treatment, DP) %>%
  top_n(1, Intensity_MZ)

# highest features in standards
lam_stds %>% 
  group_by(DP) %>% 
  top_n(1, Intensity_MZ)


# the feature RT=14.2 min DT=4.1 ms m/z=527.16 represents the canonical product (DP3) [M+Na]
# the feature RT=17.8 min DT=5.3 ms m/z=689.22 represents the canonical product (DP4) [M+Na]
# the feature RT=25.7 min DT=6.5 ms m/z=851.28 represents the canonical product (DP5) [M+Na]

# remove non-increasing features ------------------------------------------

wo_stds<-wo_stds[wo_stds$RTandDT!="RT=12.6 min DT=6 ms",]

wo_stds<-wo_stds[wo_stds$RTandDT!="RT=15.5 min DT=6.5 ms",]

# check features statistically --------------------------------------------


wo_stds$Position<-NA
wo_stds$Position[wo_stds$Treatment=="none"]<-NA
wo_stds$Position[wo_stds$Treatment=="GH16" & wo_stds$Incubation_time_min==0]<-1
wo_stds$Position[wo_stds$Treatment=="GH16" & wo_stds$Incubation_time_min==10]<-10
wo_stds$Position[wo_stds$Treatment=="GH16" & wo_stds$Incubation_time_min==20]<-20
wo_stds$Position[wo_stds$Treatment=="GH16" & wo_stds$Incubation_time_min==30]<-30



# Linear model canonical product ------------------------------------------

canonical_prods<-wo_stds[which(wo_stds$RTandDT %in% c("RT=14.2 min DT=4.1 ms", 
                                                      "RT=17.8 min DT=5.3 ms",
                                                      "RT=25.7 min DT=6.5 ms")),]

lm_canonical_DP3<-lm(log(Intensity_MZ)~Position, canonical_prods[canonical_prods$DP=="3",])
lm_canonical_DP4<-lm(log(Intensity_MZ)~Position, canonical_prods[canonical_prods$DP=="4",])


# natural log transform (for first order products)


wo_stds$DPlabel<-paste("DP", wo_stds$DP, sep="")
wo_stds$DPlabel<-factor(wo_stds$DPlabel, levels = c("DP5","DP4","DP3"))

wo_stds$RTandDT<-factor(wo_stds$RTandDT, levels = c())

png(filename="C:/Users/admin/ownCloud/Deep-water kelp/Presentations/Thesis committee meeting 11May2021/Fig3a_sideprods_06May2021.png",
    height = 12, width=18.9, units = "cm", res = 200)
ggplot()+
  geom_point(data=wo_stds,aes(x=Position,y=log(Intensity_MZ),fill=RTandDT,), pch=21,size=3, alpha=0.5)+
  scale_fill_discrete(name="Features", guide=guide_legend(reverse=T))+
#  scale_y_continuous(trans = scales::log_trans())+

  stat_cor(data=wo_stds[which(wo_stds$RTandDT %in% c("RT=14.2 min DT=4.1 ms", 
                                                     "RT=17.8 min DT=5.3 ms",
                                                     "RT=25.7 min DT=6.5 ms")),],
           mapping=aes(x=Position,y=log(Intensity_MZ), label = paste(..rr.label..)),
           r.accuracy = 0.01, 
             label.y = 18,
             label.x = 0,
             size = 3.5,
             family = font_family) + 
    stat_cor(data=wo_stds[which(wo_stds$RTandDT %in% c("RT=14.2 min DT=4.1 ms", 
                                                       "RT=17.8 min DT=5.3 ms",
                                                       "RT=25.7 min DT=6.5 ms")),],
             mapping=aes(x=Position,y=log(Intensity_MZ), label = paste(..p.label..)),
             p.accuracy = 0.001,
             label.y = 17,
             label.x = 0,
             size = 3.5,
             family = font_family) + #add p value
    geom_smooth(data=wo_stds[which(wo_stds$RTandDT %in% c("RT=14.2 min DT=4.1 ms", 
                                                          "RT=17.8 min DT=5.3 ms",
                                                          "RT=25.7 min DT=6.5 ms")),],
                mapping=aes(x=Position,y=log(Intensity_MZ)),
                method = "lm",
                se = FALSE,
                colour = "#616468",
                alpha = 0.7,
                size = 1) +
  xlab("Incubation time (min)")+
  ylab("ln(Intensity)")+
  facet_grid(cols=vars(DPlabel))+
  guides(ncol=1)+
  theme(axis.text.x= element_text(colour="#000000", size = font_size, family = font_family),
        axis.text.y = element_text(colour ="#000000", size = font_size, family = font_family),
        axis.title.y = element_text(colour="#000000", size = font_size_titles, family = font_family),
        axis.title.x = element_text(colour="#000000", size = font_size_titles, family = font_family),
        legend.title= element_text(colour="#000000", size = font_size_titles, family = font_family,
                                   hjust = 0.5),
        legend.text = element_text(colour="#000000", size = font_size, family = font_family),
        legend.key = element_blank(),
        legend.box.background = element_blank(),
        panel.background = element_rect(colour="#000000", fill=NA),
        plot.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour="#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(colour="#000000", size=font_size, family = font_family))
dev.off()


#### figure 3 again, wider but shorter -----

png(filename="C:/Users/admin/ownCloud/Deep-water kelp/Presentations/Thesis committee meeting 11May2021/Fig3a_sideprods_06May2021_shorter.png",
    height = 8, width=18.9, units = "cm", res = 200)
ggplot()+
  geom_point(data=wo_stds,aes(x=Position,y=log(Intensity_MZ),fill=RTandDT,), pch=21,size=2.5, alpha=0.5)+
  scale_fill_discrete(name="Features", guide=guide_legend(reverse=T))+
  scale_y_continuous(limits = c(7,17))+
  
  stat_cor(data=wo_stds[which(wo_stds$RTandDT %in% c("RT=14.2 min DT=4.1 ms", 
                                                     "RT=17.8 min DT=5.3 ms",
                                                     "RT=25.7 min DT=6.5 ms")),],
           mapping=aes(x=Position,y=log(Intensity_MZ), label = paste(..rr.label..)),
           r.accuracy = 0.01, 
           label.y = 16.7,
           label.x = 11,
           size = 3,
           family = font_family) + # r squared
  stat_cor(data=wo_stds[which(wo_stds$RTandDT %in% c("RT=14.2 min DT=4.1 ms", 
                                                     "RT=17.8 min DT=5.3 ms",
                                                     "RT=25.7 min DT=6.5 ms")),],
           mapping=aes(x=Position,y=log(Intensity_MZ), label = paste(..p.label..)),
           p.accuracy = 0.001,
           label.y = 16,
           label.x = 11,
           size = 3,
           family = font_family) + #add p value
  geom_smooth(data=wo_stds[which(wo_stds$RTandDT %in% c("RT=14.2 min DT=4.1 ms", 
                                                        "RT=17.8 min DT=5.3 ms",
                                                        "RT=25.7 min DT=6.5 ms")),],
              mapping=aes(x=Position,y=log(Intensity_MZ)),
              method = "lm",
              se = FALSE,
              colour = "#616468",
              alpha = 0.7,
              size = 1) +
  xlab("Incubation time (min)")+
  ylab("ln(Intensity)")+
  facet_grid(cols=vars(DPlabel))+
  theme(axis.text.x= element_text(colour="#000000", size = font_size, family = font_family),
        axis.text.y = element_text(colour ="#000000", size = font_size, family = font_family),
        axis.title.y = element_text(colour="#000000", size = font_size_titles, family = font_family),
        axis.title.x = element_text(colour="#000000", size = font_size_titles, family = font_family),
        legend.title= element_text(colour="#000000", size = font_size_titles, family = font_family,
                                   hjust = 0.5),
        legend.text = element_text(colour="#000000", size = font_size, family = font_family),
        legend.key = element_blank(),
        legend.box.background = element_blank(),
        panel.background = element_rect(colour="#000000", fill=NA),
        plot.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour="#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(colour="#000000", size=font_size_titles, family = font_family))+
  guides(fill=guide_legend(ncol=2))
dev.off()

#double transform natural log (for second order products)

Fig_samples_points_second<-ggplot()+
  geom_point(data=wo_stds,aes(x=Position,y=log(log(Intensity_MZ)),fill=RTandDT), pch=21,size=3, alpha=0.5)+
  #  scale_y_continuous(trans = scales::log_trans())+
  
  stat_cor(data=wo_stds[which(wo_stds$RTandDT %!in% c("RT=14.2 min DT=4.1 ms", 
                                                     "RT=17.8 min DT=5.3 ms",
                                                     "RT=25.7 min DT=6.5 ms")),],
           mapping=aes(x=Position,y=log(log(Intensity_MZ)), label = paste(..rr.label..), group=RTandDT),
           r.accuracy = 0.01, 
           label.y = c(3, 2.95, 2.9, 2.85, 3, 2.95, 2.9, 2.85, 2.8, 3, 2.95, 2.9, 2.85, 2.8, 2.75),
           label.x = 0,
           size = 3.5,
           position = position_dodge2(),
           family = "Avenir") + 
  stat_cor(data=wo_stds[which(wo_stds$RTandDT %!in% c("RT=14.2 min DT=4.1 ms", 
                                                     "RT=17.8 min DT=5.3 ms",
                                                     "RT=25.7 min DT=6.5 ms")),],
           mapping=aes(x=Position,y=log(log(Intensity_MZ)), label = paste(..p.label..), group=RTandDT),
           p.accuracy = 0.001,
           label.y = c(3, 2.95, 2.9, 2.85, 3, 2.95, 2.9, 2.85, 2.8, 3, 2.95, 2.9, 2.85, 2.8, 2.75),
           label.x = 20,
           size = 3.5,
           family = "Avenir") + #add p value
  geom_smooth(data=wo_stds[which(wo_stds$RTandDT %!in% c("RT=14.2 min DT=4.1 ms", 
                                                        "RT=17.8 min DT=5.3 ms",
                                                        "RT=25.7 min DT=6.5 ms")),],
              mapping=aes(x=Position,y=log(log(Intensity_MZ)), group=RTandDT),
              method = "lm",
              se = FALSE,
              colour = "#616468",
              alpha = 0.7,
              size = 1) +
  facet_grid(cols=vars(DP))+
  guides(ncol=1)+
  theme_light()

  
  geom_point(data=wo_stds,aes(x=Position,y=log(log(Intensity_MZ)),fill=RTandDT,),pch=21,size=3, alpha=0.5)+
  scale_y_continuous(trans = scales::log_trans())+
  #  scale_x_continuous(trans = "sqrt")+
  facet_wrap(facets=wo_stds$DP,nrow=1)+
  guides(ncol=1)+
  theme_light()


Fig_samples_lines<-ggplot()+
  geom_line(data=wo_stds,aes(x=Position,y=Intensity_MZ,color=RTandDT),alpha=0.8,size=1.2)+
  scale_y_continuous(trans = scales::log_trans())+
  #  scale_x_continuous(trans = "sqrt")+
  facet_wrap(facets=wo_stds$DP,nrow=1)+
  guides(ncol=1)



# assessment multivariate -------------------------------------------------

wo_stds_cast<-dcast(wo_stds, Position ~ RTandDT, value.var = "Intensity_MZ")

# MANOVA test
res.man <- manova(cbind(wo_stds_cast[,4], wo_stds_cast[,3]) ~ wo_stds_cast[,1])
summary(res.man)

# I actually think this problem is not multivaryate but many individual test are
# necessary

lm_1<-lm(cbind(sqrt(wo_stds_cast[,5]),sqrt(wo_stds_cast[,6])) ~sqrt(wo_stds_cast[,1]))
summary(lm_1)
plot(lm_1)

# Trying with polynominal regression
#better: find proper function that describes enzyme reactions, fit nls (non linear least squares)

beta_1_3_DP3<-wo_stds[wo_stds$RT_min==14.2 & wo_stds$Position>=0,]

beta_1_3_DP3$Position2 = beta_1_3_DP3$Position^2
beta_1_3_DP3$Position3 = beta_1_3_DP3$Position^3
poly_regressor = lm(formula = Intensity_IM ~ Position+Position2,data = beta_1_3_DP3)

summary(poly_regressor)

X_grid = seq(min(beta_1_3_DP3$Position), max(beta_1_3_DP3$Position), 0.1)
ggplot() +
  geom_point(aes(x = beta_1_3_DP3$Position, y = beta_1_3_DP3$Intensity_IM),colour = 'black') +
  geom_line(aes(x = X_grid, y = predict(poly_regressor, newdata = data.frame(Position = X_grid, Position2=X_grid^2))),colour = 'red')+
  ggtitle('Polynomial Regression')

beta_1_3_DP4<-wo_stds[wo_stds$RT_min==17.8 & wo_stds$Position>=0,]

beta_1_3_DP4$Position2 = beta_1_3_DP4$Position^2
beta_1_3_DP4$Position3 = beta_1_3_DP4$Position^3
poly_regressor = lm(formula = Intensity_IM ~ Position+Position2,data = beta_1_3_DP4)

summary(poly_regressor)

X_grid = seq(min(beta_1_3_DP4$Position), max(beta_1_3_DP4$Position), 0.1)
ggplot() +
  geom_point(aes(x = beta_1_3_DP4$Position, y = beta_1_3_DP4$Intensity_IM),colour = 'black') +
  geom_line(aes(x = X_grid, y = predict(poly_regressor, newdata = data.frame(Position = X_grid, Position2=X_grid^2))),colour = 'red')+
  ggtitle('Polynomial Regression')


beta_1_3_DP5<-wo_stds[wo_stds$RT_min==25.9 & wo_stds$Position>=0,]

beta_1_3_DP5$Position2 = beta_1_3_DP5$Position^2
beta_1_3_DP5$Position3 = beta_1_3_DP5$Position^3
poly_regressor = lm(formula = Intensity_IM ~ Position+Position3,data = beta_1_3_DP5)

summary(poly_regressor)

X_grid = seq(min(beta_1_3_DP5$Position), max(beta_1_3_DP5$Position), 0.1)
ggplot() +
  geom_point(aes(x = beta_1_3_DP5$Position, y = beta_1_3_DP5$Intensity_IM),colour = 'black') +
  geom_line(aes(x = X_grid, y = predict(poly_regressor, newdata = data.frame(Position = X_grid, Position2=X_grid^2))),colour = 'red')+
  ggtitle('Polynomial Regression')


### FIGURE FOR POMPU meeting in Berlin November 2019

wo_stds_filtered<-wo_stds_filtered[order(wo_stds_filtered$DT_ms),]

Fig_samples_filtered_DP3<-ggplot()+
  geom_line(data=wo_stds_filtered[wo_stds_filtered$DP==3,],aes(x=Position,y=Intensity_IM,color=RTandDT),alpha=0.8, size=1.5)+
  scale_y_log10(name="Intensity (A.U.)",limits=c(300,5000000))+
  scale_color_manual(values=wes_palette("Darjeeling1", n=5) )+
  scale_x_continuous(name="Incubation time",breaks=c(-10,0,10,20,30),labels=c("not\ndigested","0 min","10 min","20 min", "30 min"))+
  guides(ncol=4)+
  theme(axis.text.x= element_text(colour="#000000", size = font_size, family = font_family),
        axis.text.y = element_text(colour ="#000000", size = font_size, family = font_family),
        axis.title.y = element_text(colour="#000000", size = font_size, family = font_family),
        axis.title.x = element_text(colour="#000000", size = font_size, family = font_family),
        legend.title= element_text(colour="#000000",face="bold", size = font_size, family = font_family),
        legend.text = element_text(colour="#000000", size = font_size, family = font_family),
        legend.position="none",
        legend.box.background = element_rect(colour="black"),
        panel.background = element_rect(colour="#000000", fill=NA),
        plot.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour="#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  ggtitle(label="Trisaccharides")

Fig_samples_filtered_DP4<-ggplot()+
  geom_line(data=wo_stds_filtered[wo_stds_filtered$DP==4,],aes(x=Position,y=Intensity_IM,color=RTandDT),alpha=0.8, size=1.5)+
  scale_y_log10(name="Intensity (A.U.)",limits=c(300,5000000))+
  scale_color_manual(values=c(wes_palette("Darjeeling1", n=1), wes_palette("Darjeeling2", n=5)))+
  guides(ncol=4)+
  scale_x_continuous(name="Incubation time",breaks=c(-10,0,10,20,30),labels=c("not\ndigested","0 min","10 min","20 min", "30 min"))+
  theme(axis.text.x= element_text(colour="#000000", size = font_size, family = font_family),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour="#000000", size = font_size, family = font_family),
        legend.title= element_text(colour="#000000",face="bold", size = font_size, family = font_family),
        legend.text = element_text(colour="#000000", size = font_size, family = font_family),
        legend.position="none",
        legend.box.background = element_rect(colour="black"),
        panel.background = element_rect(colour="#000000", fill=NA),
        plot.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour="#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  ggtitle(label="Tetrasaccharides")

Fig_samples_filtered_DP5<-ggplot()+
  geom_point(data=wo_stds_filtered[wo_stds_filtered$DP==5,],aes(x=Position,y=Intensity_IM,fill=RTandDT),pch=21, size=3)+
  scale_y_log10(name="Intensity (A.U.)",limits=c(300,5000000))+
  scale_fill_manual(values=c("#74A089","#FDDDA0","#0B775E", "#FAEFD1","#F2300F", "#899DA4",  "#FAD510","#C93312"))+
  scale_x_continuous(name="Incubation time",breaks=c(-10,0,10,20,30),labels=c("not\ndigested","0 min","10 min","20 min", "30 min"))+
  theme(axis.text.x= element_text(colour="#000000", size = font_size, family = font_family),
        axis.text.y = element_text(colour ="#000000", size = font_size, family = font_family),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour="#000000", size = font_size, family = font_family),
        legend.title= element_text(colour="#000000",face="bold", size = font_size, family = font_family),
        legend.text = element_text(colour="#000000", size = font_size, family = font_family),
        legend.position="none",
        legend.box.background = element_rect(colour="black"),
        panel.background = element_rect(colour="#000000", fill=NA),
        plot.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour="#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  ggtitle(label = "Pentasaccharides")



DP7_diversification<-plot_grid(Fig_samples_filtered_DP3,Fig_samples_filtered_DP4, align="h",nrow=1,rel_widths = c(1.2,1))
save_plot(filename = "C:/Users/admin/ownCloud/Deep-water kelp/Presentations/SAB2021_pitch/Fig_DP7_diversification.tiff",DP7_diversification,base_height = 5,base_width = 12)          


# all lengths in one

# order features by mean int
wo_stds_filtered$RTandDT<-fct_reorder(.f=wo_stds_filtered$RTandDT,.x=wo_stds_filtered$Intensity_IM, .fun = max)

Fig_samples_filtered_allDPs<-ggplot()+
  geom_line(data=wo_stds_filtered,aes(x=Position,y=Intensity_IM,color=RTandDT),alpha=0.6, size=1)+
  scale_y_log10(name="Intensity (A.U.)",limits=c(300,5000000))+
  scale_color_manual(values=c(rep("white", 16), rep("black",3)))+
  scale_x_continuous(name="Incubation time",breaks=c(-10,0,10,20,30),labels=c("not\ndigested","0 min","10 min","20 min", "30 min"))+
  guides(ncol=4)+
  theme(axis.text.x= element_text(colour="#000000", size = font_size, family = font_family),
        axis.text.y = element_text(colour ="#000000", size = font_size, family = font_family),
        axis.title.y = element_text(colour="#000000", size = font_size, family = font_family),
        axis.title.x = element_text(colour="#000000", size = font_size, family = font_family),
        legend.title= element_text(colour="#000000",face="bold", size = font_size, family = font_family),
        legend.text = element_text(colour="#000000", size = font_size, family = font_family),
        legend.position="none",
        legend.box.background = element_rect(colour="black"),
        panel.background = element_rect(colour="#000000", fill=NA),
        plot.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour="#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())


### Figure without log 10

Fig_samples_filtered_DP3_nolog10<-ggplot()+
  geom_point(data=wo_stds_filtered[wo_stds_filtered$DP==3,],aes(x=Position,y=Intensity_IM,fill=RTandDT),pch=21, size=3)+
  scale_y_continuous(name="Intensity (A.U.)",limits=c(300,5000000))+
  scale_fill_manual(values=c("#FDDDA0", "#F2300F","#C93312","#FAD510","#0B775E" ))+
  scale_x_continuous(name="Incubation time",breaks=c(-10,0,10,20,30),labels=c("not\ndigested","0 min","10 min","20 min", "30 min"))+
  guides(ncol=4)+
  theme(axis.text.x= element_text(colour="#000000", size = font_size, family = font_family),
        axis.text.y = element_text(colour ="#000000", size = font_size, family = font_family),
        axis.title.y = element_text(colour="#000000", size = font_size, family = font_family),
        axis.title.x = element_text(colour="#000000", size = font_size, family = font_family),
        legend.title= element_text(colour="#000000",face="bold", size = font_size, family = font_family),
        legend.text = element_text(colour="#000000", size = font_size, family = font_family),
        legend.position="none",
        legend.box.background = element_rect(colour="black"),
        panel.background = element_rect(colour="#000000", fill=NA),
        plot.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour="#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  ggtitle(label="Trisaccharides")

Fig_samples_filtered_DP4_nolog10<-ggplot()+
  geom_point(data=wo_stds_filtered[wo_stds_filtered$DP==4,],aes(x=Position,y=Intensity_IM,fill=RTandDT),pch=21, size=3)+
  scale_y_continuous(name="Intensity (A.U.)",limits=c(300,5000000))+
  scale_fill_manual(values=c("#74A089","#FDDDA0","#0B775E", "#FAD510" ,"#F2300F","#C93312"))+
  guides(ncol=4)+
  scale_x_continuous(name="Incubation time",breaks=c(-10,0,10,20,30),labels=c("not\ndigested","0 min","10 min","20 min", "30 min"))+
  theme(axis.text.x= element_text(colour="#000000", size = font_size, family = font_family),
        axis.text.y = element_text(colour ="#000000", size = font_size, family = font_family),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour="#000000", size = font_size, family = font_family),
        legend.title= element_text(colour="#000000",face="bold", size = font_size, family = font_family),
        legend.text = element_text(colour="#000000", size = font_size, family = font_family),
        legend.position="none",
        legend.box.background = element_rect(colour="black"),
        panel.background = element_rect(colour="#000000", fill=NA),
        plot.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour="#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  ggtitle(label="Tetrasaccharides")

Fig_samples_filtered_DP5_nolog10<-ggplot()+
  geom_point(data=wo_stds_filtered[wo_stds_filtered$DP==5,],aes(x=Position,y=Intensity_IM,fill=RTandDT),pch=21, size=3)+
  scale_y_continuous(name="Intensity (A.U.)",limits=c(300,5000000))+
  scale_fill_manual(values=c("#74A089","#FDDDA0","#0B775E", "#FAEFD1","#F2300F", "#899DA4",  "#FAD510","#C93312"))+
  scale_x_continuous(name="Incubation time",breaks=c(-10,0,10,20,30),labels=c("not\ndigested","0 min","10 min","20 min", "30 min"))+
  theme(axis.text.x= element_text(colour="#000000", size = font_size, family = font_family),
        axis.text.y = element_text(colour ="#000000", size = font_size, family = font_family),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour="#000000", size = font_size, family = font_family),
        legend.title= element_text(colour="#000000",face="bold", size = font_size, family = font_family),
        legend.text = element_text(colour="#000000", size = font_size, family = font_family),
        legend.position="none",
        legend.box.background = element_rect(colour="black"),
        panel.background = element_rect(colour="#000000", fill=NA),
        plot.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour="#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  ggtitle(label = "Pentasaccharides")



DP7_diversification_nolog10<-plot_grid(Fig_samples_filtered_DP3_nolog10,Fig_samples_filtered_DP4_nolog10,Fig_samples_filtered_DP5_nolog10, align="h",nrow=1,rel_widths = c(1.05,1,1))
save_plot(filename = "C:/Users/hbuck/Dropbox/Deep-water kelp/Nantes2019/Fig_DP7_diversification_nolog10.tiff",DP7_diversification_nolog10,base_height = 5,base_width = 12)          

