setwd("/Users/margotbligh/Google_Drive/MPI_PhD/Lab-things/alpha-mannan/Polaribacter_Hel1-33-78_enzymes/202109_GH99_PLx_GH92")

#1: Install packages --------------------------------------------------------
library(BiocStyle)
library(xcms)
library(ggplot2)
library(tidyverse)
library(scales)
library(data.table)
library(MSnbase)
library(CAMERA)
library(plyr)
library(viridis)
library(patchwork)
library(grid)
library(systemfonts)
library(wesanderson)
library(ggsci)

#2. Import and inspect MS data --------------------------------------------------------
#get file paths to mzML files
fp <- dir(path = "direct-injection", 
          all.files = FALSE, 
          full.names = TRUE)

#create phenodata data.frame
pd <- data.frame(name = basename(fp) %>% 
                     sub(".*313.*", "k-i-k 6", .) %>% 
                     sub(".*391.*", "k 6", .) %>% 
                     sub("20210907_DI1_neg_", "", .) %>% 
                     gsub("_1000xdilute|.mzML", "", .),
                 ion = basename(fp) %>% 
                     sub(".*313.*", 313, .) %>% 
                     sub(".*391.*", 391, .) %>% 
                     sub("20210907_DI1_neg_.*", "NA", .),
                 stringsAsFactors = FALSE)

#read in data
data <- readMSData(files = fp, 
                   pdata = new("NAnnotatedDataFrame", 
                               pd), 
                   mode = "onDisk")

#3. Plot MS2 for 313 ion -----
kik6.m4h <- 313.0234948299089

#subset data
data.ms2 <- data[data@featureData@data$msLevel == 2]
data_313.ms2 <- filterFile(data.ms2,
                           file =  data$ion == 313)

chr313 <- chromatogram(data_313.ms2,
                       mz = c(kik6.m4h - 0.001,
                              kik6.m4h + 0.001),
                       msLevel = 2) 
plot(chr313, xlim = c(146, 150))


data_313.ms2.rt <- filterRt(data_313.ms2, rt = c(146, 150))

#extract data for retention time 149.0215
mz <-data_313.ms2.rt[[14]]@mz
raw_intensity <- data_313.ms2.rt[[14]]@intensity
intensity <- raw_intensity / sum(raw_intensity) * 100
labels <- as.character(round(mz, 3))
labels[intensity < 1.8] <- ""


tiff(filename = "./analysis/ms2_plots/kik6_313_directinjection_v1.tiff",
     height = 6, width = 6, units = "in", res = 600)
ggplot() +
    geom_segment(aes(x= mz, 
                     xend= mz, 
                     y=0, 
                     yend=intensity),
                 lwd = 1) +
    geom_text(aes(x = mz,
                  y = intensity,
                  label = labels),
              angle = 90,
              nudge_y = 5,
              size = 5.5,
              family = "Avenir") +
    scale_x_continuous(breaks = seq(0,1000, by = 100),
                       name = expression(italic(m/z))) +
    scale_y_continuous(name = "Relative intensity (%)",
                       expand = expansion(mult = c(0, 0.2)),
                       breaks = seq(0, 100, by = 20)) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank())
dev.off()


