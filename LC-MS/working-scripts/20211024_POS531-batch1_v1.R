setwd("/Users/margotbligh/ownCloud/marglyco/Data/LC-MS_data/POS531/MS31_20211020")

#1: Install packages 
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

#functions
scientific_function <- function(x) {
    text <- gsub("E0", "", gsub("e\\+0", "E", scales::scientific_format()(x)))
    text
}


#2. Import and inspect MS data --------------------------------------------------------
#get file paths to mzML files
fp <- dir(path = "mzML", all.files = FALSE, full.names = TRUE)

#create phenodata data.frame
pd <- data.frame(name = basename(fp) %>%
                     gsub("MS31_20211020_|.mzML", "", .),
                 sampletype = basename(fp) %>%
                     sub(".*POS531.*", "extract", .) %>% 
                     sub(".*procedure.*", "procedureblank", .) %>% 
                     sub(".*QCmix.*", "QCmix", .) %>% 
                     sub(".*Qcspike.*", "QCspikemix", .),
                 station = basename(fp) %>%
                     gsub(".*POS531_|_Ex_.*", "", .) %>% 
                     sub("MS31_.*", "NA", .),
                 timepoint = basename(fp) %>%
                     gsub(".*Ex_|[ABCDEF]_\\d{2}.mzML", "", .) %>% 
                     sub("MS31_.*", "NA", .),
                 group = basename(fp) %>% 
                     gsub("MS31_20211020_|[ABCDEF]_\\d{2}.mzML", "", .) %>% 
                     sub(".*procedure.*", "procedureblank", .) %>% 
                     sub(".*QCmix.*", "QCmix", .) %>% 
                     sub(".*Qcspike.*", "QCspikemix", .),
                 replicate = basename(fp) %>%
                     gsub(".*Ex_\\d{1,2}|_\\d{2}.mzML", "", .) %>% 
                     sub(".*mix|.*spike-|.*blank", "", .),
                 run = basename(fp) %>% 
                     gsub(".*Ex_\\d{1,2}[ABCDEF]_|.mzML", "",.) %>% 
                     sub(".*[1-9]_|.*spike-\\d{2,3}_", "", .) %>% 
                     sub(".*[[:alpha:]]_", "", .) %>% as.numeric(),
                 file = fp,
                 stringsAsFactors = FALSE)
pd <- pd[pd$run <=51,]
fp <- fp[fp %in% pd$file]
pd$file <- NULL

#read in data
all_data <- readMSData(files = fp, 
                       pdata = new("NAnnotatedDataFrame", 
                                   pd), 
                       mode = "onDisk")

#split MS1 and MS2
data <- all_data[all_data@featureData@data$msLevel == 1]
data_ms2 <- all_data

#3: Plot TIC ----
pal_group <- hcl.colors(n = length(unique(pd$group)),
                        palette = "Dark3")
names(pal_group) <- unique(pd$group)

#plot the tic as boxplot
tc <- split(tic(all_data), 
            f = fromFile(all_data))
boxplot(tc, 
        col = pal_group[all_data$group],
        ylab = "intensity", 
        main = "total ion current",
        names = all_data$name,
        las=2,
        cex.axis = 0.8)

#4: Extract chromatograms for stds -----
stds.masses <- fread("standards-masses.csv")
stds.masses <- stds.masses[complete.cases(stds.masses),]
stds.mz <- c(stds.masses$`Reduced [M+H]+ ion m/z`, stds.masses$`[M+H]+ ion m/z`)
names(stds.mz) <- c(stds.masses$Compound, stds.masses$Compound)
stds.mz <- unique(stds.mz)

#hexose monomer----
hexose.chr <- chromatogram(data,
                           mz = c(400.244-0.001, 400.244+0.001),
                           rt = c(1*60, 15*60))
red_hexose.chr <- chromatogram(data,
                               mz = c(398.229-0.001, 398.229+0.001),
                               rt = c(1*60, 15*60))
chr_hex.df <- data.frame(sample = as.character(),
                         group = as.character(),
                         rt = as.numeric(),
                         intensity = as.numeric(), 
                         reduced = as.character())

for (i in 1:length(pd$name)){
    rt = hexose.chr[[i]]@rtime/60
    intensity = hexose.chr[[i]]@intensity
    sample = rep(hexose.chr$name[i], length(rt))
    group = rep(hexose.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "n")
    chr_hex.df <- rbind(chr_hex.df,
                        temp)
}
for (i in 1:length(pd$name)){
    rt = red_hexose.chr[[i]]@rtime/60
    intensity = red_hexose.chr[[i]]@intensity
    sample = rep(red_hexose.chr$name[i], length(rt))
    group = rep(red_hexose.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "y")
    chr_hex.df <- rbind(chr_hex.df,
                        temp)
}

chr_hex.df[is.na(chr_hex.df)] <- 0


ggplot() +
    geom_line(mapping = aes(rt, intensity, colour = reduced, group = sample),
              data = chr_hex.df,
              lwd = 1) +
    scale_colour_manual(values = c("#FEC000", "#2BB6AF"), name = "Reduced") +
    labs(x= "Retention time (min)", y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group), scales = "free_y") +
    xlim(0, 15) +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, angle = 360, hjust = 0,
                                      family = "Avenir LT 65 Medium"),
          strip.background = element_blank(),
          axis.text = element_text(size = 12, family = "Avenir"),
          axis.title = element_text(size = 12, family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          legend.position = "top",
          legend.text = element_text(size = 12, family = "Avenir"),
          axis.line = element_blank())





#pentose monomer----
pentose.chr <- chromatogram(data,
                           mz = c(370.234-0.001, 370.234+0.001),
                           rt = c(1*60, 15*60))
red_pentose.chr <- chromatogram(data,
                               mz = c(368.218-0.001, 368.218+0.001),
                               rt = c(1*60, 15*60))
chr_pent.df <- data.frame(sample = as.character(),
                         group = as.character(),
                         rt = as.numeric(),
                         intensity = as.numeric(), 
                         reduced = as.character())

for (i in 1:length(pd$name)){
    rt = pentose.chr[[i]]@rtime/60
    intensity = pentose.chr[[i]]@intensity
    sample = rep(pentose.chr$name[i], length(rt))
    group = rep(pentose.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "n")
    chr_pent.df <- rbind(chr_pent.df,
                        temp)
}
for (i in 1:length(pd$name)){
    rt = red_pentose.chr[[i]]@rtime/60
    intensity = red_pentose.chr[[i]]@intensity
    sample = rep(red_pentose.chr$name[i], length(rt))
    group = rep(red_pentose.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "y")
    chr_pent.df <- rbind(chr_pent.df,
                        temp)
}

chr_pent.df[is.na(chr_pent.df)] <- 0


ggplot() +
    geom_line(mapping = aes(rt, intensity, colour = reduced, group = sample),
              data = chr_pent.df,
              lwd = 1) +
    scale_colour_manual(values = c("#FEC000", "#2BB6AF"), name = "Reduced") +
    labs(x= "Retention time (min)", y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group), scales = "free_y") +
    xlim(0, 15) +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, angle = 360, hjust = 0,
                                      family = "Avenir LT 65 Medium"),
          strip.background = element_blank(),
          axis.text = element_text(size = 12, family = "Avenir"),
          axis.title = element_text(size = 12, family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          legend.position = "top",
          legend.text = element_text(size = 12, family = "Avenir"),
          axis.line = element_blank())





#hexose monomer----
hexose.chr <- chromatogram(data,
                           mz = c(400.244-0.001, 400.244+0.001),
                           rt = c(1*60, 15*60))
red_hexose.chr <- chromatogram(data,
                               mz = c(398.229-0.001, 398.229+0.001),
                               rt = c(1*60, 15*60))
chr_hex.df <- data.frame(sample = as.character(),
                         group = as.character(),
                         rt = as.numeric(),
                         intensity = as.numeric(), 
                         reduced = as.character())

for (i in 1:length(pd$name)){
    rt = hexose.chr[[i]]@rtime/60
    intensity = hexose.chr[[i]]@intensity
    sample = rep(hexose.chr$name[i], length(rt))
    group = rep(hexose.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "n")
    chr_hex.df <- rbind(chr_hex.df,
                        temp)
}
for (i in 1:length(pd$name)){
    rt = red_hexose.chr[[i]]@rtime/60
    intensity = red_hexose.chr[[i]]@intensity
    sample = rep(red_hexose.chr$name[i], length(rt))
    group = rep(red_hexose.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "y")
    chr_hex.df <- rbind(chr_hex.df,
                        temp)
}

chr_hex.df[is.na(chr_hex.df)] <- 0


ggplot() +
    geom_line(mapping = aes(rt, intensity, colour = reduced, group = sample),
              data = chr_hex.df,
              lwd = 1) +
    scale_colour_manual(values = c("#FEC000", "#2BB6AF"), name = "Reduced") +
    labs(x= "Retention time (min)", y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group), scales = "free_y") +
    xlim(0, 15) +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, angle = 360, hjust = 0,
                                      family = "Avenir LT 65 Medium"),
          strip.background = element_blank(),
          axis.text = element_text(size = 12, family = "Avenir"),
          axis.title = element_text(size = 12, family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          legend.position = "top",
          legend.text = element_text(size = 12, family = "Avenir"),
          axis.line = element_blank())





#hexose dimer----
hexose2.chr <- chromatogram(data,
                           mz = c(562.298-0.001, 562.298+0.001),
                           rt = c(10*60, 25*60))
red_hexose2.chr <- chromatogram(data,
                               mz = c(560.282-0.001, 560.282+0.001),
                               rt = c(10*60, 25*60))
chr_hex2.df <- data.frame(sample = as.character(),
                         group = as.character(),
                         rt = as.numeric(),
                         intensity = as.numeric(), 
                         reduced = as.character())

for (i in 1:length(pd$name)){
    rt = hexose2.chr[[i]]@rtime/60
    intensity = hexose2.chr[[i]]@intensity
    sample = rep(hexose2.chr$name[i], length(rt))
    group = rep(hexose2.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "n")
    chr_hex2.df <- rbind(chr_hex2.df,
                        temp)
}
for (i in 1:length(pd$name)){
    rt = red_hexose2.chr[[i]]@rtime/60
    intensity = red_hexose2.chr[[i]]@intensity
    sample = rep(red_hexose2.chr$name[i], length(rt))
    group = rep(red_hexose2.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "y")
    chr_hex2.df <- rbind(chr_hex2.df,
                        temp)
}

chr_hex2.df[is.na(chr_hex2.df)] <- 0


ggplot() +
    geom_line(mapping = aes(rt, intensity, colour = reduced, group = sample),
              data = chr_hex2.df,
              lwd = 1) +
    scale_colour_manual(values = c("#FEC000", "#2BB6AF"), name = "Reduced") +
    labs(x= "Retention time (min)", y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group), scales = "free_y") +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, angle = 360, hjust = 0,
                                      family = "Avenir LT 65 Medium"),
          strip.background = element_blank(),
          axis.text = element_text(size = 12, family = "Avenir"),
          axis.title = element_text(size = 12, family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          legend.position = "top",
          legend.text = element_text(size = 12, family = "Avenir"),
          axis.line = element_blank())






#pentose dimer----
pentose2.chr <- chromatogram(data,
                            mz = c(502.27645-0.001, 502.27645+0.001),
                            rt = c(1*60, 15*60))
red_pentose2.chr <- chromatogram(data,
                                mz = c(500.2608-0.001, 500.2608+0.001),
                                rt = c(1*60, 15*60))
chr_pent2.df <- data.frame(sample = as.character(),
                          group = as.character(),
                          rt = as.numeric(),
                          intensity = as.numeric(), 
                          reduced = as.character())

for (i in 1:length(pd$name)){
    rt = pentose2.chr[[i]]@rtime/60
    intensity = pentose2.chr[[i]]@intensity
    sample = rep(pentose2.chr$name[i], length(rt))
    group = rep(pentose2.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "n")
    chr_pent2.df <- rbind(chr_pent2.df,
                         temp)
}
for (i in 1:length(pd$name)){
    rt = red_pentose2.chr[[i]]@rtime/60
    intensity = red_pentose2.chr[[i]]@intensity
    sample = rep(red_pentose2.chr$name[i], length(rt))
    group = rep(red_pentose2.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "y")
    chr_pent2.df <- rbind(chr_pent2.df,
                         temp)
}

chr_pent2.df[is.na(chr_pent2.df)] <- 0


ggplot() +
    geom_line(mapping = aes(rt, intensity, colour = reduced, group = sample),
              data = chr_pent2.df,
              lwd = 1) +
    scale_colour_manual(values = c("#FEC000", "#2BB6AF"), name = "Reduced") +
    labs(x= "Retention time (min)", y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group), scales = "free_y") +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, angle = 360, hjust = 0,
                                      family = "Avenir LT 65 Medium"),
          strip.background = element_blank(),
          axis.text = element_text(size = 12, family = "Avenir"),
          axis.title = element_text(size = 12, family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          legend.position = "top",
          legend.text = element_text(size = 12, family = "Avenir"),
          axis.line = element_blank())





#hexose deoxyhexose dimer----
hexose2deoxy1.chr <- chromatogram(data,
                            mz = c(546.30267-0.001, 546.30267+0.001),
                            rt = c(7*60, 15*60))
red_hexose2deoxy1.chr <- chromatogram(data,
                                mz = c(544.28702-0.001, 544.28702+0.001),
                                rt = c(7*60, 15*60))
chr_hex2deoxy1.df <- data.frame(sample = as.character(),
                          group = as.character(),
                          rt = as.numeric(),
                          intensity = as.numeric(), 
                          reduced = as.character())

for (i in 1:length(pd$name)){
    rt = hexose2deoxy1.chr[[i]]@rtime/60
    intensity = hexose2deoxy1.chr[[i]]@intensity
    sample = rep(hexose2deoxy1.chr$name[i], length(rt))
    group = rep(hexose2deoxy1.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "n")
    chr_hex2deoxy1.df <- rbind(chr_hex2deoxy1.df,
                         temp)
}
for (i in 1:length(pd$name)){
    rt = red_hexose2deoxy1.chr[[i]]@rtime/60
    intensity = red_hexose2deoxy1.chr[[i]]@intensity
    sample = rep(red_hexose2deoxy1.chr$name[i], length(rt))
    group = rep(red_hexose2deoxy1.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "y")
    chr_hex2deoxy1.df <- rbind(chr_hex2deoxy1.df,
                         temp)
}

chr_hex2deoxy1.df[is.na(chr_hex2deoxy1.df)] <- 0


ggplot() +
    geom_line(mapping = aes(rt, intensity, colour = reduced, group = sample),
              data = chr_hex2deoxy1.df,
              lwd = 1) +
    scale_colour_manual(values = c("#FEC000", "#2BB6AF"), name = "Reduced") +
    labs(x= "Retention time (min)", y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group), scales = "free_y") +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, angle = 360, hjust = 0,
                                      family = "Avenir LT 65 Medium"),
          strip.background = element_blank(),
          axis.text = element_text(size = 12, family = "Avenir"),
          axis.title = element_text(size = 12, family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          legend.position = "top",
          legend.text = element_text(size = 12, family = "Avenir"),
          axis.line = element_blank())


#hexose pentose dimer-----
hex1pent1.chr <- chromatogram(data,
                             mz = c(532.28702-0.001, 532.28702+0.001),
                             rt = c(9*60, 17*60))
red_hex1pent1.chr <- chromatogram(data,
                                 mz = c(530.27137-0.001, 530.27137+0.001),
                                 rt = c(9*60, 17*60))
chr_hex1pent1.df <- data.frame(sample = as.character(),
                           group = as.character(),
                           rt = as.numeric(),
                           intensity = as.numeric(), 
                           reduced = as.character())

for (i in 1:length(pd$name)){
    rt = hex1pent1.chr[[i]]@rtime/60
    intensity = hex1pent1.chr[[i]]@intensity
    sample = rep(hex1pent1.chr$name[i], length(rt))
    group = rep(hex1pent1.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "n")
    chr_hex1pent1.df <- rbind(chr_hex1pent1.df,
                          temp)
}
for (i in 1:length(pd$name)){
    rt = red_hex1pent1.chr[[i]]@rtime/60
    intensity = red_hex1pent1.chr[[i]]@intensity
    sample = rep(red_hex1pent1.chr$name[i], length(rt))
    group = rep(red_hex1pent1.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "y")
    chr_hex1pent1.df <- rbind(chr_hex1pent1.df,
                          temp)
}

chr_hex1pent1.df[is.na(chr_hex1pent1.df)] <- 0


ggplot() +
    geom_line(mapping = aes(rt, intensity, colour = reduced, group = sample),
              data = chr_hex1pent1.df,
              lwd = 1) +
    scale_colour_manual(values = c("#FEC000", "#2BB6AF"), name = "Reduced") +
    labs(x= "Retention time (min)", y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group), scales = "free_y") +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, angle = 360, hjust = 0,
                                      family = "Avenir LT 65 Medium"),
          strip.background = element_blank(),
          axis.text = element_text(size = 12, family = "Avenir"),
          axis.title = element_text(size = 12, family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          legend.position = "top",
          legend.text = element_text(size = 12, family = "Avenir"),
          axis.line = element_blank())



#hexose trimer----
hexose3.chr <- chromatogram(data,
                            mz = c(724.350-0.001, 724.350+0.001),
                            rt = c(15*60, 30*60))
red_hexose3.chr <- chromatogram(data,
                                mz = c(722.335-0.001, 722.335+0.001),
                                rt = c(15*60, 30*60))
chr_hex3.df <- data.frame(sample = as.character(),
                          group = as.character(),
                          rt = as.numeric(),
                          intensity = as.numeric(), 
                          reduced = as.character())

for (i in 1:length(pd$name)){
    rt = hexose3.chr[[i]]@rtime/60
    intensity = hexose3.chr[[i]]@intensity
    sample = rep(hexose3.chr$name[i], length(rt))
    group = rep(hexose3.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "n")
    chr_hex3.df <- rbind(chr_hex3.df,
                         temp)
}
for (i in 1:length(pd$name)){
    rt = red_hexose3.chr[[i]]@rtime/60
    intensity = red_hexose3.chr[[i]]@intensity
    sample = rep(red_hexose3.chr$name[i], length(rt))
    group = rep(red_hexose3.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "y")
    chr_hex3.df <- rbind(chr_hex3.df,
                         temp)
}

chr_hex3.df[is.na(chr_hex3.df)] <- 0


ggplot() +
    geom_line(mapping = aes(rt, intensity, colour = reduced, group = sample),
              data = chr_hex3.df,
              lwd = 1) +
    scale_colour_manual(values = c("#FEC000", "#2BB6AF"), name = "Reduced") +
    labs(x= "Retention time (min)", y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group), scales = "free_y") +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, angle = 360, hjust = 0,
                                      family = "Avenir LT 65 Medium"),
          strip.background = element_blank(),
          axis.text = element_text(size = 12, family = "Avenir"),
          axis.title = element_text(size = 12, family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          legend.position = "top",
          legend.text = element_text(size = 12, family = "Avenir"),
          axis.line = element_blank())






#pentose trimer----
pentose3.chr <- chromatogram(data,
                             mz = c(634.31871-0.001, 634.31871+0.001),
                             rt = c(10*60, 20*60))
red_pentose3.chr <- chromatogram(data,
                                 mz = c(632.30306-0.001, 632.30306+0.001),
                                 rt = c(10*60, 20*60))
chr_pent3.df <- data.frame(sample = as.character(),
                           group = as.character(),
                           rt = as.numeric(),
                           intensity = as.numeric(), 
                           reduced = as.character())

for (i in 1:length(pd$name)){
    rt = pentose3.chr[[i]]@rtime/60
    intensity = pentose3.chr[[i]]@intensity
    sample = rep(pentose3.chr$name[i], length(rt))
    group = rep(pentose3.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "n")
    chr_pent3.df <- rbind(chr_pent3.df,
                          temp)
}
for (i in 1:length(pd$name)){
    rt = red_pentose3.chr[[i]]@rtime/60
    intensity = red_pentose3.chr[[i]]@intensity
    sample = rep(red_pentose3.chr$name[i], length(rt))
    group = rep(red_pentose3.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "y")
    chr_pent3.df <- rbind(chr_pent3.df,
                          temp)
}

chr_pent3.df[is.na(chr_pent3.df)] <- 0


ggplot() +
    geom_line(mapping = aes(rt, intensity, colour = reduced, group = sample),
              data = chr_pent3.df,
              lwd = 1) +
    scale_colour_manual(values = c("#FEC000", "#2BB6AF"), name = "Reduced") +
    labs(x= "Retention time (min)", y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group), scales = "free_y") +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, angle = 360, hjust = 0,
                                      family = "Avenir LT 65 Medium"),
          strip.background = element_blank(),
          axis.text = element_text(size = 12, family = "Avenir"),
          axis.title = element_text(size = 12, family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          legend.position = "top",
          legend.text = element_text(size = 12, family = "Avenir"),
          axis.line = element_blank())






#hexose tetramer----
hexose4.chr <- chromatogram(data,
                            mz = c(886.40324-0.001, 886.40324+0.001),
                            rt = c(20*60, 35*60))
red_hexose4.chr <- chromatogram(data,
                                mz = c(884.38759-0.001, 884.38759+0.001),
                                rt = c(20*60, 35*60))
chr_hex4.df <- data.frame(sample = as.character(),
                          group = as.character(),
                          rt = as.numeric(),
                          intensity = as.numeric(), 
                          reduced = as.character())

for (i in 1:length(pd$name)){
    rt = hexose4.chr[[i]]@rtime/60
    intensity = hexose4.chr[[i]]@intensity
    sample = rep(hexose4.chr$name[i], length(rt))
    group = rep(hexose4.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "n")
    chr_hex4.df <- rbind(chr_hex4.df,
                         temp)
}
for (i in 1:length(pd$name)){
    rt = red_hexose4.chr[[i]]@rtime/60
    intensity = red_hexose4.chr[[i]]@intensity
    sample = rep(red_hexose4.chr$name[i], length(rt))
    group = rep(red_hexose4.chr$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity,
                       reduced = "y")
    chr_hex4.df <- rbind(chr_hex4.df,
                         temp)
}

chr_hex4.df[is.na(chr_hex4.df)] <- 0


ggplot() +
    geom_line(mapping = aes(rt, intensity, colour = reduced, group = sample),
              data = chr_hex4.df,
              lwd = 1) +
    scale_colour_manual(values = c("#FEC000", "#2BB6AF"), name = "Reduced") +
    labs(x= "Retention time (min)", y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group), scales = "free_y") +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, angle = 360, hjust = 0,
                                      family = "Avenir LT 65 Medium"),
          strip.background = element_blank(),
          axis.text = element_text(size = 12, family = "Avenir"),
          axis.title = element_text(size = 12, family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          legend.position = "top",
          legend.text = element_text(size = 12, family = "Avenir"),
          axis.line = element_blank())





