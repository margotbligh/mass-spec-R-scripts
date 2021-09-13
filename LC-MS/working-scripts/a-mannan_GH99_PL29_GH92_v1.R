setwd("/Users/margotbligh/Google_Drive/MPI_PhD/Lab-things/alpha-mannan/Polaribacter_Hel1-33-78_enzymes/202108_GH99_PLx_GH92")
load("./analysis/RData/RData_20210903.RData")


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
fp <- dir(path = "MS31_20210901", 
          all.files = FALSE, 
          full.names = TRUE)

#create phenodata data.frame
#each sample must have a unique name!
pd <- data.frame(name = basename(fp) %>%
                     sub(".*mannose-3-sulfate_1mg_perL_17.*", "mannose-3-sulfate_2", .) %>% 
                     sub(".*mannose-3-sulfate_1mg.*", "mannose-3-sulfate_1", .) %>% 
                     sub("MS31_20210901_", "", .) %>%
                     sub("_1mg.*", "", .) %>% 
                     sub("k_01.*", "k_1", .) %>% 
                     sub("k_03.*", "k_2", .) %>% 
                     sub("k_12.*", "k_3", .) %>% 
                     sub("_\\d\\d.mzML", "", .),
                 sample_type = basename(fp) %>% 
                     sub(".*_[2345][ab]_.*", "digest", .) %>% 
                     sub(".*_1[ab]_.*", "neg", .) %>% 
                     sub(".*_blank_.*", "solvent blank", .) %>% 
                     sub(".*mannose.*", "standard", .) %>% 
                     sub(".*HILIC.*", "HILIC standard", .),
                 enzyme = basename(fp) %>% 
                     sub(".*_2[ab]_.*", "PLx", .) %>% 
                     sub(".*_3[ab]_.*", "GH99-sulphatase", .) %>% 
                     sub(".*_4[ab]_.*", "GH92", .) %>% 
                     sub(".*_5[ab]_.*", "GH92 + GH99-sulphatase", .) %>%
                     sub("MS31.*", NA, .),
                 group = basename(fp) %>% 
                     sub(".*_1[ab]_.*", "no digest", .) %>% 
                     sub(".*_2[ab]_.*", "PLx digest", .) %>% 
                     sub(".*_3[ab]_.*", "GH99-sulphatase digest", .) %>% 
                     sub(".*_4[ab]_.*", "GH92 digest", .) %>% 
                     sub(".*_5[ab]_.*", "GH92 + GH99-sulphatase digest", .) %>% 
                     gsub("MS31_20210901_|_1mg.*", "", .) %>% 
                     sub(".*blank.*", "solvent blank", .) %>% 
                     sub(".*HILIC.*", "HILIC standard", .),
                 stringsAsFactors = FALSE)

#read in data
all_data <- readMSData(files = fp, 
                   pdata = new("NAnnotatedDataFrame", 
                               pd), 
                   mode = "onDisk")

#split MS1 and MS2
data <- all_data[all_data@featureData@data$msLevel == 1]
data_ms2 <- all_data
    
#3: Create initial output directories -------------------------------------
dir.create("./analysis",
           showWarnings = FALSE)
dir.create("./analysis/RData",
           showWarnings = FALSE)
dir.create("./analysis/analysis_plots",
           showWarnings = FALSE)
dir.create("./analysis/analysis_tables",
           showWarnings = FALSE)
dir.create("./analysis/ms2_plots",
           showWarnings = FALSE)
dir.create("./analysis/processing_plots",
           showWarnings = FALSE)

#save RData object
save(data, 
     file = "./analysis/RData/data.RData")
save(data_ms2, 
     file = "./analysis/RData/data_ms2.RData")


#4: Plot TIC ----
dir.create("./analysis/processing_plots/tic",
           showWarnings = FALSE)

pal_group <- hcl.colors(n = length(unique(pd$group)),
                        palette = "Dark3")
names(pal_group) <- unique(pd$group)

#plot the tic as boxplot
tc <- split(tic(all_data), 
            f = fromFile(all_data))
cairo_pdf("./analysis/processing_plots/tic/tic_boxplot.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mar=c(9,5,1,1))
boxplot(tc, 
        col = pal_group[all_data$group],
        ylab = "intensity", 
        main = "total ion current",
        names = all_data$name,
        las=2,
        cex.axis = 0.8)
dev.off()

#plot as a chromatogram
tic <- chromatogram(all_data, aggregationFun = "sum")
tic.df <- data.frame(sample = as.character(),
                     group = as.character(),
                     rt = as.numeric(),
                     intensity = as.numeric())
for (i in 1:length(pd$name)){
    rt = tic[[i]]@rtime / 60
    intensity = tic[[i]]@intensity
    sample = rep(all_data$name[i], length(rt))
    group = rep(all_data$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity)
    tic.df <- rbind(temp,
                    tic.df)
    
}

cairo_pdf("./analysis/processing_plots/tic/tic_chromatogram.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
ggplot() +
    geom_line(mapping = aes(rt,
                            intensity,
                            colour = group,
                            group = sample),
              data = tic.df,
              lwd = 1.2) +
    scale_colour_manual(values = pal_group) +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group),
               scales = "free_y") +
    theme_classic() +
    theme(strip.text = element_blank(),
          axis.text = element_text(size = 12,
                                   family = "Avenir"),
          axis.title = element_text(size = 14, 
                                    family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          legend.text = element_text(size = 6, 
                                     family = "Avenir"),
          axis.line = element_blank())
dev.off()

#5: Peak picking (CentWave) ---------------------------
#set parameters
cwp<-CentWaveParam()
cwp@ppm<-2.5
cwp@peakwidth<-c(10,45)
cwp@snthresh<-20
cwp@noise <- 5000
cwp@prefilter <- c(5, 1000)
cwp@mzdiff <- 0.001

#check with chromatograms for sulphated mannoase
dir.create("./analysis/processing_plots/peakpicking",
           showWarnings = FALSE)

error = 0.0025
chr1 <- chromatogram(data,
                     rt = c(50, 400),
                     mz = c(259.01293 - error,
                            259.01293 + error))

pal1 <- hcl.colors(n = length(pd$name),
                   palette = "Dark3")
names(pal1) <- pd$name

chr1_cwp <- findChromPeaks(chr1, param = cwp)

pal1.std <- pal1[which(pd$sample_type == "standard")]
cairo_pdf("./analysis/processing_plots/peakpicking/peakpicking_2.2ppm_20sn_pw10to45_noise5000_prefilter3-1000_mzdiff0.001_sulphatedmannose-standards.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mfrow = c(3,2))
for (i in 1:length(pal1.std)){
    plot(chr1_cwp[,pd$name == names(pal1.std)[i]],
         lwd = 2,
         cex.main = 1,
         peakCol = "black",
         peakType = "rectangle",
         main = names(pal1.std)[i])
}
dev.off()

pal1.digest <- pal1[which(pd$sample_type == "digest")]
cairo_pdf("./analysis/processing_plots/peakpicking/peakpicking_2.2ppm_20sn_pw10to45_noise5000_prefilter3-1000_mzdiff0.001_sulphatedmannose-digests.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mfrow = c(4,2))
for (i in 1:length(pal1.digest)){
    plot(chr1_cwp[,pd$name == names(pal1.digest)[i]],
         lwd = 2,
         cex.main = 1,
         peakCol = "black",
         peakType = "rectangle",
         main = names(pal1.digest)[i])
}
dev.off()

#pick peaks
data_peaks<-findChromPeaks(data, 
                           param=cwp)


#check picked peaks match chromatograms
error = 0.0025
chr2 <- chromatogram(data_peaks,
                     rt = c(50, 400),
                     mz = c(259.01293 - error,
                            259.01293 + error))

cairo_pdf("./analysis/processing_plots/peakpicking/peakpicking_2.5ppm_20sn_pw7to45_noise5000_prefilter3-1000_mzdiff0.001_sulphatedmannose-standards_picked.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mfrow = c(3,2))
for (i in 1:length(pal1.std)){
    plot(chr2[,pd$name == names(pal1.std)[i]],
         lwd = 2,
         cex.main = 1,
         peakCol = "black",
         peakType = "rectangle",
         main = names(pal1.std)[i])
}
dev.off()

cairo_pdf("./analysis/processing_plots/peakpicking/peakpicking_2.5ppm_20sn_pw7to45_noise5000_prefilter3-1000_mzdiff0.001_sulphatedmannose-digests_picked.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mfrow = c(4,2))
for (i in 1:length(pal1.digest)){
    plot(chr2[,pd$name == names(pal1.digest)[i]],
         lwd = 2,
         cex.main = 1,
         peakCol = "black",
         peakType = "rectangle",
         main = names(pal1.digest)[i])
}
dev.off()

##1.3 ppm; peakwidth 7-45; snthresh 20; noise 5000; prefilter 3, 1000
#on chromatogram:
    #standards: everything picked + additional noise
    #digests: everything picked - 5a only half of peak picked
#on data:
    #standards: smaller earlier peaks missed for mannose-3-sulphate, less noise 
    #digests: nothing picked!?

##ppm changed to 2
#on chromatogram:
    #standards: everything picked + additional noise (same as 1.3)
    #digests: everything picked - 5a only half of peak picked (same as 1.3)
#on data:
    #standards: smaller earlier peaks missed for mannose-3-sulphate,
            #less noise, smaller early peak for mannose-6-sulphate split
    #digests: peaks picked for 4b, 5a and 5b as for chromatogram, 4a missing

##ppm changed to 2.5
#on chromatogram:
    #standards: everything picked + additional noise (same as 1.3)
    #digests: everything picked - 5a only half of peak picked (same as 1.3)
#on data:
    #standards: smaller earlier peaks missed for mannose-3-sulphate,
        #noisy, smaller early peak for mannose-6-sulphate split
    #digests: peaks picked for all, but 4a peak split into 3

##ppm changed to 2.2, peak width 10 to 45, noise 2000, prefilter 5 1000
#on chromatogram:
    #standards: everything picked + additional noise (same as 1.3)
    #digests: everything picked nicely!
#on data:
    #standards: smaller earlier peaks missed for mannose-3-sulphate,
        #mannose-2-sulphate missed entirely!
        #noisy, smaller early peak for mannose-6-sulphate split
    #digests: peaks picked for all, but 4a peak split into 2

##changed mzdiff to 0.001
#only thing changed is that now 4a peak not split into 2 :) 

##changed ppm to 1.3 again
#mannose-2-sulphate still not picked in data
#4a and 4b now not picked in data

##changed ppm to 2 again
#mannose-2-sulphate still not picked in data
#4a not picked in data

##changed ppm to 2.2 again
#mannose-2-sulphate still not picked in data
#all digests nicely picked :) 

##changed ppm to 2.5 again
#everything picked so nicely!!!!

#pick peaks for ms2
data_ms2<-findChromPeaks(data_ms2, 
                         param=cwp)

#save as RData objects
save(data_peaks, 
     file = "./analysis/RData/data_peaks.RData")
save(data_ms2, 
     file = "./analysis/RData/data_ms2.RData")

#6: Group peaks to create "features"---------
#parameters
pdp <- PeakDensityParam(sampleGroups = data$group,
                        binSize = 0.005,
                        bw = 6,
                        minFraction = 0.5) 

#check parameters
#extract and plot chromatograms to test settings
dir.create("./analysis/processing_plots/peakgrouping",
           showWarnings = FALSE)

cairo_pdf("./analysis/processing_plots/peakgrouping/peakgrouping_binsize0.005_minFraction0.5_bw6_sulphatedmannose.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mar=c(4,4,3,11))
plotChromPeakDensity(chr2, 
                     col = pal1, 
                     param = pdp,
                     peakBg = pal1[chromPeaks(chr2)[, "sample"]],
                     peakCol = pal1[chromPeaks(chr2)[, "sample"]],
                     peakPch = 16)
legend("topright",
       legend = paste0(seq(1,length(pal1),1),
                       "=",
                       names(pal1)),
       inset=c(-0.5,0),
       fill = pal1,
       pt.cex = 0.3,
       cex = 0.5,
       bty = "n",
       horiz = FALSE,
       xpd=TRUE,
       ncol = 2)
dev.off()


#group peaks
data_peaks_grouped <- groupChromPeaks(data_peaks, param = pdp)

save(data_peaks_grouped, 
     file = "./analysis/RData/data_peaks_grouped.RData")

#7: Fill in missing peaks----------
fpp <- FillChromPeaksParam()
data_peaks_grouped_filled <- fillChromPeaks(data_peaks_grouped)

save(data_peaks_grouped_filled, 
     file = "./analysis/RData/data_peaks_grouped_filled.RData")

#9: PCA plot: filled vs not filled features ----
dir.create("./analysis/processing_plots/pca",
           showWarnings = FALSE)

##NOT FILLED
#get intensity values
ft_ints <- featureValues(data_peaks_grouped)
ft_ints <- as.data.frame(ft_ints)
names(ft_ints) <- pd$name
#log transform
ft_ints <- log2(ft_ints)
#replace all NA values with zero
ft_ints[is.na(ft_ints)] <- 0
#perform PCA with the intensities mean centered
pc1 <- prcomp(t(ft_ints), center = TRUE)
#plot
tiff("./analysis/processing_plots/pca/pca-unfilled-featureIntensities.tiff",
     res = 300,
     height = 9,
     width = 12,
     units = "in")
pcSummary1 <- summary(pc1)
par(family = "Avenir",
    mar=c(5,4,2,9))
plot(pc1$x[, 1], 
     pc1$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary1$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary1$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "black", 
     lwd = 1.3,
     bg = pal_group[pd$group], 
     cex = 2)
grid()
legend("topleft",
       legend = names(pal_group) %>% unique(),
       pt.bg = pal_group %>% unique(),
       col = "black",
       pch = 21,
       pt.lwd = 1.3,
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       ncol = 1,
       x.intersp = 0.8,
       y.intersp = 1.2,
       inset=c(1,0),
       xpd=TRUE)
dev.off()

#seperate into digests (including neg control) and standards
#HILIC standard in the middle 
#solvent blanks in middle
#no separation of replicates


##FILLED
#get intensity values
ft_ints.fl <- featureValues(data_peaks_grouped_filled)
ft_ints.fl <- as.data.frame(ft_ints.fl)
names(ft_ints.fl) <- pd$name
#log transform
ft_ints.fl <- log2(ft_ints.fl)
#replace all NA values with zero
ft_ints.fl[is.na(ft_ints.fl)] <- 0
#perform PCA with the intensities mean centered
pc2 <- prcomp(t(ft_ints.fl), center = TRUE)
#plot
tiff("./analysis/processing_plots/pca/pca-unfilled-featureIntensities.tiff",
     res = 300,
     height = 9,
     width = 12,
     units = "in")
pcSummary2 <- summary(pc2)
par(family = "Avenir",
    mar=c(5,4,2,9))
plot(pc2$x[, 1], 
     pc2$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary2$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary2$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "black", 
     lwd = 1.3,
     bg = pal_group[pd$group], 
     cex = 2)
grid()
legend("topleft",
       legend = names(pal_group) %>% unique(),
       pt.bg = pal_group %>% unique(),
       col = "black",
       pch = 21,
       pt.lwd = 1.3,
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       ncol = 1,
       x.intersp = 0.8,
       y.intersp = 1.2,
       inset=c(1,0),
       xpd=TRUE)
dev.off()

#HILIC standard VERY separated - can't see separation of others
#do again without HILIC standard

#get intensity values
ft_ints.fl <- featureValues(data_peaks_grouped_filled)
ft_ints.fl <- as.data.frame(ft_ints.fl)
names(ft_ints.fl) <- pd$name
#log transform
ft_ints.fl <- log2(ft_ints.fl)
#replace all NA values with zero
ft_ints.fl[is.na(ft_ints.fl)] <- 0
#remove HILIC standard
ft_ints.fl <- subset(ft_ints.fl, select=-c(HILICstd))
#perform PCA with the intensities mean centered
pc2 <- prcomp(t(ft_ints.fl), center = TRUE)
#plot
tiff("./analysis/processing_plots/pca/pca-unfilled-featureIntensities.tiff",
     res = 300,
     height = 9,
     width = 12,
     units = "in")
pcSummary2 <- summary(pc2)
par(family = "Avenir",
    mar=c(5,4,2,9))
plot(pc2$x[, 1], 
     pc2$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary2$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary2$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "black", 
     lwd = 1.3,
     bg = pal_group[pd$group], 
     cex = 2)
grid()
legend("topleft",
       legend = names(pal_group) %>% unique(),
       pt.bg = pal_group %>% unique(),
       col = "black",
       pch = 21,
       pt.lwd = 1.3,
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       ncol = 1,
       x.intersp = 0.8,
       y.intersp = 1.2,
       inset=c(1,0),
       xpd=TRUE)
dev.off()

#looks better - continue with filled data
res <- data_peaks_grouped_filled
rm(data, data_peaks_grouped, data_peaks, data_peaks_grouped_filled)


#10: Save diffreport of xdata -----
xset <- as(res, "xcmsSet")
sampnames(xset) <- pData(res)$name
sampclass(xset) <- pData(res)$group

#11. CAMERA workflow: isotope picking / adduct detection----
an <- xsAnnotate(xset)
#group: retention time
an <- groupFWHM(an, 
                 perfwhm = 0.6)
#annotate isotopes
an <- findIsotopes(an, 
                    mzabs=0.01)
#group: correlation info
an <- groupCorr(an, 
                  cor_eic_th=0.75)
#find adducts
an <- findAdducts(an, 
                    polarity="negative")

#12. Peak list filtering and formatting----
#get peak list
pl <-getPeaklist(an)

#filter by blank exclusion (detected peaks)
pl_be <-pl[pl$no.digest==0 & 
               pl$solvent.blank==0,]

#make rownames from rt and mz of features
rownames(pl_be)<-paste(round(pl_be$rt,1),
                       round(pl_be$mz,3),
                       sep="_")
#change NA to 0
pl_be[is.na(pl_be)] <- 0

#filter for isotopes or adducts
pl_be_isoadd <- pl_be[pl_be$isotopes!=""|pl_be$adduct!="",]



#13: PCA plot: peaK list filtering  ----
#do all without HILIC standard
group <- pd$group[pd$group!= "HILIC standard"]
pal_group <- pal_group[grep("HILIC",names(pal_group), invert = T)]

sampleColNames <- pd$name %>% 
    sub("^1", "X1", .) %>% 
    sub("^2", "X2", .) %>% 
    sub("^3", "X3", .) %>% 
    sub("^4", "X4", .) %>% 
    sub("^5", "X5", .) %>% 
    gsub("-", ".", .) %>% 
    sub("sulfate$", "sulfate.1", .)
sampleColNames <- sampleColNames[sampleColNames != "HILICstd"]


##NO FILTERING
#get intensity values
pl_ints1 <- pl[,sampleColNames]
#log transform
pl_ints1 <- log2(pl_ints1)
#replace all NA values with zero
pl_ints1[is.na(pl_ints1)] <- 0
#perform PCA with the intensities mean centered
pc1 <- prcomp(t(pl_ints1), center = TRUE)
#plot
tiff("./analysis/processing_plots/pca/pca-unfilled-peaklist-unfiltered.tiff",
     res = 300,
     height = 9,
     width = 12,
     units = "in")
pcSummary1 <- summary(pc1)
par(family = "Avenir",
    mar=c(5,4,2,9))
plot(pc1$x[, 1], 
     pc1$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary1$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary1$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "black", 
     lwd = 1.3,
     bg = pal_group[group], 
     cex = 2)
grid()
legend("topleft",
       legend = names(pal_group) %>% unique(),
       pt.bg = pal_group %>% unique(),
       col = "black",
       pch = 21,
       pt.lwd = 1.3,
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       ncol = 1,
       x.intersp = 0.8,
       y.intersp = 1.2,
       inset=c(1,0),
       xpd=TRUE)
dev.off()

#exact same as filled features

##BLANK EXCLUSION
#get intensity values
pl_ints2 <- pl_be[,sampleColNames]
#set 0 to NA
pl_ints2[pl_ints2 == 0] <- NA
#log transform
pl_ints2 <- log2(pl_ints2)
#replace all NA values with zero
pl_ints2[is.na(pl_ints2)] <- 0
#perform PCA with the intensities mean centered
pc2 <- prcomp(t(pl_ints2), center = TRUE)
#plot
tiff("./analysis/processing_plots/pca/pca-unfilled-peaklist-blankexclusion.tiff",
     res = 300,
     height = 9,
     width = 12,
     units = "in")
pcSummary2 <- summary(pc2)
par(family = "Avenir",
    mar=c(5,4,2,9))
plot(pc2$x[, 1], 
     pc2$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary2$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary2$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "black", 
     lwd = 1.3,
     bg = pal_group[pd$group], 
     cex = 2)
grid()
legend("topleft",
       legend = names(pal_group) %>% unique(),
       pt.bg = pal_group %>% unique(),
       col = "black",
       pch = 21,
       pt.lwd = 1.3,
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       ncol = 1,
       x.intersp = 0.8,
       y.intersp = 1.2,
       inset=c(1,0),
       xpd=TRUE)
dev.off()

##BLANK EXCLUSION AND ISOTOPE OR ADDUCT FILTERED
#get intensity values
pl_ints3 <- pl_be_isoadd[,sampleColNames]
#set 0 to NA
pl_ints3[pl_ints3 == 0] <- NA
#log transform
pl_ints3 <- log2(pl_ints3)
#replace all NA values with zero
pl_ints3[is.na(pl_ints3)] <- 0
#perform PCA with the intensities mean centered
pc3 <- prcomp(t(pl_ints3), center = TRUE)
#plot
tiff("./analysis/processing_plots/pca/pca-unfilled-peaklist-blankexclusion-iso_add.tiff",
     res = 300,
     height = 9,
     width = 12,
     units = "in")
pcSummary3 <- summary(pc3)
par(family = "Avenir",
    mar=c(5,4,2,9))
plot(pc3$x[, 1], 
     pc3$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary3$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary3$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "black", 
     lwd = 1.3,
     bg = pal_group[pd$group], 
     cex = 2)
grid()
legend("topleft",
       legend = names(pal_group) %>% unique(),
       pt.bg = pal_group %>% unique(),
       col = "black",
       pch = 21,
       pt.lwd = 1.3,
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       ncol = 1,
       x.intersp = 0.8,
       y.intersp = 1.2,
       inset=c(1,0),
       xpd=TRUE)
dev.off()

#WITH MORE FILTERING DIGESTS (APART FROM GH92 AND GH93+GH99-SULPHATASE)
#GET CLOSER TO SOLVENT BLANKS...
#not ideal but maybe reflects reality????

#14: Collapse features with multiple isotopes -----
#for now work with blank exclusion list
peaks <- pl_be
setDT(peaks)
#split out features without an isotope detected
peaks_noiso <- peaks[peaks$isotopes=="",]
peaks_iso <- peaks[!peaks$isotopes=="",]
#make column for the isotope group
peaks_iso$isotope_group <- peaks_iso$isotopes %>% 
    sub("\\[M.*", "", .)
#order isotopes within each group correctly
peaks_iso$isotope_number <- peaks_iso$isotopes %>% 
    sub(".*\\[M\\].*", "0", .) %>% 
    sub(".*\\[M\\+", "", .) %>% 
    sub("\\].*", "", .) %>% 
    as.numeric()
peaks_iso <- peaks_iso[order(isotope_group, 
                             isotope_number),]
#get concatenated list of isotopes per group
iso_concat <- peaks_iso[, 
                        list(isotopes = paste(isotopes, 
                                              collapse = ', ')), 
                        by = isotope_group]
#remove duplicates within each isotope group (will keep [M] isotope)
#because of ordering
peaks_iso <- unique(peaks_iso, 
                    by = "isotope_group")
#merge to get concatenated isotope lists
peaks_iso <- merge(peaks_iso,
                   iso_concat,
                   by = "isotope_group")
#clean up df
peaks_iso <- peaks_iso %>% 
    select(-c("isotope_group",
              "isotope_number",
              "isotopes.x"))
names(peaks_iso)[names(peaks_iso) == 'isotopes.y'] <- 'isotopes'

#merge features with and without isotopes
peaks <- rbind.fill(peaks_noiso,
                    peaks_iso)

#15: Annotate features based on predictions----
#import prediction table
#see https://github.com/margotbligh/sugarMassesPredict
#created with: ssugarMassesPredict.py -dp 1 8 -p 0 -m unsaturated carboxyl deoxy nacetyl sulphate -n 2 -ds 1 -i neg -s 175 1400
mz_predicted <- fread("predicted_sugars.txt")

#remove "extra" columns
extraCol <- c("level_0",
              "index",
              'mass',
              'formula')

mz_predicted <- mz_predicted %>% 
    select(-all_of(extraCol))

#make long format
predicted <- mz_predicted %>% 
    gather(key = "ion",
           value = "mz",
           -name,
           -dp)

#make data.table
setDT(predicted)
setDT(peaks)
#create interval to overlap with (same width as for peak grouping)
predicted$mz <- as.numeric(predicted$mz)
predicted$mzmin <- predicted$mz-0.001
predicted$mzmax <- predicted$mz+0.001

#remove NA rows
predicted <- na.omit(predicted)

#match using foverlaps from data.table (very fast)
setkey(predicted, mzmin, mzmax)
peaks_nopred <- peaks
peaks <- foverlaps(peaks,
                   predicted)

#change NA values created during matching (features with no match) to be blank
#remove extra columns
peaks$mzmin <- NULL
peaks$mzmax <- NULL

peaks <- peaks %>% 
    replace_na(list("name"="",
                    "ion"= "", 
                    "mz" = "",
                    "dp" = ""))


#only keep matched features
peaks_matched <- peaks[!peaks$name=="",]
peaks_unmatched <- peaks[peaks$name=="",]

#make annotation column
peaks_matched$theoretical <- paste0(peaks_matched$name, 
                                    ":",
                                    peaks_matched$ion)

#aggregate so that if there are multiple predictions for one feature
#they are shown in the same row

#set column names to aggregate on
colNames <- setdiff(names(peaks_matched), 
                    names(predicted))
colNames <- colNames[!colNames == "theoretical"]
#aggregate
setDT(peaks_matched)
peaks_matched <- peaks_matched[, 
                               list(annotation = paste(theoretical, 
                                                       collapse = ', ')), 
                               by = colNames]

#rename columns in matched
names(peaks_matched)[names(peaks_matched) == 'i.mz'] <- 'mz'
names(peaks_matched)[names(peaks_matched) == 'i.mzmin'] <- 'mzmin'
names(peaks_matched)[names(peaks_matched) == 'i.mzmax'] <- 'mzmax'

#add / delete /rename coluumns of unmatched to be the same
peaks_unmatched <- peaks_unmatched %>% dplyr::select(all_of(colNames))
names(peaks_unmatched)[names(peaks_unmatched) == 'i.mz'] <- 'mz'
names(peaks_unmatched)[names(peaks_unmatched) == 'i.mzmin'] <- 'mzmin'
names(peaks_unmatched)[names(peaks_unmatched) == 'i.mzmax'] <- 'mzmax'
peaks_unmatched$annotation <- paste0("unknown_mz",
                                     round(peaks_unmatched$mz, 3))

#combine matched and unmatched
peaks_all <- rbind(peaks_matched,
                   peaks_unmatched)

#write to table
fwrite(peaks_matched,
       "./analysis/analysis_tables/matched-peaks-v1.txt",
       sep = "\t")

fwrite(peaks_all,
       "./analysis/analysis_tables/all-peaks-v1.txt",
       sep = "\t")

#16: Plot EIC for sulphated mannose ----
chr_259 <- chromatogram(res,
                        mz = c(259.0129306 - 0.001,
                               259.0129306 + 0.001))
chr_259.df <- data.frame(sample = as.character(),
                         group = as.character(),
                         rt = as.numeric(),
                         intensity = as.numeric())

for (i in 1:length(chr_259$name)){
    rt = chr_259[[i]]@rtime/60
    intensity = chr_259[[i]]@intensity
    sample = rep(chr_259$name[i], length(rt))
    group = rep(chr_259$group[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity)
    chr_259.df <- rbind(chr_259.df,
                        temp)
}

chr_259.df <- chr_259.df %>% filter(group != "HILIC standard")

#set NA to 0
chr_259.df[is.na(chr_259.df)] <- 0

#change groups
groups2 <- chr_259.df$group %>% 
    sub("PLx.*|^GH99.*", "PLx and GH99-sulphatase digests", .) %>% 
    sub("^no.*|solvent.*", "negative controls", .) %>% 
    unique()

chr_259.df$group <- chr_259.df$group %>% 
    sub("PLx.*|^GH99.*", "PLx and GH99-sulphatase digests", .) %>% 
    sub("^no.*|solvent.*", "negative controls", .)
chr_259.df$group <- factor(chr_259.df$group,
                           levels = unique(chr_259.df$group))
pal_group2 <- hcl.colors(n = length(groups2),
                         palette = "Dark3")
names(pal_group2) <- groups2

#make directory
dir.create("./analysis/analysis_plots/mannose-sulphate-eic",
           showWarnings = FALSE)

#make labelling function
scientific_function <- function(x) {
    text <- gsub("E0", "", gsub("e\\+0", "E", scales::scientific_format()(x)))
    text
}

#format
chr_259.df$group2 <- chr_259.df$group %>% 
    sub("mannose.*", "standards", .) %>% 
    sub(".*digest.*", "digests", .)

chr_259.df$group2 <- factor(chr_259.df$group2, 
                            levels = unique(chr_259.df$group2))

chr_259.df$group <- chr_259.df$group %>% 
    sub("sulfate", "sulphate", .)

chr_259.df$group <- factor(chr_259.df$group,
                           levels = unique(chr_259.df$group))

chr_259.df$sample <- factor(chr_259.df$sample,
                           levels = rev(unique(chr_259.df$sample)))

#plot eic
tiff("./analysis/analysis_plots/mannose-sulphate-eic/facet-twogroups-v2.tiff",
     width = 9,
     height = 6,
     res = 600,
     units = "in")
ggplot() +
    geom_line(mapping = aes(rt,
                            intensity,
                            colour = group,
                            group = sample),
              data = chr_259.df[chr_259.df$group!="negative controls" &
                                    chr_259.df$group!="PLx and GH99-sulphatase digests",],
              lwd = 1) +
    scale_colour_manual(values = c(rep("#F9A100", 2),
                                   rep("#46494B", 4)),
                        name = "") +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group),
               scales = "free_y") +
    xlim(0, 8) +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, 
                                      family = "Avenir LT 65 Medium",
                                      angle = 360,
                                      hjust = 0),
          strip.background = element_blank(),
          axis.text = element_text(size = 12,
                                   family = "Avenir"),
          axis.title = element_text(size = 12, 
                                    family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          legend.position = "none",
          legend.text = element_text(size = 12, 
                                     family = "Avenir"),
          axis.line = element_blank())
dev.off()

#17: Plot MS2 for sulphated mannose ----
#subset ms2 data to remove samples with only ms1
data_ms2 <- filterFile(data_ms2,
                       file = grep("HILIC", data_ms2$name, invert = T))

#assign peak info (didn't pick peaks on ms2 data)
#HILIC sample is number 14 - need to remove peaks from this sample
#and reorder sample numbering from that
cp <- chromPeaks(res)
cp <- as.data.frame(cp)
cp <- cp %>% filter(sample != 14)
cp$sample.new <- cp$sample
for (i in 15:19){
    cp$sample.new[cp$sample == i] <- i-1
}
cp$sample <- cp$sample.new
cp$sample.new <- NULL
cp <- cp %>% 
    filter(between(mz, 259.0129306-0.0005, 259.0129306+0.0005))
cp <- as.matrix(cp)
chromPeaks(data_ms2) <- cp

##extract ms2 associated with chromPeaks
ms2_259 <- chromPeakSpectra(data_ms2,
                            msLevel = 2,
                            expandRt = 2.5,
                            expandMz = 0.005,
                            skipFilled = FALSE,
                            method = "all",
                            return.type = "Spectra")


#remove zero intensity masses
ms2_259@listData <- lapply(ms2_259, 
                           clean, 
                           all = TRUE)

#combine spectra by peak within each file
ms2_259_comb <- combineSpectra(ms2_259,
                               fcol = 'peak_id',
                               mzd = 0.005,
                               intensityFun = mean)

#see how many spectra per sample
table(fromFile(ms2_259_comb))

#normalise with respect to ion with highest intensity
ms2_259_comb_old <- ms2_259_comb
ms2_259_comb@listData <- lapply(ms2_259_comb_old, 
                                normalise, 
                                method = "max")


#plot: mannose-4-sulphate
mz <- ms2_259_comb[["CP25849.F18.S00567"]]@mz
intensity <- ms2_259_comb[["CP25849.F18.S00567"]]@intensity*100
labels <- as.character(round(mz, 3))
labels[intensity < 5] <- ""
rt <- round(ms2_259_comb[["CP25849.F18.S00567"]]@rt/60, 1)


tiff(filename = paste0("./analysis/ms2_plots/mannose-4-sulphate_",
                       rt,
                       "min.tiff"),
     height = 6, width = 6, units = "in", res = 300)
cairo_pdf(filename = paste0("./analysis/ms2_plots/mannose-4-sulphate_",
                            rt,
                            "min.pdf"),
          height = 6, width = 6)
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
              nudge_y = 10,
              size = 5.5,
              family = "Avenir") +
    scale_x_continuous(breaks = seq(0,300, by = 50),
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

#plot: mannose-3-sulphate
mz <- ms2_259_comb[["CP89336.F17.S01081"]]@mz
intensity <- ms2_259_comb[["CP89336.F17.S01081"]]@intensity*100
labels <- as.character(round(mz, 3))
labels[intensity < 5] <- ""
rt <- round(ms2_259_comb[["CP89336.F17.S01081"]]@rt/60, 1)


tiff(filename = paste0("./analysis/ms2_plots/mannose-3-sulphate_",
                       rt,
                       "min.tiff"),
     height = 6, width = 6, units = "in", res = 300)
cairo_pdf(filename = paste0("./analysis/ms2_plots/mannose-3-sulphate_",
                            rt,
                            "min.pdf"),
          height = 6, width = 6)

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
              nudge_y = 10,
              size = 5.5,
              family = "Avenir") +
    scale_x_continuous(breaks = seq(0,300, by = 50),
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

#plot: mannose-2-sulphate
mz <- ms2_259_comb[["CP20439.F15.S01236"]]@mz
intensity <- ms2_259_comb[["CP20439.F15.S01236"]]@intensity*100
labels <- as.character(round(mz, 3))
labels[intensity < 3] <- ""
rt <- round(ms2_259_comb[["CP20439.F15.S01236"]]@rt/60, 1)


tiff(filename = paste0("./analysis/ms2_plots/mannose-2-sulphate_",
                       rt,
                       "min.tiff"),
     height = 6, width = 6, units = "in", res = 300)
cairo_pdf(filename = paste0("./analysis/ms2_plots/mannose-2-sulphate_",
                            rt,
                            "min.pdf"),
          height = 6, width = 6)
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
              nudge_y = 10,
              size = 5.5,
              family = "Avenir") +
    scale_x_continuous(breaks = seq(0,300, by = 50),
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

#plot: mannose-6-sulphate
mz <- ms2_259_comb[["CP27931.F19.S01477"]]@mz
intensity <- ms2_259_comb[["CP27931.F19.S01477"]]@intensity*100
labels <- as.character(round(mz, 3))
labels[intensity < 5] <- ""
rt <- round(ms2_259_comb[["CP27931.F19.S01477"]]@rt/60, 1)


tiff(filename = paste0("./analysis/ms2_plots/mannose-6-sulphate_",
                       rt,
                       "min.tiff"),
     height = 6, width = 6, units = "in", res = 300)

cairo_pdf(filename = paste0("./analysis/ms2_plots/mannose-4-sulphate_",
                            rt,
                            "min.pdf"),
          height = 6, width = 6)

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
              nudge_y = 10,
              size = 5.5,
              family = "Avenir") +
    scale_x_continuous(breaks = seq(0,300, by = 50),
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

#plot: digest 5a
mz <- ms2_259_comb[["CP11441.F09.S01572"]]@mz
intensity <- ms2_259_comb[["CP11441.F09.S01572"]]@intensity*100
labels <- as.character(round(mz, 3))
labels[intensity < 7.5] <- ""
rt <- round(ms2_259_comb[["CP11441.F09.S01572"]]@rt/60, 1)


tiff(filename = paste0("./analysis/ms2_plots/digest-5a_sulphatedmannose_",
                       rt,
                       "min.tiff"),
     height = 6, width = 6, units = "in", res = 300)
cairo_pdf(filename = paste0("./analysis/ms2_plots/digest-5a_sulphatedmannose_",
                            rt,
                            "min.pdf"),
          height = 6, width = 6)
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
              nudge_y = 10,
              size = 5.5,
              family = "Avenir") +
    scale_x_continuous(breaks = seq(0,300, by = 50),
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

#plot: digest 5b
mz <- ms2_259_comb[["CP12834.F10.S01451"]]@mz
intensity <- ms2_259_comb[["CP12834.F10.S01451"]]@intensity*100
labels <- as.character(round(mz, 3))
labels[intensity < 8] <- ""
rt <- round(ms2_259_comb[["CP12834.F10.S01451"]]@rt/60, 1)


tiff(filename = paste0("./analysis/ms2_plots/digest-5b_sulphatedmannose_",
                       rt,
                       "min.tiff"),
     height = 6, width = 6, units = "in", res = 300)
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
              nudge_y = 10,
              size = 5.5,
              family = "Avenir") +
    scale_x_continuous(breaks = seq(0,300, by = 50),
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

#17: Screening for GH99 digest products -----
gh99.df1 <- peaks_all %>% filter(GH99.sulphatase.digest >= 1 & 
                         no.digest == 0 & 
                         mannose.2.sulfate == 0 & 
                         mannose.3.sulfate == 0 & 
                         mannose.4.sulfate == 0 & 
                         mannose.6.sulfate == 0 & 
                         GH92.digest == 0 & 
                         GH92...GH99.sulphatase.digest == 0)

gh99.df1 <- cbind(rt_minutes = round(gh99.df1$rt/60, 1),
                  gh99.df1)

gh99.mz <- round(gh99.df1$mz, 3) %>% unique()
gh99.df1$mz_round <- round(gh99.df1$mz, 3)
gh99.df2 <- gh99.df1[gh99.df1$mz_round %in% gh99.mz,] %>% 
    distinct(mz_round, .keep_all = T) %>% 
    dplyr::select(c(annotation, mz_round, rt))
gh99.ion <- gh99.df2$annotation
gh99.rt <- gh99.df2$rt/60

gh99_data <- filterFile(res,
                        file = which(res$group == "GH99-sulphatase digest"))

error = 0.0025
gh99_chrlist <- list()
for (i in 1:length(gh99.mz)){
    mzr = c(gh99.mz[i] - error,
            gh99.mz[i] + error)
    gh99_chrlist[[i]] <- chromatogram(gh99_data, 
                                      mz = mzr)
}

#save list
save(gh99_chrlist,
     file = "./analysis/RData/gh99_chrlist.RData")

#extract intensity and rt values
gh99_int_list <- list()
for (i in 1:length(gh99_data$name)){
    gh99_int_list[[i]] <- lapply(gh99_chrlist, function(x) {
        x[[i]]@intensity
    }) 
}

gh99_rt_list <- list()
for (i in 1:length(gh99_data$name)){
    gh99_rt_list[[i]] <- lapply(gh99_chrlist, function(x) {
        x[[i]]@rtime
    }) 
}

#build data frame (long format)
gh99_chr.df <- data.frame(sample = as.character(),
                          ion = as.character(),
                          mz = as.numeric(),
                          rt = as.numeric(),
                          intensity = as.numeric())

for (i in 1:length(gh99_data$name)){
    for (j in 1:length(gh99.mz)){
        rt = gh99_rt_list[[i]][[j]]/60
        intensity = gh99_int_list[[i]][[j]]
        sample = rep(gh99_data$name[i], length(rt))
        ion = rep(gh99.ion[j], length(rt))
        mz = rep(gh99.mz[j], length(rt))
        temp <- data.frame(sample = sample,
                           ion = ion,
                           mz = mz,
                           rt = rt,
                           intensity = intensity)
        gh99_chr.df <- rbind(gh99_chr.df,
                         temp)
    }
}



#set NA to 0
gh99_chr.df[is.na(gh99_chr.df)] <- 0

dir.create("./analysis/analysis_plots/eic_checking",
           showWarnings = FALSE)
dir.create("./analysis/analysis_plots/eic_checking/gh99",
           showWarnings = FALSE)

#make palette
gh99_pal <- c("#F9A100",
              "#46494B")

#plot for each ion the full and restricted rt chromatograms by group
for (i in 1:length(gh99.mz)){
    ion = gh99.ion[i]
    mz = gh99.mz[i]
    rt = gh99.rt[i]
    rtmin = rt - 2.5
    if (rtmin < 0) {rtmin <- 0}
    rtmax = rt + 2.5
    
    df1 <- gh99_chr.df %>% 
        filter(between(rt, rtmin, rtmax)) %>% 
        filter(ion == !!ion)
    
    df2 <- gh99_chr.df %>% 
        filter(ion == !!ion)
    
    #plot restricted retention time
    p1 <- ggplot() +
        geom_line(mapping = aes(rt,
                                intensity,
                                colour = sample,
                                group = sample),
                  data = df1,
                  lwd = 1.2) +
        scale_colour_manual(values = gh99_pal,
                            name = "") +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(sample),
                   scales = "free_y") +
        theme_classic() +
        theme(strip.text = element_blank(),
              axis.text = element_text(size = 12,
                                       family = "Avenir"),
              axis.title = element_text(size = 14, 
                                        family = "Avenir LT 65 Medium"),
              panel.border = element_rect(colour = "#848587",
                                          size = 0.5,
                                          fill = NA),
              legend.text = element_text(size = 12, 
                                         family = "Avenir"),
              axis.line = element_blank())
    
    #plot full retention time
    p2 <- ggplot() +
        geom_line(mapping = aes(rt,
                                intensity,
                                colour = sample,
                                group = sample),
                  data = df2,
                  lwd = 1.2) +
        scale_colour_manual(values = gh99_pal,
                            name = "") +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(sample)) +
        theme_classic() +
        theme(strip.text = element_blank(),
              axis.text = element_text(size = 12,
                                       family = "Avenir"),
              axis.title = element_text(size = 14, 
                                        family = "Avenir LT 65 Medium"),
              panel.border = element_rect(colour = "#848587",
                                          size = 0.5,
                                          fill = NA),
              legend.text = element_text(size = 12, 
                                         family = "Avenir"),
              axis.line = element_blank())
    
    #combine plots
    combined <- p1 + p2 & 
        theme(legend.position = "bottom",
              plot.title = element_text(size = 12, 
                                        family = "Avenir LT 65 Medium",
                                        hjust = 0.5))
    
    #save as pdf
    cairo_pdf(paste0("./analysis/analysis_plots/eic_checking/gh99/",
                     ion, "_mz", mz,
                     ".pdf"),
              width = 12,
              height = 9)
    print(combined +
              plot_layout(ncol=2, guides = "collect") + 
              plot_annotation(title = paste0("ion = ",
                                             ion,
                                             "; m/z = ",
                                             mz)))
    
    dev.off()
}




#18: Screening for PLx digest products -----
PLx.df1 <- peaks_all %>% filter(PLx.digest >= 1 & 
                                     no.digest == 0 & 
                                     mannose.2.sulfate == 0 & 
                                     mannose.3.sulfate == 0 & 
                                     mannose.4.sulfate == 0 & 
                                     mannose.6.sulfate == 0 & 
                                     GH92.digest == 0 & 
                                     GH92...GH99.sulphatase.digest == 0)

PLx.df1 <- cbind(rt_minutes = round(PLx.df1$rt/60, 1),
                 PLx.df1)

PLx.mz <- round(PLx.df1$mz, 3) %>% unique()
PLx.df1$mz_round <- round(PLx.df1$mz, 3)
PLx.df2 <- PLx.df1[PLx.df1$mz_round %in% PLx.mz,] %>% 
    distinct(mz_round, .keep_all = T) %>% 
    dplyr::select(c(annotation, mz_round, rt))
PLx.ion <- PLx.df2$annotation
PLx.rt <- PLx.df2$rt/60

PLx_data <- filterFile(res,
                        file = which(res$group == "PLx digest"))

error = 0.0025
PLx_chrlist <- list()
for (i in 1:length(PLx.mz)){
    mzr = c(PLx.mz[i] - error,
            PLx.mz[i] + error)
    PLx_chrlist[[i]] <- chromatogram(PLx_data, 
                                      mz = mzr)
}

#save list
save(PLx_chrlist,
     file = "./analysis/RData/PLx_chrlist.RData")

#extract intensity and rt values
PLx_int_list <- list()
for (i in 1:length(PLx_data$name)){
    PLx_int_list[[i]] <- lapply(PLx_chrlist, function(x) {
        x[[i]]@intensity
    }) 
}

PLx_rt_list <- list()
for (i in 1:length(PLx_data$name)){
    PLx_rt_list[[i]] <- lapply(PLx_chrlist, function(x) {
        x[[i]]@rtime
    }) 
}

#build data frame (long format)
PLx_chr.df <- data.frame(sample = as.character(),
                         ion = as.character(),
                         mz = as.numeric(),
                         rt = as.numeric(),
                         intensity = as.numeric())

for (i in 1:length(PLx_data$name)){
    for (j in 1:length(PLx.mz)){
        rt = PLx_rt_list[[i]][[j]]/60
        intensity = PLx_int_list[[i]][[j]]
        sample = rep(PLx_data$name[i], length(rt))
        ion = rep(PLx.ion[j], length(rt))
        mz = rep(PLx.mz[j], length(rt))
        temp <- data.frame(sample = sample,
                           ion = ion,
                           mz = mz,
                           rt = rt,
                           intensity = intensity)
        PLx_chr.df <- rbind(PLx_chr.df,
                            temp)
    }
}



#set NA to 0
PLx_chr.df[is.na(PLx_chr.df)] <- 0

dir.create("./analysis/analysis_plots/eic_checking/PLx",
           showWarnings = FALSE)

#make palette
PLx_pal <- c("#F9A100",
              "#46494B")

#plot for each ion the full and restricted rt chromatograms by group
for (i in 1:length(PLx.mz)){
    ion = PLx.ion[i]
    mz = PLx.mz[i]
    rt = PLx.rt[i]
    rtmin = rt - 2.5
    if (rtmin < 0) {rtmin <- 0}
    rtmax = rt + 2.5
    
    df1 <- PLx_chr.df %>% 
        filter(between(rt, rtmin, rtmax)) %>% 
        filter(ion == !!ion)
    
    df2 <- PLx_chr.df %>% 
        filter(ion == !!ion)
    
    #plot restricted retention time
    p1 <- ggplot() +
        geom_line(mapping = aes(rt,
                                intensity,
                                colour = sample,
                                group = sample),
                  data = df1,
                  lwd = 1.2) +
        scale_colour_manual(values = PLx_pal,
                            name = "") +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(sample),
                   scales = "free_y") +
        theme_classic() +
        theme(strip.text = element_blank(),
              axis.text = element_text(size = 12,
                                       family = "Avenir"),
              axis.title = element_text(size = 14, 
                                        family = "Avenir LT 65 Medium"),
              panel.border = element_rect(colour = "#848587",
                                          size = 0.5,
                                          fill = NA),
              legend.text = element_text(size = 12, 
                                         family = "Avenir"),
              axis.line = element_blank())
    
    #plot full retention time
    p2 <- ggplot() +
        geom_line(mapping = aes(rt,
                                intensity,
                                colour = sample,
                                group = sample),
                  data = df2,
                  lwd = 1.2) +
        scale_colour_manual(values = PLx_pal,
                            name = "") +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(sample)) +
        theme_classic() +
        theme(strip.text = element_blank(),
              axis.text = element_text(size = 12,
                                       family = "Avenir"),
              axis.title = element_text(size = 14, 
                                        family = "Avenir LT 65 Medium"),
              panel.border = element_rect(colour = "#848587",
                                          size = 0.5,
                                          fill = NA),
              legend.text = element_text(size = 12, 
                                         family = "Avenir"),
              axis.line = element_blank())
    
    #combine plots
    combined <- p1 + p2 & 
        theme(legend.position = "bottom",
              plot.title = element_text(size = 12, 
                                        family = "Avenir LT 65 Medium",
                                        hjust = 0.5))
    
    #save as pdf
    cairo_pdf(paste0("./analysis/analysis_plots/eic_checking/PLx/",
                     ion, "_mz", mz,
                     ".pdf"),
              width = 12,
              height = 9)
    print(combined +
              plot_layout(ncol=2, guides = "collect") + 
              plot_annotation(title = paste0("ion = ",
                                             ion,
                                             "; m/z = ",
                                             mz)))
    
    dev.off()
}




#19: Re-pick peaks with less stringent settings----
#set parameters
cwp<-CentWaveParam()
cwp@ppm<-5
cwp@peakwidth<-c(10,55)
cwp@snthresh<-5
cwp@noise <- 5000
cwp@prefilter <- c(5, 1000)
cwp@mzdiff <- 0.001

#pick peaks
data_peaks2 <- findChromPeaks(res, param = cwp)

save(data_peaks2,
     file = "./analysis/RData/data_peaks2.RData")

#20: Re-group peaks to create "features"---------
#parameters
pdp <- PeakDensityParam(sampleGroups = data_peaks2$group,
                        binSize = 0.005,
                        bw = 6,
                        minFraction = 0.5) 
#group peaks
data_peaks_grouped2 <- groupChromPeaks(data_peaks2, param = pdp)

save(data_peaks_grouped2, 
     file = "./analysis/RData/data_peaks_grouped2.RData")

#21: Fill in missing peaks----------
fpp <- FillChromPeaksParam()
data_peaks_grouped_filled2 <- fillChromPeaks(data_peaks_grouped2)

save(data_peaks_grouped_filled2, 
     file = "./analysis/RData/data_peaks_grouped_filled2.RData")

#22: Save diffreport of xdata -----
xset2 <- as(data_peaks_grouped_filled2, "xcmsSet")
sampnames(xset2) <- pData(data_peaks_grouped_filled2)$name
sampclass(xset2) <- pData(data_peaks_grouped_filled2)$group

#23. CAMERA workflow: isotope picking / adduct detection----
an2 <- xsAnnotate(xset2)
#group: retention time
an2 <- groupFWHM(an2, 
                 perfwhm = 0.6)
#annotate isotopes
an2 <- findIsotopes(an2, 
                    mzabs=0.01)
#group: correlation info
an2 <- groupCorr(an2, 
                 cor_eic_th=0.75)
#find adducts
an2 <- findAdducts(an2, 
                   polarity="negative")


#24. Peak list filtering and formatting----
#get peak list
pl <-getPeaklist(an2)

#filter by blank exclusion (detected peaks)
pl_be <-pl[pl$no.digest==0 & 
               pl$solvent.blank==0,]

#make rownames from rt and mz of features
rownames(pl_be)<-paste(round(pl_be$rt,1),
                       round(pl_be$mz,3),
                       sep="_")
#change NA to 0
pl_be[is.na(pl_be)] <- 0

#filter for isotopes or adducts
pl_be_isoadd <- pl_be[pl_be$isotopes!=""|pl_be$adduct!="",]

#25: PCA plot: peaK list filtering  ----
#do all without HILIC standard
sampleColNames <- pd$name %>% 
    sub("^1", "X1", .) %>% 
    sub("^2", "X2", .) %>% 
    sub("^3", "X3", .) %>% 
    sub("^4", "X4", .) %>% 
    sub("^5", "X5", .) %>% 
    gsub("-", ".", .) %>% 
    sub("sulfate$", "sulfate.1", .)
sampleColNames <- sampleColNames[sampleColNames != "HILICstd"]


##NO FILTERING
#get intensity values
pl_ints1 <- pl[,sampleColNames]
#log transform
pl_ints1 <- log2(pl_ints1)
#replace all NA values with zero
pl_ints1[is.na(pl_ints1)] <- 0
#perform PCA with the intensities mean centered
pc1 <- prcomp(t(pl_ints1), center = TRUE)
#plot
tiff("./analysis/processing_plots/pca/pca-unfilled-peaklist-unfiltered_v2.tiff",
     res = 300,
     height = 9,
     width = 12,
     units = "in")
pcSummary1 <- summary(pc1)
par(family = "Avenir",
    mar=c(5,4,2,9))
plot(pc1$x[, 1], 
     pc1$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary1$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary1$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "black", 
     lwd = 1.3,
     bg = pal_group[group], 
     cex = 2)
grid()
legend("topleft",
       legend = names(pal_group) %>% unique(),
       pt.bg = pal_group %>% unique(),
       col = "black",
       pch = 21,
       pt.lwd = 1.3,
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       ncol = 1,
       x.intersp = 0.8,
       y.intersp = 1.2,
       inset=c(1,0),
       xpd=TRUE)
text(pc1$x[, 1], 
     pc1$x[,2],
     labels = rownames(pc1$x),
     pos = 1,
     col = "black",
     cex = 0.8)
dev.off()

#looks quite ok...
#slightly more separation in digests

##BLANK EXCLUSION
#get intensity values
pl_ints2 <- pl_be[,sampleColNames]
#set 0 to NA
pl_ints2[pl_ints2 == 0] <- NA
#log transform
pl_ints2 <- log2(pl_ints2)
#replace all NA values with zero
pl_ints2[is.na(pl_ints2)] <- 0
#perform PCA with the intensities mean centered
pc2 <- prcomp(t(pl_ints2), center = TRUE)
#plot
tiff("./analysis/processing_plots/pca/pca-unfilled-peaklist-blankexclusion_v2.tiff",
     res = 300,
     height = 9,
     width = 12,
     units = "in")
pcSummary2 <- summary(pc2)
par(family = "Avenir",
    mar=c(5,4,2,9))
plot(pc2$x[, 1], 
     pc2$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary2$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary2$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "black", 
     lwd = 1.3,
     bg = pal_group[group], 
     cex = 2)
grid()
legend("topleft",
       legend = names(pal_group) %>% unique(),
       pt.bg = pal_group %>% unique(),
       col = "black",
       pch = 21,
       pt.lwd = 1.3,
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       ncol = 1,
       x.intersp = 0.8,
       y.intersp = 1.2,
       inset=c(1,0),
       xpd=TRUE)
text(pc2$x[, 1], 
     pc2$x[,2],
     labels = rownames(pc2$x),
     pos = 1,
     col = "black",
     cex = 0.8)
dev.off()

##BLANK EXCLUSION AND ISOTOPE OR ADDUCT FILTERED
#get intensity values
pl_ints3 <- pl_be_isoadd[,sampleColNames]
#set 0 to NA
pl_ints3[pl_ints3 == 0] <- NA
#log transform
pl_ints3 <- log2(pl_ints3)
#replace all NA values with zero
pl_ints3[is.na(pl_ints3)] <- 0
#perform PCA with the intensities mean centered
pc3 <- prcomp(t(pl_ints3), center = TRUE)
#plot
tiff("./analysis/processing_plots/pca/pca-unfilled-peaklist-blankexclusion-iso_add_v2.tiff",
     res = 300,
     height = 9,
     width = 12,
     units = "in")
pcSummary3 <- summary(pc3)
par(family = "Avenir",
    mar=c(5,4,2,9))
plot(pc3$x[, 1], 
     pc3$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary3$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary3$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "black", 
     lwd = 1.3,
     bg = pal_group[group], 
     cex = 2)
grid()
legend("topleft",
       legend = names(pal_group) %>% unique(),
       pt.bg = pal_group %>% unique(),
       col = "black",
       pch = 21,
       pt.lwd = 1.3,
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       ncol = 1,
       x.intersp = 0.8,
       y.intersp = 1.2,
       inset=c(1,0),
       xpd=TRUE)
text(pc3$x[, 1], 
     pc3$x[,2],
     labels = rownames(pc3$x),
     pos = 1,
     col = "black",
     cex = 0.8)
dev.off()

#WITH MORE FILTERING DIGESTS (APART FROM GH92 AND GH93+GH99-SULPHATASE)
#GET CLOSER TO SOLVENT BLANKS...
#not ideal but maybe reflects reality????

#26: Collapse features with multiple isotopes -----
#for now work with blank exclusion list
peaks <- pl_be
setDT(peaks)
#split out features without an isotope detected
peaks_noiso <- peaks[peaks$isotopes=="",]
peaks_iso <- peaks[!peaks$isotopes=="",]
#make column for the isotope group
peaks_iso$isotope_group <- peaks_iso$isotopes %>% 
    sub("\\[M.*", "", .)
#order isotopes within each group correctly
peaks_iso$isotope_number <- peaks_iso$isotopes %>% 
    sub(".*\\[M\\].*", "0", .) %>% 
    sub(".*\\[M\\+", "", .) %>% 
    sub("\\].*", "", .) %>% 
    as.numeric()
peaks_iso <- peaks_iso[order(isotope_group, 
                             isotope_number),]
#get concatenated list of isotopes per group
iso_concat <- peaks_iso[, 
                        list(isotopes = paste(isotopes, 
                                              collapse = ', ')), 
                        by = isotope_group]
#remove duplicates within each isotope group (will keep [M] isotope)
#because of ordering
peaks_iso <- unique(peaks_iso, 
                    by = "isotope_group")
#merge to get concatenated isotope lists
peaks_iso <- merge(peaks_iso,
                   iso_concat,
                   by = "isotope_group")
#clean up df
peaks_iso <- peaks_iso %>% 
    dplyr::select(-c("isotope_group",
                     "isotope_number",
                     "isotopes.x"))
names(peaks_iso)[names(peaks_iso) == 'isotopes.y'] <- 'isotopes'

#merge features with and without isotopes
peaks <- rbind.fill(peaks_noiso,
                    peaks_iso)

#27: Annotate features based on predictions----
#make data.table
setDT(peaks)

#match using foverlaps from data.table (very fast)
setkey(predicted, mzmin, mzmax)
peaks_nopred <- peaks
peaks <- foverlaps(peaks,
                   predicted)

#change NA values created during matching (features with no match) to be blank
#remove extra columns
peaks$mzmin <- NULL
peaks$mzmax <- NULL

peaks <- peaks %>% 
    replace_na(list("name"="",
                    "ion"= "", 
                    "mz" = "",
                    "dp" = ""))


#only keep matched features
peaks_matched <- peaks[!peaks$name=="",]
peaks_unmatched <- peaks[peaks$name=="",]

#make annotation column
peaks_matched$theoretical <- paste0(peaks_matched$name, 
                                    ":",
                                    peaks_matched$ion)

#aggregate so that if there are multiple predictions for one feature
#they are shown in the same row

#set column names to aggregate on
colNames <- setdiff(names(peaks_matched), 
                    names(predicted))
colNames <- colNames[!colNames == "theoretical"]
#aggregate
setDT(peaks_matched)
peaks_matched <- peaks_matched[, 
                               list(annotation = paste(theoretical, 
                                                       collapse = ', ')), 
                               by = colNames]

#rename columns in matched
names(peaks_matched)[names(peaks_matched) == 'i.mz'] <- 'mz'
names(peaks_matched)[names(peaks_matched) == 'i.mzmin'] <- 'mzmin'
names(peaks_matched)[names(peaks_matched) == 'i.mzmax'] <- 'mzmax'

#add / delete /rename coluumns of unmatched to be the same
peaks_unmatched <- peaks_unmatched %>% dplyr::select(all_of(colNames))
names(peaks_unmatched)[names(peaks_unmatched) == 'i.mz'] <- 'mz'
names(peaks_unmatched)[names(peaks_unmatched) == 'i.mzmin'] <- 'mzmin'
names(peaks_unmatched)[names(peaks_unmatched) == 'i.mzmax'] <- 'mzmax'
peaks_unmatched$annotation <- paste0("unknown_mz",
                                     round(peaks_unmatched$mz, 3))

#combine matched and unmatched
peaks_all <- rbind(peaks_matched,
                   peaks_unmatched)

#write to table
fwrite(peaks_matched,
       "./analysis/analysis_tables/matched-peaks-v2.txt",
       sep = "\t")

fwrite(peaks_all,
       "./analysis/analysis_tables/all-peaks-v2.txt",
       sep = "\t")

#17: Screening for GH99 digest products -----
gh99.old.df1 <- gh99.df1
gh99.old.df2 <- gh99.df2
gh99_chr.old.df <- gh99_chr.df

gh99.df1 <- peaks_all %>% filter(GH99.sulphatase.digest >= 1 & 
                                     no.digest == 0 & 
                                     mannose.2.sulfate == 0 & 
                                     mannose.3.sulfate == 0 & 
                                     mannose.4.sulfate == 0 & 
                                     mannose.6.sulfate == 0 & 
                                     GH92.digest == 0 & 
                                     GH92...GH99.sulphatase.digest == 0)

gh99.df1 <- cbind(rt_minutes = round(gh99.df1$rt/60, 1),
                  gh99.df1)

gh99.mz <- round(gh99.df1$mz, 3) %>% unique()
gh99.df1$mz_round <- round(gh99.df1$mz, 3)
gh99.df2 <- gh99.df1[gh99.df1$mz_round %in% gh99.mz,] %>% 
    distinct(mz_round, .keep_all = T) %>% 
    dplyr::select(c(annotation, mz_round, rt))
gh99.ion <- gh99.df2$annotation
gh99.rt <- gh99.df2$rt/60

gh99_data <- filterFile(res,
                        file = which(res$group == "GH99-sulphatase digest"))

error = 0.0025
gh99_chrlist <- list()
for (i in 1:length(gh99.mz)){
    mzr = c(gh99.mz[i] - error,
            gh99.mz[i] + error)
    gh99_chrlist[[i]] <- chromatogram(gh99_data, 
                                      mz = mzr)
}

#save list
save(gh99_chrlist,
     file = "./analysis/RData/gh99_chrlist.RData")

#extract intensity and rt values
gh99_int_list <- list()
for (i in 1:length(gh99_data$name)){
    gh99_int_list[[i]] <- lapply(gh99_chrlist, function(x) {
        x[[i]]@intensity
    }) 
}

gh99_rt_list <- list()
for (i in 1:length(gh99_data$name)){
    gh99_rt_list[[i]] <- lapply(gh99_chrlist, function(x) {
        x[[i]]@rtime
    }) 
}

#build data frame (long format)
gh99_chr.df <- data.frame(sample = as.character(),
                          ion = as.character(),
                          mz = as.numeric(),
                          rt = as.numeric(),
                          intensity = as.numeric())

for (i in 1:length(gh99_data$name)){
    for (j in 1:length(gh99.mz)){
        rt = gh99_rt_list[[i]][[j]]/60
        intensity = gh99_int_list[[i]][[j]]
        sample = rep(gh99_data$name[i], length(rt))
        ion = rep(gh99.ion[j], length(rt))
        mz = rep(gh99.mz[j], length(rt))
        temp <- data.frame(sample = sample,
                           ion = ion,
                           mz = mz,
                           rt = rt,
                           intensity = intensity)
        gh99_chr.df <- rbind(gh99_chr.df,
                             temp)
    }
}



#set NA to 0
gh99_chr.df[is.na(gh99_chr.df)] <- 0

#make directory
dir.create("./analysis/analysis_plots/eic_checking/gh99",
           showWarnings = FALSE)

#make palette
gh99_pal <- c("#F9A100",
              "#46494B")

#plot for each ion the full and restricted rt chromatograms by group
for (i in 1:length(gh99.mz)){
    ion = gh99.ion[i]
    mz = gh99.mz[i]
    rt = gh99.rt[i]
    rtmin = rt - 2.5
    if (rtmin < 0) {rtmin <- 0}
    rtmax = rt + 2.5
    
    df1 <- gh99_chr.df %>% 
        filter(between(rt, rtmin, rtmax)) %>% 
        filter(ion == !!ion)
    
    df2 <- gh99_chr.df %>% 
        filter(ion == !!ion)
    
    #plot restricted retention time
    p1 <- ggplot() +
        geom_line(mapping = aes(rt,
                                intensity,
                                colour = sample,
                                group = sample),
                  data = df1,
                  lwd = 1.2) +
        scale_colour_manual(values = gh99_pal,
                            name = "") +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(sample),
                   scales = "free_y") +
        theme_classic() +
        theme(strip.text = element_blank(),
              axis.text = element_text(size = 12,
                                       family = "Avenir"),
              axis.title = element_text(size = 14, 
                                        family = "Avenir LT 65 Medium"),
              panel.border = element_rect(colour = "#848587",
                                          size = 0.5,
                                          fill = NA),
              legend.text = element_text(size = 12, 
                                         family = "Avenir"),
              axis.line = element_blank())
    
    #plot full retention time
    p2 <- ggplot() +
        geom_line(mapping = aes(rt,
                                intensity,
                                colour = sample,
                                group = sample),
                  data = df2,
                  lwd = 1.2) +
        scale_colour_manual(values = gh99_pal,
                            name = "") +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(sample)) +
        theme_classic() +
        theme(strip.text = element_blank(),
              axis.text = element_text(size = 12,
                                       family = "Avenir"),
              axis.title = element_text(size = 14, 
                                        family = "Avenir LT 65 Medium"),
              panel.border = element_rect(colour = "#848587",
                                          size = 0.5,
                                          fill = NA),
              legend.text = element_text(size = 12, 
                                         family = "Avenir"),
              axis.line = element_blank())
    
    #combine plots
    combined <- p1 + p2 & 
        theme(legend.position = "bottom",
              plot.title = element_text(size = 12, 
                                        family = "Avenir LT 65 Medium",
                                        hjust = 0.5))
    
    #save as pdf
    cairo_pdf(paste0("./analysis/analysis_plots/eic_checking/gh99/",
                     ion, "_mz", mz,
                     ".pdf"),
              width = 12,
              height = 9)
    print(combined +
              plot_layout(ncol=2, guides = "collect") + 
              plot_annotation(title = paste0("ion = ",
                                             ion,
                                             "; m/z = ",
                                             mz)))
    
    dev.off()
}


save(gh99_chrlist,
     file = "./analysis/RData/gh99_chrlist.RData")

#mz to check
gh99.mz.check <- c(177.019, 482.559, 208.960, 258.764, 330.918, 397.218, 
                   406.629, 480.562, 688.834, 708.353, 408.626, 184.915, 
                   216.853, 258.921, 349.899, 380.994, 436.887, 484.871,
                   500.845, 713.456, 228.903)

sgh99_chrlist <- list()
for (i in 1:length(gh99.mz.check)){
    mzr = c(gh99.mz.check[i] - error,
            gh99.mz.check[i] + error)
    gh99_chrlist[[i]] <- chromatogram(data_peaks_grouped_filled2, 
                                      mz = mzr)
}

#extract intensity and rt values
gh99_int_list <- list()
for (i in 1:length(data_peaks_grouped_filled2$name)){
    gh99_int_list[[i]] <- lapply(gh99_chrlist, function(x) {
        x[[i]]@intensity
    }) 
}

gh99_rt_list <- list()
for (i in 1:length(data_peaks_grouped_filled2$name)){
    gh99_rt_list[[i]] <- lapply(gh99_chrlist, function(x) {
        x[[i]]@rtime
    }) 
}

#build data frame (long format)
gh99_chr.df <- data.frame(sample = as.character(),
                          group = as.character(),
                          mz = as.numeric(),
                          rt = as.numeric(),
                          intensity = as.numeric())

for (i in 1:length(data_peaks_grouped_filled2$name)){
    for (j in 1:length(gh99.mz.check)){
        rt = gh99_rt_list[[i]][[j]]/60
        intensity = gh99_int_list[[i]][[j]]
        sample = rep(data_peaks_grouped_filled2$name[i], length(rt))
        group = rep(data_peaks_grouped_filled2$group[i], length(rt))
        mz = rep(gh99.mz.check[j], length(rt))
        temp <- data.frame(sample = sample,
                           group = group,
                           mz = mz,
                           rt = rt,
                           intensity = intensity)
        gh99_chr.df <- rbind(gh99_chr.df,
                             temp)
    }
}



#set NA to 0
gh99_chr.df[is.na(gh99_chr.df)] <- 0

#remove HILIC standard
gh99_chr.df <- gh99_chr.df[gh99_chr.df$group !="HILIC standard",]

#make new groups to plot together
gh99_chr.df$group2 <- gh99_chr.df$group %>% 
    sub("GH92 digest", "digests 4 and 5", .) %>%
    sub(".*\\+.*", "digests 4 and 5", .) %>% 
    sub("^mannose.*", "standards", .)

#make directory
dir.create("./analysis/analysis_plots/eic_checking/gh99",
           showWarnings = FALSE)


#plot for each ion the full and restricted rt chromatograms by group
for (i in 1:length(gh99.mz.check)){
    mz = gh99.mz.check[i]
    rt = gh99.rt[i]
    rtmin = rt - 2.5
    if (rtmin < 0) {rtmin <- 0}
    rtmax = rt + 2.5
    
    df1 <- gh99_chr.df %>% 
        filter(between(rt, rtmin, rtmax)) %>% 
        filter(ion == !!ion)
    
    df2 <- gh99_chr.df %>% 
        filter(ion == !!ion)
    
    #plot restricted retention time
    p1 <- ggplot() +
        geom_line(mapping = aes(rt,
                                intensity,
                                colour = sample,
                                group = sample),
                  data = df1,
                  lwd = 1.2) +
        scale_colour_manual(values = gh99_pal,
                            name = "") +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(sample),
                   scales = "free_y") +
        theme_classic() +
        theme(strip.text = element_blank(),
              axis.text = element_text(size = 12,
                                       family = "Avenir"),
              axis.title = element_text(size = 14, 
                                        family = "Avenir LT 65 Medium"),
              panel.border = element_rect(colour = "#848587",
                                          size = 0.5,
                                          fill = NA),
              legend.text = element_text(size = 12, 
                                         family = "Avenir"),
              axis.line = element_blank())
    
    #plot full retention time
    p2 <- ggplot() +
        geom_line(mapping = aes(rt,
                                intensity,
                                colour = sample,
                                group = sample),
                  data = df2,
                  lwd = 1.2) +
        scale_colour_manual(values = gh99_pal,
                            name = "") +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(sample)) +
        theme_classic() +
        theme(strip.text = element_blank(),
              axis.text = element_text(size = 12,
                                       family = "Avenir"),
              axis.title = element_text(size = 14, 
                                        family = "Avenir LT 65 Medium"),
              panel.border = element_rect(colour = "#848587",
                                          size = 0.5,
                                          fill = NA),
              legend.text = element_text(size = 12, 
                                         family = "Avenir"),
              axis.line = element_blank())
    
    #combine plots
    combined <- p1 + p2 & 
        theme(legend.position = "bottom",
              plot.title = element_text(size = 12, 
                                        family = "Avenir LT 65 Medium",
                                        hjust = 0.5))
    
    #save as pdf
    cairo_pdf(paste0("./analysis/analysis_plots/eic_checking/gh99/",
                     ion, "_mz", mz,
                     ".pdf"),
              width = 12,
              height = 9)
    print(combined +
              plot_layout(ncol=2, guides = "collect") + 
              plot_annotation(title = paste0("ion = ",
                                             ion,
                                             "; m/z = ",
                                             mz)))
    
    dev.off()
}




#18: Screening for PLx digest products -----
PLx.df1 <- peaks_all %>% filter(PLx.digest >= 1 & 
                                    no.digest == 0 & 
                                    mannose.2.sulfate == 0 & 
                                    mannose.3.sulfate == 0 & 
                                    mannose.4.sulfate == 0 & 
                                    mannose.6.sulfate == 0 & 
                                    GH92.digest == 0 & 
                                    GH92...GH99.sulphatase.digest == 0)

PLx.df1 <- cbind(rt_minutes = round(PLx.df1$rt/60, 1),
                 PLx.df1)

PLx.mz <- round(PLx.df1$mz, 3) %>% unique()
PLx.df1$mz_round <- round(PLx.df1$mz, 3)
PLx.df2 <- PLx.df1[PLx.df1$mz_round %in% PLx.mz,] %>% 
    distinct(mz_round, .keep_all = T) %>% 
    dplyr::select(c(annotation, mz_round, rt))
PLx.ion <- PLx.df2$annotation

PLx.rt <- PLx.df2$rt/60

PLx_data <- filterFile(res,
                       file = which(res$group == "PLx digest"))

error = 0.0025
PLx_chrlist <- list()
for (i in 1:length(PLx.mz)){
    mzr = c(PLx.mz[i] - error,
            PLx.mz[i] + error)
    PLx_chrlist[[i]] <- chromatogram(PLx_data, 
                                     mz = mzr)
}

#save list
save(PLx_chrlist,
     file = "./analysis/RData/PLx_chrlist.RData")

#extract intensity and rt values
PLx_int_list <- list()
for (i in 1:length(PLx_data$name)){
    PLx_int_list[[i]] <- lapply(PLx_chrlist, function(x) {
        x[[i]]@intensity
    }) 
}

PLx_rt_list <- list()
for (i in 1:length(PLx_data$name)){
    PLx_rt_list[[i]] <- lapply(PLx_chrlist, function(x) {
        x[[i]]@rtime
    }) 
}

#build data frame (long format)
PLx_chr.df <- data.frame(sample = as.character(),
                         ion = as.character(),
                         mz = as.numeric(),
                         rt = as.numeric(),
                         intensity = as.numeric())

for (i in 1:length(PLx_data$name)){
    for (j in 1:length(PLx.mz)){
        rt = PLx_rt_list[[i]][[j]]/60
        intensity = PLx_int_list[[i]][[j]]
        sample = rep(PLx_data$name[i], length(rt))
        ion = rep(PLx.ion[j], length(rt))
        mz = rep(PLx.mz[j], length(rt))
        temp <- data.frame(sample = sample,
                           ion = ion,
                           mz = mz,
                           rt = rt,
                           intensity = intensity)
        PLx_chr.df <- rbind(PLx_chr.df,
                            temp)
    }
}



#set NA to 0
PLx_chr.df[is.na(PLx_chr.df)] <- 0

dir.create("./analysis/analysis_plots/eic_checking/PLx",
           showWarnings = FALSE)

#make palette
PLx_pal <- c("#F9A100",
             "#46494B")

#plot for each ion the full and restricted rt chromatograms by group
for (i in 1:length(PLx.mz)){
    ion = PLx.ion[i]
    mz = PLx.mz[i]
    rt = PLx.rt[i]
    rtmin = rt - 2.5
    if (rtmin < 0) {rtmin <- 0}
    rtmax = rt + 2.5
    
    df1 <- PLx_chr.df %>% 
        filter(between(rt, rtmin, rtmax)) %>% 
        filter(ion == !!ion)
    
    df2 <- PLx_chr.df %>% 
        filter(ion == !!ion)
    
    #plot restricted retention time
    p1 <- ggplot() +
        geom_line(mapping = aes(rt,
                                intensity,
                                colour = sample,
                                group = sample),
                  data = df1,
                  lwd = 1.2) +
        scale_colour_manual(values = PLx_pal,
                            name = "") +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(sample),
                   scales = "free_y") +
        theme_classic() +
        theme(strip.text = element_blank(),
              axis.text = element_text(size = 12,
                                       family = "Avenir"),
              axis.title = element_text(size = 14, 
                                        family = "Avenir LT 65 Medium"),
              panel.border = element_rect(colour = "#848587",
                                          size = 0.5,
                                          fill = NA),
              legend.text = element_text(size = 12, 
                                         family = "Avenir"),
              axis.line = element_blank())
    
    #plot full retention time
    p2 <- ggplot() +
        geom_line(mapping = aes(rt,
                                intensity,
                                colour = sample,
                                group = sample),
                  data = df2,
                  lwd = 1.2) +
        scale_colour_manual(values = PLx_pal,
                            name = "") +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(sample)) +
        theme_classic() +
        theme(strip.text = element_blank(),
              axis.text = element_text(size = 12,
                                       family = "Avenir"),
              axis.title = element_text(size = 14, 
                                        family = "Avenir LT 65 Medium"),
              panel.border = element_rect(colour = "#848587",
                                          size = 0.5,
                                          fill = NA),
              legend.text = element_text(size = 12, 
                                         family = "Avenir"),
              axis.line = element_blank())
    
    #combine plots
    combined <- p1 + p2 & 
        theme(legend.position = "bottom",
              plot.title = element_text(size = 12, 
                                        family = "Avenir LT 65 Medium",
                                        hjust = 0.5))
    
    #save as pdf
    cairo_pdf(paste0("./analysis/analysis_plots/eic_checking/PLx/",
                     ion, "_mz", mz,
                     ".pdf"),
              width = 12,
              height = 9)
    print(combined +
              plot_layout(ncol=2, guides = "collect") + 
              plot_annotation(title = paste0("ion = ",
                                             ion,
                                             "; m/z = ",
                                             mz)))
    
    dev.off()
}

#mz to check
PLx.mz.check <- c(258.764, 262.759, 283.264, 330.699, 338.688, 392.652, 
                  412.621, 510.967, 533.8, 556.492, 560.818, 628.428,
                  708.353, 560.487, 184.831, 242.79, 258.921, 702.361,
                  )

#16: Compare samples with heatmap ----
#all done without HILIC standard

##ALL PEAKS 
#add unique rownames
setDF(peaks_all)
x <- paste0("rt_",
            round(peaks_all$rt/60,1),
            "min_",
            peaks_all$annotation)
#which(duplicated(x)) #find which ones are duplicated and add mz values
#95 and 96 are the same
#230 and 232 are the same
x[95] <- paste0("rt_",
                round(peaks_all$rt[95]/60,1),
                "min_",
                peaks_all$annotation[95],
                "_mz",
                round(peaks_all$mz[95],3))
x[96] <- paste0("rt_",
                round(peaks_all$rt[96]/60,1),
                "min_",
                peaks_all$annotation[96],
                "_mz",
                round(peaks_all$mz[96],3))
x[230] <- paste0("rt_",
                round(peaks_all$rt[230]/60,1),
                "min_",
                peaks_all$annotation[230],
                "_mz",
                round(peaks_all$mz[230],3))
x[232] <- paste0("rt_",
                round(peaks_all$rt[232]/60,1),
                "min_",
                peaks_all$annotation[232],
                "_mz",
                round(peaks_all$mz[232],3))

rownames(peaks_all) <- x

#subset to only have intensity counts
counts_all <- peaks_all %>% 
    select(sampleColNames)

#set factor level 
group<-factor(pd$group[pd$group!= "HILIC standard"])

#DGEList:Creates a DGEList object from a table of counts 
#(rows=features, columns=samples), 
#group indicator for each column, 
#library size (optional) and a table of feature annotation (optional).
y_n_all <- DGEList(counts=counts_all,
               group=group)
#filterByExpr {edgeR}
#determine which features have sufficiently large counts to be retained for stats
#output is a logical vector
keep_n_all <- filterByExpr(y_n_all)
y_n_all <- y_n_all[keep_n_all,,keep.lib.sizes=FALSE]

#Calculate normalisation factors to scale the raw library sizes
y_n_all <- calcNormFactors(y_n_all)

#creates a design (or model) matrix, e.g., by expanding factors to a 
#set of dummy variables (depending on the contrasts) and 
#expanding interactions similarly.
design <- model.matrix(~group)

#estimate disparity
y_n_all <- estimateDisp(y_n_all,design)
y_n_all <- estimateCommonDisp(y_n_all)

#test difference
#output:
#log2-fold-change (logFC), 
#the average log2-counts-per-million (logCPM), 
#and the two-sided p-value (PValue).
tested_n_all <-exactTest(y_n_all)
hist(tested_n_all$table[,"PValue"], breaks=50)

#extract most different
result_n_all <- topTags(tested_n_all, 
                    n=nrow(tested_n_all$table))

#log transform
dge_n_all <- log2(y_n_all$counts + 1)

#subset
selY_n_all <- dge_n_all[rownames(result_n_all$table)[result_n_all$table$FDR<0.05 & 
                                                         result_n_all$table$logFC < -2 |
                                                         result_n_all$table$FDR<0.05 & 
                                                         result_n_all$table$logFC > 2],]


#make heatmap
cimColour <- rev(viridis(1000))
# cimColurCols <- c(rep("#DA95CD", 3),
#                   rep("#ADA4E2", 3))
# 
# cimColourRows <- rownames(volcanoData_n)[volcanoData_n$diff!="stable"]
# cimColourRows[cimColourRows %in% 
#                   rownames(volcanoData_n)[
#                       volcanoData_n$label_diff=="unknownDOWN"]] <- "#ADA4E2"
# cimColourRows[cimColourRows %in% 
#                   rownames(volcanoData_n)[
#                       volcanoData_n$label_diff=="matchedDOWN"]] <- "#5D3D7F"
# cimColourRows[cimColourRows %in% 
#                   rownames(volcanoData_n)[
#                       volcanoData_n$label_diff=="unknownUP"]] <- "#DA95CD"
# cimColourRows[cimColourRows %in% 
#                   rownames(volcanoData_n)[
#                       volcanoData_n$label_diff=="matchedUP"]] <- "#8D0141"    


# tiff("./analysis/heatmap/heatmap_matched_unmatched.tiff",
#      units = "in",
#      res = 300,
#      width = 12,
#      height = 16)

# svg("./analysis/heatmap/heatmap_matched_unmatched.svg",
#     width = 12,
#     height = 9)

par(mar= c(10, 15, 15, 40))
cim(selY_n_all, 
    color = cimColour,
    symkey=FALSE,
    mar=c(5,20),
    #row.sideColors = cimColourRows,
    row.names = TRUE,
    #row.cex = 0.8,
    #keysize = c(0.1, 0.1)
    #save = "tiff",
    #name.save = "./analysis/heatmap/heatmap_matched_unmatched.tiff"
)
dev.off()

rownames(selY_n) <- rownames(selY_n) %>% 
    sub("_unknown", "", .)


#plot with heatmap.2

tiff("./analysis/heatmap/heatmap_matched_unmatched_2.tiff",
     units = "in",
     res = 300,
     width = 12,
     height = 16)

svg("./analysis/heatmap/heatmap_matched_unmatched_2.svg",
    width = 12,
    height = 16)


par(family = "Avenir")
h<- heatmap.2(selY_n,
              #colour param
              col = cimColour, # color palette defined earlier
              RowSideColors = cimColourRows, # Colour row for rows
              trace = "none", # controls trace lines inside the heat map
              dendrogram = "both",
              density.info = "none",
              rowsep = 1:nrow(selY_n),
              sepcolor = "white",
              sepwidth = c(0.01, 0.01),
              
              
              #colour key param
              key = TRUE, # show the colour key
              key.title = NA,
              key.xlab = "log2(intensity)",
              key.ylab = NA,
              
              #layout param
              margins=c(10,58), # height, width margins around plot
              lwid= c(0.08,0.005, 0.3),
              lhei = c(1,15),
              lmat = rbind(c(5,5,4), c(3,1,2)),
              cexRow = 1.5,# text size rows
              cexCol = 2, # text size cols
              srtCol = 45
              
)

dev.off()




#16: Extract MS2 associated with annotated features in final peak list-----
#get peak info
cp <- chromPeaks(data_ms2)

#format final feature list
peaksFiltered <- peaks_matched
peaksFiltered$mz <- NULL
names(peaksFiltered)[names(peaksFiltered) == "i.mz"] <- "mz"
names(peaksFiltered)[names(peaksFiltered) == "i.mzmin"] <- "mzmin"
names(peaksFiltered)[names(peaksFiltered) == "i.mzmax"] <- "mzmax"
peaksFiltered <- peaksFiltered[,c("mz", "mzmin", "mzmax")]

#filter peaks to match identified features
cp <- as.data.frame(cp)
setDT(cp); setDT(peaksFiltered)
setkey(cp, mzmin, mzmax)
setkey(peaksFiltered, mzmin, mzmax)
cp.matched <- foverlaps(cp,
                        peaksFiltered)
cp.filtered <- cp.matched[cp.matched$mz!=""]

cp.filtered$mz <- NULL
cp.filtered$mzmin <- NULL
cp.filtered$mzmax <- NULL
names(cp.filtered)[names(cp.filtered) == "i.mz"] <- "mz"
names(cp.filtered)[names(cp.filtered) == "i.mzmin"] <- "mzmin"
names(cp.filtered)[names(cp.filtered) == "i.mzmax"] <- "mzmax"

#assign chromPeaks tp ms2 file
cp.new <- as.matrix(cp.filtered)
data_ms2_old <- data_ms2 #keep to be safe
chromPeaks(data_ms2) <- cp.new

#filter file to only contain samples with ms2 data
data_ms2_all <- data_ms2 #keep to be safe
data_ms2 <- filterFile(data_ms2_all,
                       file = which(!grepl("standard|Solvent",
                                           data_ms2_all$name)))

#extract ms2 associated with chromPeaks
ms2_features <- chromPeakSpectra(data_ms2,
                                 msLevel = 2,
                                 expandRt = 2.5,
                                 expandMz = 0.005,
                                 skipFilled = FALSE,
                                 method = "all",
                                 return.type = "Spectra")


#remove zero intensity masses
ms2_features@listData <- lapply(ms2_features, 
                                clean, 
                                all = TRUE)

#combine spectra by peak within each file
ms2_features_comb <- combineSpectra(ms2_features,
                                    fcol = 'peak_id',
                                    mzd = 0.005,
                                    intensityFun = mean)

#see how many spectra per sample
table(fromFile(ms2_features_comb))

#normalise with respect to ion with highest intensity
ms2_features_comb_old <- ms2_features_comb
ms2_features_comb@listData <- lapply(ms2_features_comb_old, 
                                     normalise, 
                                     method = "max")

#see which features have ms2 spectra
precursors <- precursorMz(ms2_features_comb)
rtime <- rtime(ms2_features_comb)
ms2_samples <- fromFile(ms2_features_comb)
names(ms2_samples) <- data_ms2$name[ms2_samples]

precursors.df <- data.frame(peakId = names(precursors),
                            precursorMz = precursors,
                            rtime = rtime,
                            sample = names(ms2_samples),
                            i.mzmin = precursors - 0.0025,
                            i.mzmax = precursors + 0.0025)

setDT(precursors.df)
setDT(peaks_matched)
setkey(precursors.df, i.mzmin, i.mzmax)
sample_peaks_matched_ms2 <- foverlaps(peaks_matched,
                                      precursors.df)
sample_peaks_matched_ms2 <- sample_peaks_matched_ms2 %>% 
    replace_na(list("peakId"="",
                    "precursorMz" = "", 
                    "rtime"= "",
                    "sample" = "",
                    "mzmin" = "", 
                    "mzmax" = ""))
sample_peaks_matched_ms2 <- sample_peaks_matched_ms2[
    sample_peaks_matched_ms2$peakId!=""]

#extract data
ms2.df <- data.frame(sample.name = as.character(),
                     sample.number = as.numeric(),
                     precursorMz = as.numeric(),
                     rt = as.numeric(),
                     mz = as.numeric(),
                     intensity = as.numeric())

for (i in 1:length(ms2_features_comb)){
    mz = sprintf("%.4f",ms2_features_comb[[i]]@mz)
    intensity = ms2_features_comb[[i]]@intensity * 100
    rt = rep(sprintf("%.4f", ms2_features_comb[[i]]@rt / 60),
             length(mz))
    precursorMz = rep(sprintf("%.4f",ms2_features_comb[[i]]@precursorMz),
                      length(mz))
    sample.number = rep(fromFile(ms2_features_comb[[i]]),
                        length(mz))
    sample.name = rep(data_ms2$name[fromFile(ms2_features_comb[[i]])],
                        length(mz))
    temp <- data.frame(sample.number = sample.number,
                       sample.name = sample.name,
                       precursorMz = precursorMz,
                       rt = rt,
                       mz = mz,
                       intensity = intensity)
    ms2.df <- rbind(temp,
                    ms2.df)
}

ms2.df$precursorMz <- as.numeric(ms2.df$precursorMz)
ms2.df$rt <- as.numeric(ms2.df$rt)
ms2.df$mz <- as.numeric(ms2.df$mz)

#17: Screen MS2 associated with features for significant ions -----
sig_ions <- data.frame(mz = c(96.96011,
                              80.96519,
                              259.0129,
                              241.0024,
                              198.9918,
                              180.9812,
                              138.9707),
                       fragment = c("HSO4-",
                                    "HSO3-",
                                    "C6H11O6SO3-",
                                    "C6H9O5SO3-",
                                    "C4H7O4SO3-",
                                    "C4H5O3SO3-",
                                    "C2H3O2SO3-"))

sig_ions$mzmin <- sig_ions$mz - 0.005
sig_ions$mzmax <- sig_ions$mz + 0.005
ms2.df$mzmin <- ms2.df$mz - 0.005
ms2.df$mzmax <- ms2.df$mz + 0.005

setDT(sig_ions); setDT(ms2.df)
setkey(sig_ions, mzmin, mzmax)
ms2.df_matched <- foverlaps(ms2.df,
                            sig_ions)

ms2.df_matched <- ms2.df_matched %>% 
    replace_na(list("mz"="",
                    "fragment" = "", 
                    "mzmin"= "",
                    "mzmax" = ""))
ms2.df_matched <- ms2.df_matched[ms2.df_matched$fragment !=""]

#18: Compute and screen differences within MS2 associated with features ----
ms2_differences.df <- data.frame(sample.name = as.character(),
                                 sample.number = as.numeric(),
                                 precursorMz = as.numeric(),
                                 rt = as.numeric(),
                                 differences = as.numeric())

for (i in 1:length(ms2_features_comb)){
    precursorMz.var = as.numeric(sprintf("%.4f",
                                         ms2_features_comb[[i]]@precursorMz))
    sample.number = fromFile(ms2_features_comb[[i]])
    sample.name = data_ms2$name[sample.number]
    rt = sprintf("%.4f", ms2_features_comb[[i]]@rt / 60)
    x <- as.numeric(sprintf("%.4f",ms2_features_comb[[i]]@mz))
    y <- data.frame(difference = c(dist(x)))
    y <- cbind(sample.name = rep(sample.name,
                                  nrow(y)),
               sample.number = rep(sample.number,
                                   nrow(y)),
               precursorMz = rep(precursorMz.var,
                                 nrow(y)),
               rt = rep(rt,
                        nrow(y)),
               y)
    ms2_differences.df <- rbind(y,
                                ms2_differences.df) 
}

#significant differences:
#18.010565 -> H2O
#60.021130 -> C2H4O2
#120.042260 -> C4H8O4
#162.052824 -> C6H10O5
#79.957 -> SO3

sig_dif <- data.frame(mz = c(18.010565,
                             60.021130,
                             120.042260,
                             162.052824,
                             79.957),
                      loss = c("H2O",
                               "C2H4O2",
                               "C4H8O4",
                               "C6H10O5",
                               "SO3"))

sig_dif$mzmin <- sig_dif$mz - 0.005
sig_dif$mzmax <- sig_dif$mz + 0.005
ms2_differences.df$mzmin <- ms2_differences.df$difference - 0.005
ms2_differences.df$mzmax <- ms2_differences.df$difference + 0.005

setDT(sig_dif); setDT(ms2_differences.df)
setkey(sig_dif, mzmin, mzmax)
ms2_differences.df_matched <- foverlaps(ms2_differences.df,
                                        sig_dif)

ms2_differences.df_matched <- ms2_differences.df_matched %>% 
    replace_na(list("mz"="",
                    "loss" = "", 
                    "mzmin"= "",
                    "mzmax" = ""))
ms2_differences.df_matched <- ms2_differences.df_matched[
    ms2_differences.df_matched$loss !=""]

#19: Screen MS2 associated with all picked peaks for significant differences ----
#filter files
data_ms2_old <- filterFile(data_ms2_old,
                       file = which(!grepl("standard|Solvent",
                                           data_ms2_old$name)))
#extract spectra
ms2_peaks <- chromPeakSpectra(data_ms2_old,
                              msLevel = 2,
                              expandRt = 2.5,
                              expandMz = 0.005,
                              skipFilled = FALSE,
                              method = "all",
                              return.type = "Spectra")

#remove zero intensity masses
ms2_peaks@listData <- lapply(ms2_peaks, 
                             clean, 
                             all = TRUE)

#combine spectra by peak within each file
ms2_peaks_comb <- combineSpectra(ms2_peaks,
                                 fcol = 'peak_id',
                                 mzd = 0.005,
                                 intensityFun = mean)

#see how many spectra per sample
table(fromFile(ms2_peaks_comb))

#normalise with respect to ion with highest intensity
ms2_peaks_comb_old <- ms2_peaks_comb
ms2_peaks_comb@listData <- lapply(ms2_peaks_comb_old, 
                                     normalise, 
                                     method = "max")

#build dataframe of differences
ms2_peaks_differences.df <- data.frame(sample.name = as.character(),
                                       sample.number = as.numeric(),
                                       precursorMz = as.numeric(),
                                       rt = as.numeric(),
                                       differences = as.numeric())

for (i in 1:length(ms2_peaks_comb)){
    precursorMz.var = as.numeric(sprintf("%.4f",
                                         ms2_peaks_comb[[i]]@precursorMz))
    sample.number = fromFile(ms2_peaks_comb[[i]])
    sample.name = data_ms2_old$name[sample.number]
    rt = sprintf("%.4f", ms2_peaks_comb[[i]]@rt / 60)
    k <- ms2_peaks_comb[[i]]@intensity
    x <- as.numeric(sprintf("%.4f",
                            ms2_peaks_comb[[i]]@mz[which(k > 0.05)]))
    y <- data.frame(difference = c(dist(x)))
    y <- cbind(sample.name = rep(sample.name,
                                 nrow(y)),
               sample.number = rep(sample.number,
                                   nrow(y)),
               precursorMz = rep(precursorMz.var,
                                 nrow(y)),
               rt = rep(rt,
                        nrow(y)),
               y)
    ms2_peaks_differences.df <- rbind(y,
                                      ms2_peaks_differences.df) 
}

ms2_peaks_differences.df$mzmin <- ms2_peaks_differences.df$difference - 0.0025
ms2_peaks_differences.df$mzmax <- ms2_peaks_differences.df$difference + 0.0025

setDT(ms2_peaks_differences.df)
ms2_peaks_differences.df_matched <- foverlaps(ms2_peaks_differences.df,
                                              sig_dif)

ms2_peaks_differences.df_matched <- ms2_peaks_differences.df_matched %>% 
    replace_na(list("mz"="",
                    "loss" = "", 
                    "mzmin"= "",
                    "mzmax" = ""))
ms2_peaks_differences.df_matched <- ms2_peaks_differences.df_matched[
    ms2_peaks_differences.df_matched$loss !=""]


#match precursors with predictions
ms2_peaks_differences.df_matched$mzmin <- NULL
ms2_peaks_differences.df_matched$mzmax <- NULL
ms2_peaks_differences.df_matched$i.mzmax <- NULL
ms2_peaks_differences.df_matched$i.mzmin <- NULL

ms2_peaks_differences.df_matched$mzmin <- ms2_peaks_differences.df_matched$precursorMz-0.0025
ms2_peaks_differences.df_matched$mzmax <- ms2_peaks_differences.df_matched$precursorMz+0.0025

ms2_peaks_differences.df_matched_nopred <- ms2_peaks_differences.df_matched
ms2_peaks_differences.df_matched <- foverlaps(ms2_peaks_differences.df_matched,
                                              predicted)

#change NA values created during matching (features with no match) to be blank
#remove extra columns
ms2_peaks_differences.df_matched$mzmin <- NULL
ms2_peaks_differences.df_matched$mzmax <- NULL

ms2_peaks_differences.df_matched <- ms2_peaks_differences.df_matched %>% 
    replace_na(list("name"="",
                    "ion"= "", 
                    "mz" = "",
                    "dp" = ""))
#only keep matched features
ms2_peaks_differences.df_matched_annot <- ms2_peaks_differences.df_matched[!ms2_peaks_differences.df_matched$name=="",]







