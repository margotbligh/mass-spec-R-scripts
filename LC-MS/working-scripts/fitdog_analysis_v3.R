setwd("~/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020")
load("./analysis/neutral/RData/RData_20210223.RData")

# 1: Install packages --------------------------------------------------------
library(BiocStyle)
#get version 1.0.4.6 of Rccp, which many parts of xcms built against
# install.packages("https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.4.6.tar.gz",
#                  repos=NULL,
#                  type="source") 
#DO NOT UPDATE RCCP WHEN ASKED IN OTHER PACKAGE INSTALLATIONS
#BiocManager::install("xcms")
library(xcms)
#BiocManager::install("faahKO")
library(faahKO)
library(pander)
library(RColorBrewer)
library(magrittr)
library(pheatmap)
library(MSnbase)
library(msdata)
library(png)
#BiocManager::install("IPO")
library(IPO)
library(tidyr)
library(detect)
library(devtools)
#install_github('WMBEdmands/MetMSLine')
library(MetMSLine)
#vignette('MetMSLineBasics')
library(pcaMethods)
#library(bisoreg) #no longer on CRAN
#BiocManager::install("statTarget")
library(statTarget)
library(randomForest)
library(rlist)
library(purrr)
library(wesanderson)
library(ggplot2)
library(reshape2)
library(extrafont)
library(Rmisc)
#BiocManager::install("edgeR")
library(edgeR)
library(limma)
#BiocManager::install("mixOmics")
library(mixOmics)
#BiocManager::install("HTSFilter")
library(HTSFilter)
library(rstatix)
library(reshape2)
library(scales)
library(data.table)
#install.packages("remotes")
library(remotes)
# remotes::install_github("lgatto/ProtGenerics")
# remotes::install_github("lgatto/MSnbase")
# remotes::install_github("EuracBiomedicalResearch/CompMetaboTools")
library(ggridges)
library(gridExtra)
library(tidyverse)
library(ggpubr)
library(viridis)
library(lemon)
library(gtable)
library(grid)
library(colorspace)

#2. Import and inspect data --------------------------------------------------------

#set working directory
setwd("~/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020")

#get file paths to mzML files
neutral_fp <- dir(path = "./neutral_fitdog/mzML-files", 
                  all.files = FALSE, 
                  full.names = TRUE)
anionic_fp <- dir(path = "./anionic_fitdog/mzML-files", 
                  all.files = FALSE, 
                  full.names = TRUE)
control_fp <- dir(path = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_mannan/orbitrap/data/controls/mzML-files/20201223", 
                  all.files = FALSE, 
                  full.names = TRUE)
control_fp <- c(control_fp,
                dir(path = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_mannan/orbitrap/data/controls/mzML-files/20201221", 
                    all.files = FALSE, 
                    full.names = TRUE))

control_fp <- control_fp[grep("neg|22_P1000|_23|_16|_12",
                              control_fp)]
#controls:
#solvent blank 16 from 20201222
#solvent blank 12 from 20201221
#standard mix
#blank acetone precipitation
#blank procainamide labelling


fp <- c(anionic_fp,
        neutral_fp,
        control_fp)

rm(anionic_fp,
   neutral_fp,
   control_fp)
#create phenodata data.frame
#every sample must have a unique name!

name <- basename(fp) %>%
    sub("MS31_20201222_|MS31_20201221_", "", .) %>% 
    sub("-procA", "", .) %>% 
    sub("-rep", "", .) %>% 
    sub("pos_", "", .) %>% 
    sub("SolventBlank_7_12", "solvent blank 1",.) %>% 
    sub("SolventBlank_16", "solvent blank 2",.) %>% 
    sub("P1000.*", "standard mix", .) %>% 
    sub("neg_", "blank ", .) %>% 
    sub("_\\d\\d.mzML|.mzML|_\\d_\\d\\d.mzML", "", .) %>% 
    sub("lambda-carra", "lambda carrageenan", .) %>% 
    sub("yeastmannan", "yeast mannan", .) %>% 
    sub("lam_", "laminarin_", .) %>%
    sub("_fitdog", " fitdog", .) %>%
    sub("_omix", " standard mix", .) %>%
    sub("_nodigest", " undigested", .) %>%
    sub("_gh|-gh", " + GH", .) %>%
    gsub("-|_", " ", .) %>%
    sub("1$", " #1", .) %>%
    sub("2$", " #2", .) %>%
    gsub("\\s\\s", " ", .)

sample_type <- basename(fp) %>% 
    sub(".*fitdog.*", "FITDOG", .) %>% 
    sub(".*heat-inact.*|.*nodigest.*", "undigested", .) %>% 
    sub(".*gh.*", "enzyme digest", .) %>% 
    sub(".*omix.*|.*P1000.*", "standard mix", .) %>% 
    sub(".*neg_.*|.*SolventBlank.*", "neg", .)

substrate <- basename(fp) %>% 
    sub(".*lambda.*", "lambda carrageenan", .) %>% 
    sub(".*lam_.*", "laminarin", .) %>% 
    sub(".*fucoidan.*", "fucoidan", .) %>% 
    sub(".*yeastmannan.*", "yeast mannan", .) %>% 
    sub(".*P1000.*", "mixed", .) %>% 
    sub(".*neg.*|.*Blank.*", "blank", ., ignore.case = TRUE)

replicate <- name %>% 
    sub(".*#", "", .) %>% 
    sub("\\D+.*", "NA", .)

method <- fp %>% 
    sub(".*neutral.*MS31_20201222.*|.*SolventBlank_16.*", 
        "pos mode inclusion list MSMS", .) %>% 
    sub(".*neutral.*MS31_20201221.*|.*anionic.*|.*SolventBlank_7_12.*|.*neg.*", 
        "neg and pos mode MS1 only", .) %>% 
    sub(".*P1000.*", 
        "pos mode inclusion list MSMS", .)

pd <- data.frame(name = name,
                 sample_type = sample_type,
                 substrate = substrate,
                 replicate = replicate,
                 method = method,
                 stringsAsFactors = FALSE)

#load data
data <- readMSData(files = fp, 
                   pdata = new("NAnnotatedDataFrame", 
                               pd), 
                   mode = "onDisk")

#separate positive and negative mode spectra
pdata <- data[data@featureData@data$polarity == 1]
ndata <- data[data@featureData@data$polarity == 0]

#keep only MS1 data
pdata <- pdata[pdata@featureData@data$msLevel == 1]
ndata <- ndata[ndata@featureData@data$msLevel == 1]

#see number of spectra for each MS level
table(msLevel(pdata))
table(msLevel(ndata))

#import mz values predicted using my script with the command:
#sugarMassesPredict-procainamide.py -dp 1 10 -p 0 -m none -n 1 -ld D -oh 0 -b 0 -i pos -s 250 2000 -o neutral-masses-predicted.txt -l procainamide

neutral.mz.ptable <- fread("./analysis/neutral/neutral-masses-predicted.txt")
neutral.mz.pvec <- neutral.mz.ptable$`[M+H]`
neutral.name.pvec <- neutral.mz.ptable$name

#3: Create initial output directories -------------------------------------
dir.create("./analysis/all_samples/RData",
           showWarnings = FALSE)
dir.create("./analysis/all_samples/processing_plots",
           showWarnings = FALSE)
dir.create("./analysis/all_samples/analysis_plots",
           showWarnings = FALSE)
dir.create("./analysis/all_samples/processing_tables",
           showWarnings = FALSE)
dir.create("./analysis/all_samples/analysis_tables",
           showWarnings = FALSE)

#4: Set phenodata colours and graphic param -------------------------------------
par(family = "Avenir")

#set colours for samples
#set same colours for replicates
pal_warm <- hcl.colors(5, "Warm")
pal_cold <- hcl.colors(5, "Cold")

pal_sample <- c(pal_warm,
                pal_cold,
                pal_cold[1:4],
                "#C6EDD0",
                "#E2D3A9",
                "#DFAF3B",
                "#C5C6C7",
                "#C5C6C7")
names(pal_sample) <- name

#show colours
show_col(pal_sample)

#set colours for sample types
pal_type <- hcl.colors(n = 5, palette = "Dynamic")
names(pal_type) <- unique(sample_type)

#show colours
show_col(pal_type)

#set colours for substrate
pal_substrate <- c(pal_type,
                   "#DFAF3B")
names(pal_substrate) <- unique(substrate)

#show colours
show_col(pal_substrate)

#set colours for methods
pal_method <- hcl.colors(5, "Pastel 1")
pal_method <- pal_method[c(2,5)]
names(pal_method) <- unique(pd$method)

#show colours
show_col(pal_method)

#5: Plot BPC and TIC ---------------------------
dir.create("./analysis/all_samples/processing_plots/bpc_tic",
           showWarnings = FALSE)

#plot the tic as boxplot
tc <- split(tic(pdata), 
            f = fromFile(pdata))
cairo_pdf("./analysis/all_samples/processing_plots/bpc_tic/pos_tic_boxplot.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mar=c(15,5,1,1))
boxplot(tc, 
        col = pal_sample[pdata$name],
        ylab = "intensity", 
        main = "total ion current",
        names = pdata$name,
        las=2,
        cex.axis = 0.8)
dev.off()

tc <- split(tic(ndata), 
            f = fromFile(ndata))
cairo_pdf("./analysis/all_samples/processing_plots/bpc_tic/neg_tic_boxplot.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mar=c(15,5,1,1))
boxplot(tc, 
        col = pal_sample[ndata$name],
        ylab = "intensity", 
        main = "total ion current",
        names = ndata$name,
        las=2,
        cex.axis = 0.8)
dev.off()

#extract BPC
bpi.pos <- chromatogram(pdata,
                        aggregationFun = "max")
bpi.neg <- chromatogram(ndata,
                        aggregationFun = "max")
#bin the BPC
bpi.pos_bin <- bin(bpi.pos,
                   binSize = 2)

bpi.neg_bin <- bin(bpi.neg,
                   binSize = 2)

#calculate correlation on the log2 transformed base peak intensities
#for some weird reason for pos mode only works if you remove last row... 
#some samples have one extra spectrum? shows up as --inf for all other samples
x<- log2(do.call(cbind, 
                 lapply(bpi.pos_bin, 
                        intensity)))
x<- x[1:nrow(x)-1, ]
cormat.pos <- cor(x)
colnames(cormat.pos) <- rownames(cormat.pos) <- pdata$name

x<- log2(do.call(cbind, 
                 lapply(bpi.neg_bin, 
                        intensity)))
cormat.neg <- cor(x)
colnames(cormat.neg) <- rownames(cormat.neg) <- ndata$name

#define which phenodata columns should be highlighted in the plot
ann <- data.frame(method = pdata$method)
rownames(ann) <- pdata$name

#perform the cluster analysis for positive 
cairo_pdf("./analysis/all_samples/processing_plots/bpc_tic/bpc.pos_heatmap_bymethod.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
pheatmap(cormat.pos, 
         annotation = ann,
         annotation_colors = list(method = pal_method),
         color = rev(hcl.colors(n=100,
                                palette = "Viridis")),
)
dev.off()

#define which phenodata columns should be highlighted in the plot
ann <- data.frame(method = pdata$substrate)
rownames(ann) <- pdata$name

#perform the cluster analysis for negative 
cairo_pdf("./analysis/all_samples/processing_plots/bpc_tic/bpc.neg_heatmap_bysubstrate.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
pheatmap(cormat.neg, 
         annotation = ann,
         annotation_colors = list(substrate = pal_substrate),
         color = rev(hcl.colors(n=100,
                                palette = "Viridis")),
)
dev.off()

#6: Peak picking (CentWave) ---------------------------
dir.create("./analysis/all_samples/processing_plots/peakpicking",
           showWarnings = FALSE)

cwp_default <- CentWaveParam() # for reference

# Object of class:  CentWaveParam 
# Parameters:
# ppm: 25 
# peakwidth: 20, 50 
# snthresh: 10 
# prefilter: 3, 100 
# mzCenterFun: wMean 
# integrate: 1 
# mzdiff: -0.001 
# fitgauss: FALSE 
# noise: 0 
# verboseColumns: FALSE 
# roiList length: 0 
# firstBaselineCheck TRUE 
# roiScales length: 0 

cwp<-CentWaveParam()
cwp@ppm<-1.6
cwp@peakwidth<-c(10,100)
cwp@snthresh<-5

#pick peaks for positive data
pdata_peaks<-findChromPeaks(pdata, 
                            param=cwp)

#check with chromatograms
error = 0.0025
chr2 <- chromatogram(pdata_peaks,
                     rt = c(770, 900),
                     mz = c(neutral.mz.pvec[2] - error,
                            neutral.mz.pvec[2] + error))

cairo_pdf("./analysis/all_samples/processing_plots/peakpicking/pos_1.6ppm_5sn_pw10to100_dp2.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mfrow=c(4,5),
    mar = c(3,2,3,2))
for (i in 1:length(name)){
    plot(chr2[,i],
         col = pal_sample[i],
         lwd = 2,
         main = pdata_peaks$name[i],
         cex.main = 1)
}
dev.off()

chr4 <- chromatogram(pdata_peaks,
                     rt = c(1000, 
                            1400),
                     mz = c(neutral.mz.pvec[4] - error,
                            neutral.mz.pvec[4] + error))

cairo_pdf("./analysis/all_samples/processing_plots/peakpicking/pos_1.6ppm_5sn_pw10to100_dp4.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mfrow=c(4,5),
    mar = c(3,2,3,2))
for (i in 1:length(name)){
    plot(chr4[,i],
         col = pal_sample[i],
         lwd = 2,
         main = pdata_peaks$name[i],
         cex.main = 1)
}
dev.off()


#pick peaks for negative data
ndata_peaks<-findChromPeaks(ndata, 
                            param=cwp)


#check with chromatograms
chr1_neg <- chromatogram(ndata_peaks,
                         rt = c(400,
                                700),
                         mz = c(398.229112 - error,
                                398.229112 + error))

cairo_pdf("./analysis/all_samples/processing_plots/peakpicking/neg_1.6ppm_5sn_pw10to100_dp1.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mfrow=c(5,4),
    mar = c(3,2,3,2))
for (i in 1:length(name)){
    plot(chr1_neg[,i],
         col = pal_sample[i],
         lwd = 2,
         main = ndata_peaks$name[i],
         cex.main = 1)
}
dev.off()


cairo_pdf("./analysis/all_samples/processing_plots/peakpicking/neg_1.6ppm_5sn_pw10to100_dp1_xic.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
ndata_peaks %>% 
    filterMz(mz = c(398.229112 - 0.01,
                    398.229112 + 0.01)) %>% 
    filterRt(rt = c(550,
                    600)) %>% 
    plot(type = "XIC")
dev.off()


chr3_neg <- chromatogram(ndata_peaks,
                         rt = c(900,
                                1100),
                         mz = c(722.334762 - error,
                                722.334762 + error))

cairo_pdf("./analysis/all_samples/processing_plots/peakpicking/neg_1.6ppm_5sn_pw10to100_dp3.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mfrow=c(5,4),
    mar = c(3,2,3,2))
for (i in 1:length(name)){
    plot(chr3_neg[,i],
         col = pal_sample[i],
         lwd = 2,
         main = ndata_peaks$name[i],
         cex.main = 1)
}
dev.off()

save(pdata, 
     file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/all_samples/RData/pdata.RData")
save(pdata_peaks, 
     file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/all_samples/RData/pdata_peaks.RData")

save(ndata, 
     file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/all_samples/RData/ndata.RData")
save(ndata_peaks, 
     file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/all_samples/RData/ndata_peaks.RData")


#7: group peaks to create "features"---------
# parameters for grouping in PeakDensityParam
# - sampleGroups
# - bw : defines smoothness of density function (default = 30)
# - binSize : m/z width of the bin/slice in which peaks are grouped (default = 0.25)
# - maxFeatures : max number of features to be defined in one bin (default = 50)
# - minFraction : min proportion of samples in at least one group which a peak has to be present (default = 0.5)
# - minSamples: min number of samples in at least one group which a peak has to be present (default = 1)
# 
# - "dry run" correspondence to check settings
# - extract bpc for specific mz and rt range using "chromatogram"
# - run: plotChromPeakDensity(bpc, PeakDensityParam, col = sample group or name colours)
# --> upper panel = chromatogram with the identified peaks
# --> lower panel 
# points = rt of identified peaks (x) per sample (y)
# black line = distribution along x axis
# grey rectangle = feature 

#parameters
pdp <- PeakDensityParam(sampleGroups = pdata$sample_type,
                        binSize = 0.005,
                        bw = 6) 


# extract and plot chromatograms to test settings
dir.create("./analysis/all_samples/processing_plots/peakgrouping",
           showWarnings = FALSE)

mzr = c(neutral.mz.pvec[2] - error,
        neutral.mz.pvec[2] + error)
chr_mzr562 <- chromatogram(pdata_peaks, 
                           mz = mzr)
cairo_pdf(paste0("./analysis/all_samples/processing_plots/peakgrouping/pos_mz562_bw6.pdf"),
          family = "Avenir",
          width = 12,
          height = 9)
par(mar=c(4,4,3,2))
plotChromPeakDensity(chr_mzr562, 
                     col = pal_sample, 
                     param = pdp,
                     peakBg = pal_sample[chromPeaks(chr_mzr562)[, "sample"]],
                     peakCol = pal_sample[chromPeaks(chr_mzr562)[, "sample"]],
                     peakPch = 16)
legend("topright",
       legend = paste0(seq(1,length(pal_sample),1),
                       "=",
                       names(pal_sample)),
       fill = pal_sample,
       pt.cex = 0.3,
       cex = 0.5,
       bty = "n",
       horiz = FALSE,
       xpd=TRUE)
dev.off()

mzr = c(neutral.mz.pvec[4] - error,
        neutral.mz.pvec[4] + error)
chr_mzr886 <- chromatogram(pdata_peaks, 
                           mz = mzr)
cairo_pdf(paste0("./analysis/all_samples/processing_plots/peakgrouping/pos_mz886_bw6.pdf"),
          family = "Avenir",
          width = 12,
          height = 9)
par(mar=c(4,4,3,2))
plotChromPeakDensity(chr_mzr886, 
                     col = pal_sample, 
                     param = pdp,
                     peakBg = pal_sample[chromPeaks(chr_mzr886)[, "sample"]],
                     peakCol = pal_sample[chromPeaks(chr_mzr886)[, "sample"]],
                     peakPch = 16)
legend("topright",
       legend = paste0(seq(1,length(pal_sample),1),
                       "=",
                       names(pal_sample)),
       fill = pal_sample,
       pt.cex = 0.3,
       cex = 0.5,
       bty = "n",
       horiz = FALSE,
       xpd=TRUE)
dev.off()


#group peaks
pdata_peaks_grouped <- groupChromPeaks(pdata_peaks, param = pdp)

pdp_neg <- PeakDensityParam(sampleGroups = ndata_peaks$sample_type,
                            binSize = 0.005,
                            bw = 6) 

ndata_peaks_grouped <- groupChromPeaks(ndata_peaks, param = pdp_neg)


save(pdata_peaks_grouped, 
     file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/all_samples/RData/pdata_peaks_grouped.RData")

save(ndata_peaks_grouped, 
     file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/all_samples/RData/ndata_peaks_grouped.RData")


#8: fill in missing peaks----------
fpp <- FillChromPeaksParam() #default param very conservative, use them
pdata_peaks_grouped_filled <- fillChromPeaks(pdata_peaks_grouped)
ndata_peaks_grouped_filled <- fillChromPeaks(ndata_peaks_grouped)

save(pdata_peaks_grouped_filled, 
     file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/all_samples/RData/pdata_peaks_grouped_filled.RData")

save(ndata_peaks_grouped_filled, 
     file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/all_samples/RData/ndata_peaks_grouped_filled.RData")

#9: plot PCA --------------
dir.create("./analysis/all_samples/processing_plots/pca",
           showWarnings = FALSE)

#POSITIVE, GROUPED BUT NOT FILLED

#get feat intensities
ft_ints.pos<-log2(featureValues(pdata_peaks_grouped, value = "into"))

##Replace all NA values with zero
ft_ints.pos[is.na(ft_ints.pos)] <- 0

## Perform the PCA with the intensities  mean centered.
pc <- prcomp(t(ft_ints.pos), center = TRUE)

## Plot the PCA
pal_method_long <- pdata_peaks_grouped$method
pal_method_long <- pal_method_long %>% 
    sub(unique(pdata_peaks_grouped$method)[1], pal_method[1], .) %>% 
    sub(unique(pdata_peaks_grouped$method)[2], pal_method[2], .)


cairo_pdf("./analysis/all_samples/processing_plots/pca/posgroupedpeaks_bymethod.pdf",
          family = "Avenir",
          width = 12,
          height = 9)

pcSummary <- summary(pc)
par(family = "Avenir")
plot(pc$x[, 1], 
     pc$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "darkgrey", 
     bg = pal_method_long, 
     cex = 2)
grid()
text(pc$x[, 1], 
     pc$x[,2], 
     labels = pdata_peaks_grouped$name, 
     col = "black",
     pos = 1, 
     cex = 0.5)
dev.off()

#RESULT:  
#not great
#mainly clusters by method...
#however, within each cluster, it looks okayyyy

#POSITIVE, GROUPED and FILLED

#get feat intensities
ft_ints.posfilled<-log2(featureValues(pdata_peaks_grouped_filled, 
                                      value = "into"))

##Replace all NA values with zero
ft_ints.posfilled[is.na(ft_ints.posfilled)] <- 0

## Perform the PCA with the intensities  mean centered.
pc <- prcomp(t(ft_ints.posfilled), center = TRUE)

## Plot the PCA
cairo_pdf("./analysis/all_samples/processing_plots/pca/posgroupedfilledpeaks_bysample.pdf",
          family = "Avenir",
          width = 12,
          height = 9)

pcSummary <- summary(pc)
par(family = "Avenir")
plot(pc$x[, 1], 
     pc$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "darkgrey", 
     bg = pal_sample, 
     cex = 2)
grid()
text(pc$x[, 1], 
     pc$x[,2], 
     labels = pdata_peaks_grouped$name, 
     col = "black",
     pos = 1, 
     cex = 0.5)
dev.off()
#RESULT:
#unsure if this is better.....
#continue with filled peaks for now but may change

pdata <- pdata_peaks_grouped_filled


#NEGATIVE, GROUPED BUT NOT FILLED

#get feat intensities
ft_ints.neg<-log2(featureValues(ndata_peaks_grouped, value = "into"))

##Replace all NA values with zero
ft_ints.neg[is.na(ft_ints.neg)] <- 0

## Perform the PCA with the intensities  mean centered.
pc <- prcomp(t(ft_ints.neg), center = TRUE)

#plot
cairo_pdf("./analysis/all_samples/processing_plots/pca/neggroupedpeaks_bysample.pdf",
          family = "Avenir",
          width = 12,
          height = 9)

pcSummary <- summary(pc)
par(family = "Avenir")
plot(pc$x[, 1], 
     pc$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "darkgrey", 
     bg = pal_sample[ndata_peaks_grouped$name], 
     cex = 2)
grid()
text(pc$x[, 1], 
     pc$x[,2], 
     labels = ndata_peaks_grouped$name, 
     col = "black",
     pos = 1, 
     cex = 0.5)
dev.off()

#RESULT:  
#looks kind of okay?
#lambda carageenan and yeast mannan enzyme reactions, and lam standard mix very different
#pretty much all other samples cluster together

#NEGATIVE, GROUPED and FILLED

#get feat intensities
ft_ints.negfilled<-log2(featureValues(ndata_peaks_grouped_filled, 
                                      value = "into"))

##Replace all NA values with zero
ft_ints.negfilled[is.na(ft_ints.negfilled)] <- 0

## Perform the PCA with the intensities  mean centered.
pc <- prcomp(t(ft_ints.negfilled), center = TRUE)

## Plot the PCA
cairo_pdf("./analysis/all_samples/processing_plots/pca/neggroupedfilledpeaks_bysample.pdf",
          family = "Avenir",
          width = 12,
          height = 9)

pcSummary <- summary(pc)
par(family = "Avenir")
plot(pc$x[, 1], 
     pc$x[,2], 
     pch = 21, 
     main = "",
     xlab = paste0("PC1: ", 
                   format(pcSummary$importance[2, 1] * 100,
                          digits = 3), 
                   " % variance"),
     ylab = paste0("PC2: ", 
                   format(pcSummary$importance[2, 2] * 100,
                          digits = 3), 
                   " % variance"),
     col = "darkgrey", 
     bg = pal_sample[ndata_peaks_grouped_filled$name], 
     cex = 2)
grid()
text(pc$x[, 1], 
     pc$x[,2], 
     labels = ndata_peaks_grouped$name, 
     col = "black",
     pos = 1, 
     cex = 0.5)
dev.off()
#RESULT:
#maybe looks better?
#blanks now quite separate
#samples with anionic and neutral substrates cluster separately

#again, continue with filled peaks


ndata <- ndata_peaks_grouped_filled

#10: save diffreport of xdata -----
xset.p <- as(pdata, "xcmsSet")
xset.n <- as(ndata, "xcmsSet")

sampnames(xset.p) <- pData(pdata)$name
sampclass(xset.p) <- pData(pdata)$sample_type

sampnames(xset.n) <- pData(ndata)$name
sampclass(xset.n) <- pData(ndata)$sample_type

#11. Isotope picking and filtering ----
##create xsannotate object
#extracts the peaktable from a provided xcmsSet,
#which is used for all further analysis

an.p <- xsAnnotate(xset.p)
an.n <- xsAnnotate(xset.n)

##Group peaks of a xsAnnotate object according to their retention time 
#Perfwhm = parameter defines the window width, which is used for matching
anF.p <- groupFWHM(an.p, 
                   perfwhm = 0.6)
anF.n <- groupFWHM(an.n, 
                   perfwhm = 0.6)

##Annotate isotope peaks
#Mzabs = the allowed m/z error
anI.p <- findIsotopes(anF.p, 
                      mzabs=0.01)

anI.n <- findIsotopes(anF.n, 
                      mzabs=0.01)

##Peak grouping after correlation information into pseudospectrum groups 
#cor_eic_th = correlation threshold for EIC correlation
anIC.p <- groupCorr(anI.p, 
                    cor_eic_th=0.75)
anIC.n <- groupCorr(anI.n, 
                    cor_eic_th=0.75)

##Find adducts
anFA.p <- findAdducts(anIC.p, 
                      polarity="positive")

anFA.n <- findAdducts(anIC.n, 
                      polarity="negative")

write.csv(getPeaklist(anFA.p), 
          file="./analysis/all_samples/processing_tables/pos-peaklist_xsannotate.csv")
pl.p <-getPeaklist(anFA.p)

write.csv(getPeaklist(anFA.n), 
          file="./analysis/all_samples/processing_tables/neg-peaklist_xsannotate.csv")
pl.n <-getPeaklist(anFA.n)


#filter by retention time 
#remove everything before 5 min and everything after 35 min
pl_rt.p <- pl.p %>% 
    filter(between(rt,
                   300,
                   2100))

pl_rt.n <- pl.n %>% 
    filter(between(rt,
                   300,
                   2100))

#Filter for isotopes
pl_rt_iso.p <-pl_rt.p[pl_rt.p$isotopes!="",]
pl_rt_iso.n <-pl_rt.n[pl_rt.n$isotopes!="",]


#Filter by blank exclusion
pl_rt_iso_be.p <-pl_rt_iso.p[pl_rt_iso.p$neg==0,]
pl_rt_iso_be.n <-pl_rt_iso.n[pl_rt_iso.n$neg==0,]

pl_rt_be.p <-pl_rt.p[pl_rt.p$neg==0,]
pl_rt_be.n <-pl_rt.n[pl_rt.n$neg==0,]

#make rownames from rt and mz of features
rownames(pl_rt_iso_be.p)<-paste(round(pl_rt_iso_be.p$rt,1),
                                round(pl_rt_iso_be.p$mz,3),
                                sep="_")
rownames(pl_rt_iso_be.n)<-paste(round(pl_rt_iso_be.n$rt,1),
                                round(pl_rt_iso_be.n$mz,3),
                                sep="_")


rownames(pl_rt_be.p)<-paste(round(pl_rt_be.p$rt,1),
                            round(pl_rt_be.p$mz,3),
                            sep="_")
rownames(pl_rt_be.n)<-paste(round(pl_rt_be.n$rt,1),
                            round(pl_rt_be.n$mz,3),
                            sep="_")

#filter for isotopes or adducts
pl_rt_be_isoadd.p <- pl_rt_be.p[pl_rt_be.p$isotopes!=""| 
                                    pl_rt_be.p$adduct!="",]
pl_rt_be_isoadd.n <- pl_rt_be.n[pl_rt_be.n$isotopes!=""| 
                                    pl_rt_be.n$adduct!="",]


#change NA to 0
pl_rt_be_isoadd.p[is.na(pl_rt_be_isoadd.p)] <- 0
pl_rt_be_isoadd.n[is.na(pl_rt_be_isoadd.n)] <- 0

pl_rt_iso_be.p[is.na(pl_rt_iso_be.p)] <- 0
pl_rt_iso_be.n[is.na(pl_rt_iso_be.n)] <- 0

#declutter enviro----
rm(list = ls(pattern = "an.*"))
rm(list = ls(pattern = "bpi.*"))
rm(list = ls(pattern = "chr.*"))
rm(list = ls(pattern = "xset.*"))
rm(list = ls(pattern = "ft_ints.*"))
rm(list = ls(pattern = "[np]data_.*"))

#12: heatmap -------
load("./analysis/all_samples/RData/RData_20210105_20hr.RData")

dir.create("./analysis/all_samples/analysis_plots/heatmap",
           showWarnings = FALSE)

#all samples ----
#subset to only have intensity counts
mask.p <- pdata$name %>% 
    gsub("\\s", ".", .) %>% 
    gsub("#|\\+", ".", .)

counts.p <- pl_rt_be_isoadd.p[,names(pl_rt_be_isoadd.p) %in% mask.p]

mask.n <- ndata$name %>% 
    gsub("\\s", ".", .) %>% 
    gsub("#|\\+", ".", .)

counts.n <- pl_rt_be_isoadd.n[,names(pl_rt_be_isoadd.n) %in% mask.n]

#set factor level 
group.p <-factor(pdata$sample_type)
group.n <-factor(ndata$sample_type)

#DGEList:Creates a DGEList object from a table of counts 
#(rows=features, columns=samples), 
#group indicator for each column, 
#library size (optional) and a table of feature annotation (optional).

y_n.p <- DGEList(counts=counts.p,
                 group=group.p)

y_n.n <- DGEList(counts=counts.n,
                 group=group.n)
#filterByExpr {edgeR}
#determine which features have sufficiently large counts to be retained for stats
#output is a logical vectir
keep_n.p <- filterByExpr(y_n.p)
keep_n.n <- filterByExpr(y_n.pn)

y_n.p <- y_n.p[keep_n.p,,keep.lib.sizes=FALSE]
y_n.n <- y_n.n[keep_n.n,,keep.lib.sizes=FALSE]

#Calculate normalisation factors to scale the raw library sizes
y_n.p <- calcNormFactors(y_n.p)
y_n.n <- calcNormFactors(y_n.n)

#creates a design (or model) matrix, e.g., by expanding factors to a 
#set of dummy variables (depending on the contrasts) and 
#expanding interactions similarly.
design.p <- model.matrix(~group.p)
design.n <- model.matrix(~group.n)

#estimate disparity
y_n.p <- estimateDisp(y_n.p,design.p)
y_n.n <- estimateDisp(y_n.n,design.n)

y_n.p<-estimateCommonDisp(y_n.p)
y_n.n<-estimateCommonDisp(y_n.n)

#test difference
#output:
#log2-fold-change (logFC), 
#the average log2-counts-per-million (logCPM), 
#and the two-sided p-value (PValue).
tested_n.p<-exactTest(y_n.p)
cairo_pdf("./analysis/all_samples/analysis_plots/heatmap/pos-pvalue_hist.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
hist(tested_n.p$table[,"PValue"], breaks=50)
dev.off()

tested_n.n<-exactTest(y_n.n)
cairo_pdf("./analysis/all_samples/analysis_plots/heatmap/neg-pvalue_hist.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
hist(tested_n.n$table[,"PValue"], breaks=50)
dev.off()

#don't have enough replicates in each groupn for any statistics
#to be meaningful...

#extract most different
result_n.p <- topTags(tested_n.p, 
                      n=nrow(tested_n.p$table))
result_n.n <- topTags(tested_n.n, 
                      n=nrow(tested_n.n$table))

#plot volcano plot
volcanoData_n.p <- cbind(result_n.p$table$logFC, 
                         -log10(result_n.p$table$FDR))
colnames(volcanoData_n.p) <- c("logFC", "negLogPval")

volcanoData_n.n <- cbind(result_n.n$table$logFC, 
                         -log10(result_n.n$table$FDR))
colnames(volcanoData_n.n) <- c("logFC", "negLogPval")

plot(volcanoData_n.p, pch=19)
abline(h=-log10(0.5), col="red") #add line for cutoff to include in plot
text(-3,-log10(0.49), "-log10(0.5)", col = "red")

plot(volcanoData_n.n, pch=19) 

#again, don't have enough replicates in each group for any statistics
#to be meaningful...


#counts per million and log2 (normalise and transform)
dge_n.p <- cpm(y_n.p, 
               log=TRUE, 
               prior.count = 1)
dge_n.n <- cpm(y_n.n, 
               log=TRUE, 
               prior.count = 1)


#subset positive. not possible to give a cutoff for negative, one sample per group...
selY_n.p <- dge_n.p[rownames(result_n.p$table)[result_n.p$table$FDR<=0.5],]
selY_n.n <- dge_n.n[rownames(result_n.n$table), ]

#make sample names nice
colnames(selY_n.p) <- colnames(selY_n.p) %>% 
    gsub("\\.{3}", " + ", .) %>% 
    gsub("\\.{2}", " #", .) %>% 
    sub("fitdog", "FITDOG", .) %>% 
    gsub("\\.", " ", .) %>% 
    sub("heat.*", "+ HI-GH82", .) %>% 
    sub("GH8\\s#2", "GH82", .)

colnames(selY_n.n) <- colnames(selY_n.n) %>% 
    gsub("\\.{3}", " + ", .) %>% 
    gsub("\\.{2}", " #", .) %>% 
    sub("fitdog", "FITDOG", .) %>% 
    gsub("\\.", " ", .) %>% 
    sub("heat.*", "+ HI-GH82", .) %>% 
    sub("GH8\\s#2", "GH82", .)


#make heatmap
cimColour <- rev(viridis(1000))

# cimColurRows.all <- c("#86B4CC",
#                       "#C95486",
#                       "#B9AACC",
#                       "#C95486",
#                       rep("#B9AACC",3),
#                       rep("#769AA3",3))
# 
# cimColurRows.all <- c(rep("#C95486",2),
#                       rep("#B9AACC",3),
#                       rep("#769AA3",2),
#                       "#DC6B1E",
#                       rep("#769AA3",2))


# cairo_pdf("./analysis/analysis_plots/heatmap/heatmap.pdf",
#           family = "Avenir",
#           width = 12,
#           height = 9)
cim(t(selY_n.p), 
    color = cimColour,
    symkey=FALSE,
    mar=c(7,12),
    #row.sideColors = cimColurRows.all
)

dev.off()

cim(t(selY_n.n), 
    color = cimColour,
    symkey=FALSE,
    mar=c(7,12),
    #row.sideColors = cimColurRows.all
)

#NEUTRAL SAMPLES ONLY -------
#1: load data and redefine sample type----
load(file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/all_samples/RData/pdata_peaks.RData")
neutral <- pdata_peaks
rm(pdata_peaks)
pd_neutral <- pd
neutral_groups <- pdata$name %>% 
    sub("laminarin\\sfitdog.*|.*mannan\\sfitdog.*",
        "neutral.fitdog", .) %>% 
    sub("laminarin\\sstandard.*|.*GH76.*",
        "neutral.pos", .) %>% 
    sub("fucoidan.*|lambda.*",
        "anionic", .) %>% 
    sub("blank\\sfitdog.*",
        "fitdog.neg", .) %>% 
    sub(".*blank.*",
        "neg", .) %>% 
    sub("standard.*",
        "pos", .) %>% 
    sub("yeast.*",
        "neutral.neg", .)
pd_neutral$sample_type <- neutral_groups
phenoData(neutral) <- new("NAnnotatedDataFrame", 
                          pd_neutral)

save(neutral,
     file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/neutral/RData/neutral.RData")

#2: group peaks to create "features"---------
#parameters - same as before
pdp <- PeakDensityParam(sampleGroups = neutral$sample_type,
                        binSize = 0.005,
                        bw = 6) 
#group peaks
neutral_grouped <- groupChromPeaks(neutral, param = pdp)

save(neutral_grouped, 
     file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/neutral/RData/neutral_grouped.RData")


#3: fill in missing peaks----------
fpp <- FillChromPeaksParam() #default param very conservative, use them
neutral_grouped_filled <- fillChromPeaks(neutral_grouped)

save(neutral_grouped_filled, 
     file = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/3_FITDOG/3b_orbitrap_dec2020/analysis/neutral/RData/neutral_grouped_filled.RData")

neutral <- neutral_grouped_filled
rm(neutral_grouped_filled)

#4: save diffreport of xdata -----
xset <- as(neutral, "xcmsSet")
sampnames(xset) <- pData(neutral)$name
sampclass(xset) <- pData(neutral)$sample_type

#5. isotope picking and filtering ----
an <- xsAnnotate(xset)
an <- groupFWHM(an, 
                perfwhm = 0.6)
an <- findIsotopes(an, 
                   mzabs=0.01)
an <- groupCorr(an, 
                cor_eic_th=0.75)
an <- findAdducts(an, 
                  polarity="positive")
pl <- getPeaklist(an)

#filter by retention time 
#remove everything before 5 min and everything after 35 min
pl_rt <- pl %>% 
    filter(between(rt,
                   300,
                   2100))

#Filter by blank exclusion
pl_rt_be <-pl_rt[pl_rt$neg==0,]

#make rownames from rt and mz of features
rownames(pl_rt_be)<-paste(round(pl_rt_be$rt,1),
                          round(pl_rt_be$mz,3),
                          sep="_")

#filter for isotopes or adducts
pl_rt_be_isoadd <- pl_rt_be[pl_rt_be$isotopes!=""| 
                                pl_rt_be$adduct!="",]

pl_rt_isoadd <- pl_rt[pl_rt$isotopes!=""| 
                          pl_rt$adduct!="",]


#change NA to 0
pl_rt_be_isoadd[is.na(pl_rt_be_isoadd)] <- 0
pl_rt_isoadd[is.na(pl_rt_isoadd)] <- 0


#6: get features----
pl_rt_be_isoadd_neu <- pl_rt_be_isoadd[pl_rt_be_isoadd$neutral.fitdog >= 1 |
                                           pl_rt_be_isoadd$neutral.pos >= 1,]

pl_rt_isoadd_neu <- pl_rt_isoadd[pl_rt_isoadd$neutral.fitdog >= 1 |
                                     pl_rt_isoadd$neutral.pos >= 1,]

neutral_peaks <- pl_rt_isoadd_neu %>% 
    filter(solvent.blank..1 < 1e+05 & 
               solvent.blank..2 < 1e+05 & 
               blank.acetone.precipitation < 1e+05 & 
               blank.procainamide.reaction < 1e+05 &
               fitdog.neg < 1e+05)
neutral_peaks <- neutral_peaks %>% 
    filter(laminarin.fitdog..1 > 1e+05 &
               laminarin.fitdog..2 > 1e+05 |
               laminarin.standard.mix..1 > 1e+05  &
               laminarin.standard.mix..2 > 1e+05  |
               yeast.mannan.fitdog..1 > 1e+05  &
               yeast.mannan.fitdog..2 > 1e+05  |
               yeast.mannan...GH76..1 > 1e+05 &
               yeast.mannan...GH76..2 > 1e+05)

neutral_peaks <- cbind(rt_round = round_any(neutral_peaks$rt, 
                                            5),
                       neutral_peaks)

##collapse features with multiple isotopes
setDT(neutral_peaks)
#split out features without an isotope detected
neutral_peaks_noiso <- neutral_peaks[neutral_peaks$isotopes=="",]
neutral_peaks_iso <- neutral_peaks[!neutral_peaks$isotopes=="",]
#make column for the isotope group
neutral_peaks_iso$isotope_group <- neutral_peaks_iso$isotopes %>% 
    sub("\\[M.*", "", .)
#order isotopes within each group correctly
neutral_peaks_iso$isotope_number <- neutral_peaks_iso$isotopes %>% 
    sub(".*\\[M\\].*", "0", .) %>% 
    sub(".*\\[M\\+", "", .) %>% 
    sub("\\].*", "", .) %>% 
    as.numeric()
neutral_peaks_iso <- neutral_peaks_iso[order(isotope_group, 
                                             isotope_number),]
#get concatenated list of isotopes per group
iso_concat <- neutral_peaks_iso[, 
                                list(isotopes = paste(isotopes, 
                                                      collapse = ', ')), 
                                by = isotope_group]
#remove duplicates within each isotope group (will keep [M] isotope)
#because of ordering
neutral_peaks_iso <- unique(neutral_peaks_iso, 
                            by = "isotope_group")
#merge to get concatenated isotope lists
neutral_peaks_iso <- merge(neutral_peaks_iso,
                           iso_concat,
                           by = "isotope_group")
#clean up df
neutral_peaks_iso <- neutral_peaks_iso %>% 
    select(-c("isotope_group",
              "isotope_number",
              "isotopes.x"))
names(neutral_peaks_iso)[names(neutral_peaks_iso) == 'isotopes.y'] <- 'isotopes'

#replace features that don't contain [M] isotope with [M] isotope
temp <- neutral_peaks_iso %>% 
    filter(!grepl("\\[M\\]", isotopes))

temp.vec <- temp$isotopes %>% 
    sub("\\[M.*", "", .)

pl_rt_isoadd_neu$isotope_group <- pl_rt_isoadd_neu$isotopes %>% 
    sub("\\[M.*", "", .)
pl_rt_isoadd_neu$isotope_number <- pl_rt_isoadd_neu$isotopes %>% 
    sub(".*\\[M\\].*", "0", .) %>% 
    sub(".*\\[M\\+", "", .) %>% 
    sub("\\].*", "", .) %>% 
    as.numeric()

temp_2 <- pl_rt_isoadd_neu[pl_rt_isoadd_neu$isotope_group %in% temp.vec,] %>% 
    filter(isotope_number == 0)

setDT(temp_2)
temp_2 <- temp_2[order(isotope_group),]

temp$isotope_group <- temp$isotopes %>% 
    sub("\\[M.*", "", .)
setDT(temp)
temp <- temp[order(isotope_group),]
temp_3 <- temp[temp$isotope_group %in% temp_2$isotope_group,]

temp_2$isotopes <- paste(temp_2$isotopes,
                         temp_3$isotopes,
                         sep= ", ")

temp_2$isotope_group <- NULL
temp_2$isotope_number <- NULL

temp_2 <- cbind(rt_round = round_any(temp_2$rt, 
                                     5),
                temp_2)

neutral_peaks_iso <- neutral_peaks_iso %>% 
    filter(grepl("\\[M\\]", isotopes))

neutral_peaks_iso <- rbind(neutral_peaks_iso, 
                           temp_2)

rm(temp,
   temp.vec,
   temp_2,
   temp_3)

#merge features with and without isotopes
neutral_peaks <- rbind.fill(neutral_peaks_noiso,
                            neutral_peaks_iso)

setDT(neutral_peaks)
neutral_peaks <- neutral_peaks[order(rt_round,
                                     mz),]

write.csv(neutral_peaks,
          file = "./analysis/neutral/analysis_tables/neutral_peaklist_2.csv") 

##annotate features with predicted ions
#import mz values predicted using my script with the command:
#sugarMassesPredict-procainamide.py -dp 1 10 -p 0 -m none -n 1 -ld D -oh 0 -b 0 -i pos -s 250 2000 -o neutral-masses-predicted.txt -l procainamide
#manually added [M+H+K]2+, [M+2H]2+, [M+H+Na]2+ [M-CH2OH+H] ions 
mz.predicted.table <- fread("./analysis/neutral/neutral-masses-predicted.txt")

#remove "extra" columns
extraCol <- c('dp',
              'formula')

mz.predicted.table <- mz.predicted.table %>% 
    select(-all_of(extraCol))

#convert to long format
mz.predicted.table <- gather(mz.predicted.table, key = "ion", value = "mz", -name, -mass)

#make data.table
setDT(mz.predicted.table)
#create interval to overlap with (same width as for peak grouping)
mz.predicted.table$mz <- as.numeric(mz.predicted.table$mz)
mz.predicted.table$mzmin <- mz.predicted.table$mz-0.005
mz.predicted.table$mzmax <- mz.predicted.table$mz+0.005
#match using foverlaps from data.table (very fast)
setkey(mz.predicted.table, mzmin, mzmax)
neutral_peaks_matched <- foverlaps(neutral_peaks,
                                   mz.predicted.table)

neutral_peaks_matched <- neutral_peaks_matched %>% 
    replace_na(list("name"="unknown",
                    "mass"="",
                    "ion" = "", 
                    "mz" = "", 
                    "mzmin" = "",
                    "mzmax"= ""))

neutral_peaks_matched$name[neutral_peaks_matched$name==""] <- "unknown"

for (i in 1:nrow(neutral_peaks_matched)){
    if (neutral_peaks_matched$mz[i] == "") {
        neutral_peaks_matched$mz[i] <-  neutral_peaks_matched$i.mz[i]
    }
}


write.csv(neutral_peaks_matched,
          file = "./analysis/neutral/analysis_tables/neutral_peaklist_matched_2.csv") 

#format ion names for plot
neutral_peaks_matched$mz <- as.numeric(neutral_peaks_matched$mz)

neutral_peaks_matched$name_ion <- paste0(neutral_peaks_matched$name,
                                         ":",
                                         neutral_peaks_matched$ion,
                                         " mz=",
                                         round(neutral_peaks_matched$mz,
                                               3))
neutral_peaks_matched <- neutral_peaks_matched[order(name),] 

neutral_peaks_matched_all <- neutral_peaks_matched
neutral_peaks_final <- neutral_peaks_matched[neutral_peaks_matched$name!="unknown"]

#7:extract and format eic -----
##extract eic 
neutral_ions <- neutral_peaks_final %>% 
    distinct(name_ion, .keep_all = TRUE)

neutral_mz.found.vector <- neutral_ions$i.mz %>% 
    round(., 3) %>% 
    unique()
neutral_ions.found.vector  <- neutral_ions$name_ion

data_neutral <- filterFile(neutral,
                           file = which(grepl("neutral", 
                                              neutral$sample_type)))

neutral.names <- data_neutral$name
neutral.groups <- data_neutral$sample_type
neutral.substrate <- data_neutral$substrate
neutral.rep <- data_neutral$rep
neutral.rep[neutral.rep=="NA"] <- 1
neutral.rep <- as.numeric(neutral.rep)

#get chromatograms
neutral_chr_list <- list()
error = 0.001

for (i in 1:length(neutral_mz.found.vector)){
    mzr = c(neutral_mz.found.vector[i] - error,
            neutral_mz.found.vector[i] + error)
    neutral_chr_list[[i]] <- chromatogram(data_neutral, 
                                          mz = mzr)
}

save(neutral_chr_list,
     file = "./analysis/neutral/RData/neutral_chr_list.RData")

#extract intensity and rt values
neutral_chr_int_list <- list()
for (i in 1:length(neutral.names)){
    neutral_chr_int_list[[i]] <- lapply(neutral_chr_list, function(x) {
        x[[i]]@intensity
    }) 
}

neutral_chr_rt_list <- list()
for (i in 1:length(neutral.names)){
    neutral_chr_rt_list[[i]] <- lapply(neutral_chr_list, function(x) {
        x[[i]]@rtime
    }) 
}

#build data frame (long format)
neutral.df <- data.frame(ion = as.character(),
                         sample = as.character(),
                         group = as.character(),
                         substrate = as.character(),
                         rep = as.numeric(),
                         rt = as.numeric(),
                         intensity = as.numeric())

for (i in 1:length(neutral.names)){
    for (j in 1:length(neutral_mz.found.vector)){
        rt = neutral_chr_rt_list[[i]][[j]]
        intensity = neutral_chr_int_list[[i]][[j]]
        sample = rep(neutral.names[i], length(rt))
        group = rep(neutral.groups[i], length(rt))
        substrate = rep(neutral.substrate[i], length(rt))
        ion = rep(neutral_ions.found.vector[j], length(rt))
        rep = rep(neutral.rep[i], length(rt))
        temp <- data.frame(ion = ion,
                           sample = sample,
                           group = group,
                           substrate = substrate,
                           rep = rep,
                           rt = rt,
                           intensity = intensity)
        neutral.df <- rbind(neutral.df,
                            temp)
    }
}

neutral.df[is.na(neutral.df)] <- 0

#set variables as factors
neutral.df$ion <-  neutral.df$ion %>% 
    sub("hex-1-procA", "monosaccharide", .) %>% 
    sub("hex-2-procA", "disaccharide", .) %>% 
    sub("hex-3-procA", "trisaccharide", .) %>%
    sub("hex-4-procA", "tetrasaccharide", .) %>%
    sub("hex-5-procA", "pentasaccharide", .) %>%
    sub("hex-6-procA", "hexaccharide", .)

neutral_ions.found.vector <-  neutral_ions.found.vector %>% 
    sub("hex-1-procA", "monosaccharide", .) %>% 
    sub("hex-2-procA", "disaccharide", .) %>% 
    sub("hex-3-procA", "trisaccharide", .) %>%
    sub("hex-4-procA", "tetrasaccharide", .) %>%
    sub("hex-5-procA", "pentasaccharide", .) %>%
    sub("hex-6-procA", "hexaccharide", .)

neutral.df$ion <- factor(neutral.df$ion,
                         levels = neutral_ions.found.vector)


neutral.df$group_fmt <- ""

neutral.df$group_fmt[neutral.df$substrate == "laminarin"] <- neutral.df$group[neutral.df$substrate == "laminarin"] %>% 
    sub("neutral.fitdog", "FITDOG", .) %>% 
    sub("neutral.pos", "standard mix", .)
neutral.df$group_fmt[neutral.df$substrate == "yeast mannan"] <- neutral.df$group[neutral.df$substrate == "yeast mannan"] %>% 
    sub("neutral.fitdog", "FITDOG", .) %>% 
    sub("neutral.pos", "GH76 enzyme digest", .) %>% 
    sub("neutral.neg", "undigested", .)
neutral.df$group_fmt <- factor(neutral.df$group_fmt,
                               levels = unique(neutral.df$group_fmt))

neutral.df$rt_min <- neutral.df$rt / 60

neutral.df$rep <- factor(neutral.df$rep,
                         levels = c(1,2))

neutral.df$substrate <- factor(neutral.df$substrate,
                               levels = unique(neutral.df$substrate))


#reformat ion names and only keep [M+H]+ ions
neutral_ions$new_name_ion <- neutral_ions$name_ion %>% 
    sub("hex-1-procA:\\[M-CH2O\\+H\\]\\+",
        "1Pentose:[M+H]+", .) %>% 
    sub("hex-1-procA:\\[M-H\\]\\+", 
        "1Hexenose:[M+H]+", .) %>% 
    sub("hex-1-procA",
        "1Hexose", .) %>% 
    sub("hex-2-procA:\\[M-CH2O\\+H\\]\\+",
        "1Hexose1Pentose:[M+H]+", .) %>% 
    sub("hex-2-procA:\\[M-H\\]\\+",
        "1Hexose1Hexenose:[M+H]+", .) %>% 
    sub("hex-2-procA",
        "2Hexose", .) %>% 
    sub("hex-3-procA:\\[M-H\\]\\+",
        "2Hexose1Hexenose:[M+H]+", .) %>% 
    sub("hex-3-procA:\\[M-CH2O\\+H\\]\\+",
        "2Hexose1Pentose:[M+H]+", .) %>% 
    sub("hex-3-procA",
        "3Hexose", .) %>% 
    sub("hex-4-procA:\\[M-H\\]\\+",
        "3Hexose1Hexenose:[M+H]+", .) %>% 
    sub("hex-4-procA",
        "4Hexose", .) %>% 
    sub("hex-5-procA:\\[M-H\\]\\+",
        "4Hexose1Hexenose:[M+H]+", .) %>% 
    sub("hex-5-procA",
        "5Hexose", .) %>% 
    sub("hex-6-procA",
        "6Hexose", .)

neutral_ions$new_name <- neutral_ions$new_name_ion %>% 
    sub(":.*", "", .)

neutral_ions$new_ion <- neutral_ions$new_name_ion %>% 
    gsub(".*:|\\smz.*", "", .)

neutral_ions_new <- neutral_ions[neutral_ions$new_ion=="[M+H]+",]

neutral.df$new_ion <-  neutral.df$ion %>% 
    sub("monosaccharide:\\[M-CH2O\\+H\\]\\+",
        "1Pentose:[M+H]+", .) %>% 
    sub("monosaccharide:\\[M-H\\]\\+", 
        "1Hexenose:[M+H]+", .) %>% 
    sub("monosaccharide",
        "1Hexose", .) %>% 
    sub("disaccharide:\\[M-CH2O\\+H\\]\\+",
        "1Hexose1Pentose:[M+H]+", .) %>% 
    sub("disaccharide:\\[M-H\\]\\+",
        "1Hexose1Hexenose:[M+H]+", .) %>% 
    sub("disaccharide",
        "2Hexose", .) %>% 
    sub("trisaccharide:\\[M-H\\]\\+",
        "2Hexose1Hexenose:[M+H]+", .) %>% 
    sub("trisaccharide:\\[M-CH2O\\+H\\]\\+",
        "2Hexose1Pentose:[M+H]+", .) %>% 
    sub("trisaccharide",
        "3Hexose", .) %>% 
    sub("tetrasaccharide:\\[M-H\\]\\+",
        "3Hexose1Hexenose:[M+H]+", .) %>% 
    sub("tetrasaccharide",
        "4Hexose", .) %>% 
    sub("pentasaccharide:\\[M-H\\]\\+",
        "4Hexose1Hexenose:[M+H]+", .) %>% 
    sub("pentasaccharide",
        "5Hexose", .) %>% 
    sub("hexaccharide",
        "6Hexose", .)

neutral.df_old <- neutral.df

neutral.df <- neutral.df[neutral.df$new_ion %in% 
                             neutral_ions_new$new_name_ion,]

#rename mz to m/z
neutral.df$new_ion <- neutral.df$new_ion %>% 
    sub("mz", "m/z", .)

neutral.df$new_ion <- factor(neutral.df$new_ion,
                             levels = c(unique(neutral.df$new_ion)))

#set variable

neutral.df$var <- "EIC"

#reset group variable
neutral.df$group_fmt_new <- paste0(neutral.df$group_fmt,
                                   ": replicate ",
                                   neutral.df$rep)

neutral.df$group_fmt_new <- neutral.df$group_fmt_new %>% 
    sub("standard mix.*", "standard mix", .) %>% 
    sub("GH76 enzyme digest.*", "GH76 enzyme digest", .)

neutral.df$group_ion <- paste0(neutral.df$group_fmt_new,
                               "_",
                               neutral.df$new_ion)

neutral.df_old2 <- neutral.df

extraCol <- c("ion",
              "group",
              "rt",
              "group_fmt")

neutral.df <- neutral.df %>% 
    select(-extraCol)

#identify ions with zero intensity in samples
neutral_ionIntensitySums <- neutral.df %>% 
    group_by(substrate, group_ion) %>% 
    summarise(sum = sum(intensity), .groups="keep")

neutral_zeroIntensityIons <- neutral.df %>% 
    group_by(substrate, group_ion) %>% 
    summarise(sum = sum(intensity), .groups="keep") %>% 
    filter(sum == 0)

lam_zeroIntensityIons <- neutral_zeroIntensityIons$group_ion[
    neutral_zeroIntensityIons$substrate == "laminarin"]

ym_zeroIntensityIons <- neutral_zeroIntensityIons$group_ion[
    neutral_zeroIntensityIons$substrate == "yeast mannan"]


#split by substrate and remove zero intensity ions
lam_neutral.df <- neutral.df[neutral.df$substrate == "laminarin",]
lam_neutral.df <- subset(lam_neutral.df, 
                         !(group_ion %in% lam_zeroIntensityIons))

ym_neutral.df <- neutral.df[neutral.df$substrate == "yeast mannan",]
ym_neutral.df <- subset(ym_neutral.df, 
                        !(group_ion %in% ym_zeroIntensityIons))

#order group and ion variable for plotting
ym_group_ion <- ym_neutral.df$group_ion %>% unique %>% sort

ym_group_ion <- c(ym_group_ion[25:38],
                  ym_group_ion[1:24])


ym_neutral.df$group_ion <- factor(ym_neutral.df$group_ion,
                                  levels = ym_group_ion)




#8: import FLD and format-----------
#read in text files
fld_fp <- dir(path = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_mannan/triple-quad/Data_final/FLD-files", 
              all.files = FALSE, 
              full.names = TRUE)
fld_fp <- fld_fp[grep("ym|lam_" ,
                      fld_fp)]
fld_list <- lapply(fld_fp, 
                   fread)

#make data frame
fld.groups <- fld_fp %>% 
    sub(".*fitdog.*", "FITDOG", .) %>% 
    sub(".*nodigest.*", "undigested", .) %>% 
    sub(".*lymGH76.*", "GH76 enzyme digest", .) %>% 
    sub(".*omix.*", "standard mix", .)

fld.sample <- basename(fld_fp)%>% 
    sub("2020.*HILIC_", "", .) %>% 
    sub("-FLD.txt", "", .)

fld.substrate <- basename(fld_fp)%>% 
    sub(".*lam.*", "laminarin", .) %>% 
    sub(".*ym.*", "yeast mannan", .)

fld.df <- data.frame(sample = as.character(),
                     group_fmt = as.character(),
                     substrate = as.character(),
                     rt = as.numeric(),
                     intensity = as.numeric())

for (i in 1:length(fld.sample)){
    rt = fld_list[[i]]$V1
    intensity = fld_list[[i]]$V2
    sample = rep(fld.sample[i], length(rt))
    group_fmt = rep(fld.groups[i], length(rt))
    substrate = rep(fld.substrate[i], length(rt))
    temp <- data.frame(sample = sample,
                       group_fmt = group_fmt,
                       substrate = substrate,
                       rt = rt,
                       intensity = intensity)
    fld.df <- rbind(fld.df,
                    temp)
}


fld.df$group_fmt <- factor(fld.df$group_fmt,
                           levels = unique(fld.df$group_fmt))
fld.df$substrate <- factor(fld.df$substrate,
                           levels = unique(fld.df$substrate))

#transform FLD on x axis - done manually to fit peaks
fld.df$rt_trans <- fld.df$rt - 0.75

#get FLD minimum and offset each sample so that minimum is zero
#rt range at which minimum is found based on initial plots
fld.df.notzeroed <- fld.df #keep to be safe

for (i in 1:length(fld.sample)){
    x <- fld.df.notzeroed[fld.df.notzeroed$sample==fld.sample[i],]
    int.min = x %>% 
        filter(between(rt, 5, 25)) %>% 
        select(intensity) %>% 
        min()
    x$intensity <- x$intensity + abs(int.min)
    fld.df$intensity[fld.df$sample==fld.sample[i]] <- x$intensity
}


#13: plot FLR ----
fld.df$group_fmt <- fld.df$group_fmt %>% 
    sub("GH76.*|standard.*", "positive control", .)


fld.df$group_fmt <- factor(fld.df$group_fmt,
                           levels =c("positive control",
                                     "FITDOG",
                                     "undigested"))


pal_fld <-  c("#182031",
              "#B02158",
              "#ACA4E2")

#make breaks
x <- seq(5, 25, 1)
major_breaks_zoom <- vector(length = length(x),
                            mode = "character")

for (i in 1:length(major_breaks_zoom)){
    if (x[i] %% 2 != 0){
        major_breaks_zoom[i] <- as.character(x[i])
    } else if (x[i] %% 2 == 0) {
        major_breaks_zoom[i] <- ""
    }
}


#subset FLD
fld.df.zoom <- fld.df %>% 
    filter(between(rt_trans, 5, 25))

#plot
tiff("./analysis/neutral/analysis_plots/test1.tiff", 
     res = 300, 
     height = 6, 
     width = 12, 
     units = "in")

svg("./analysis/neutral/analysis_plots/fitdog_flr_summary_v2.svg", 
    height = 6, 
    width = 12)

ggplot() +
    geom_line(mapping = aes(rt_trans,
                            intensity,
                            colour = group_fmt,
                            group = group_fmt),
              data = fld.df.zoom,
              lwd = 1.2) +
    scale_colour_manual(name = "",
                        values = pal_fld) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 12),
          legend.position = "none") +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    facet_grid(rows = vars(substrate),
               #scales = "free_y"
    ) +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x = element_blank(),
          axis.text = element_text(size = 12)) + 
    scale_x_continuous(breaks = seq(5, 25, 1),
                       labels = major_breaks_zoom,
                       limits = c(5, 25),
                       expand = c(0,0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                       labels = scales::scientific)

dev.off()


#14: make and format tables of ions ----

neutral_peaks_final$new_name_ion <- neutral_peaks_final$name_ion %>% 
    sub("hex-1-procA:\\[M-CH2O\\+H\\]\\+",
        "1Pentose:[M+H]+", .) %>% 
    sub("hex-1-procA:\\[M-H\\]\\+", 
        "1Hexenose:[M+H]+", .) %>% 
    sub("hex-1-procA",
        "1Hexose", .) %>% 
    sub("hex-2-procA:\\[M-CH2O\\+H\\]\\+",
        "1Hexose1Pentose:[M+H]+", .) %>% 
    sub("hex-2-procA:\\[M-H\\]\\+",
        "1Hexose1Hexenose:[M+H]+", .) %>% 
    sub("hex-2-procA",
        "2Hexose", .) %>% 
    sub("hex-3-procA:\\[M-H\\]\\+",
        "2Hexose1Hexenose:[M+H]+", .) %>% 
    sub("hex-3-procA:\\[M-CH2O\\+H\\]\\+",
        "2Hexose1Pentose:[M+H]+", .) %>% 
    sub("hex-3-procA",
        "3Hexose", .) %>% 
    sub("hex-4-procA:\\[M-H\\]\\+",
        "3Hexose1Hexenose:[M+H]+", .) %>% 
    sub("hex-4-procA",
        "4Hexose", .) %>% 
    sub("hex-5-procA:\\[M-H\\]\\+",
        "4Hexose1Hexenose:[M+H]+", .) %>% 
    sub("hex-5-procA",
        "5Hexose", .) %>% 
    sub("hex-6-procA",
        "6Hexose", .)

neutral_peaks_final$name <- sub(":.*", 
                                "", 
                                neutral_peaks_final$new_name_ion)

neutral_peaks_final$ion <- gsub(".*:|\\sm.*", 
                                "", 
                                neutral_peaks_final$new_name_ion)

neutral_peaks_final$calculated_mz <- round(neutral_peaks_final$mz, 3)
neutral_peaks_final$observed_mz <- round(neutral_peaks_final$i.mz, 3)
neutral_peaks_final$rt_min<- round(neutral_peaks_final$rt/60, 1)

res <- neutral_peaks_final

keepCol <- c("name", "ion", "calculated_mz", "observed_mz", "rt_min", "isotopes", 
             names(neutral_peaks_final)[
                       grep("laminarin|yeast.mannan",
                            names(neutral_peaks_final))])

res <- res %>% 
    select(keepCol)

sampleCol <- names(neutral_peaks_final)[
    grep("laminarin|yeast.mannan",
         names(neutral_peaks_final))]

setDF(res)

res[sampleCol][res[sampleCol] > 0] <- 1

fwrite(res,
       "./analysis/neutral/analysis_tables/neutralFITDOG_ionTable.txt",
       quote = FALSE,
       sep = ";")

#make feature names
res2 <- res
res2$calculated_mz <- NULL
res2$observed_mz <- NULL
res2$isotopes <- NULL
res2$feature <- paste0(res2$name,
                       " (",
                       res2$rt_min,
                       " min)")
#make long format
res2_wide <- res2
res2 <- gather(res2_wide,
               key = "sample",
               value = "count",
               -name,
               -ion,
               -rt_min,
               -feature)

#gather ions together and make it so that count per sample is maximum one
res2$count <- as.numeric(res2$count)
res2 <- res2 %>% 
    group_by(feature, sample, rt_min) %>% 
    dplyr::summarise(count = sum(count))
res2$count[res2$count > 1] <- 1

#get counts per sample type
res2$sample_type <- res2$sample %>% 
    sub("\\.\\.[12]$", "", .)
res2 <- res2 %>% 
    group_by(feature, 
             sample_type,
             rt_min) %>% 
    dplyr::summarise(count_sum=sum(count))

#15: plotting representation of results summary ----

#order by retention time
res2 <- res2[order(res2$rt_min),]

#set variables to factors
res2$feature <- factor(res2$feature,
                      levels= unique(res2$feature))

res2$type <- res2$sample_type %>% 
    sub(".*fitdog*", "FITDOG", .) %>% 
    sub(".*standard.*|.*GH76.*", "positive control", .)
res2$type <- factor(res2$type, 
                   levels = c("FITDOG",
                              "positive control"))

res2$substrate <- res2$sample_type %>% 
    sub("laminarin.*", "laminarin", .) %>% 
    sub("yeast.*", "yeast mannan", .)

res2$substrate <- factor(res2$substrate, 
                        levels = c("laminarin",
                                   "yeast mannan"))

#replace 0 with NA
res2$count_sum[res2$count_sum== 0] <- NA


#res2$count_sum_char <- as.character(res2$count_sum)

tiff("./analysis/neutral/test3.tiff",
     units = "in",
     res = 300,
     width = 8,
     height = 6)

svg("./analysis/neutral/peakcount_summary_v1.svg",
    width = 8,
    height = 6)


ggplot(res2,
       aes(x = substrate,
           y = feature,
           colour = type,
           size = count_sum,
           group = type)) +
    scale_size_continuous(breaks = c(1,2)) +
    scale_colour_manual(values = c("#182031",
                                   "#B02158"),
                        name = "Sample type",
                        labels = c("FITDOG", "positive control")) +
    geom_point(stat = "identity",
               position = position_dodge(width = 0.3)) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          axis.text = element_text(size = 12)) +
    labs(x= "Substrate",
         y = "Annotated features")
dev.off()    



#FLR/EIC GH76 ONLY -------

full_gh76 <- full_fld_eic.df[full_fld_eic.df$group_fmt == "GH76 enzyme digest",]
zoom_gh76 <- fld_eic.df[fld_eic.df$group_fmt == "GH76 enzyme digest",]

gh76_ions <- c("FLR",
               "monosaccharide:[M+H]+ mz=400.245",
               "disaccharide:[M+H]+ mz=562.298",
               "trisaccharide:[M+H]+ mz=724.35",
               "tetrasaccharide:[M+H]+ mz=886.403")


full_gh76 <- full_gh76[full_gh76$ion %in% gh76_ions, ]
zoom_gh76 <- zoom_gh76[zoom_gh76$ion %in% gh76_ions, ]

full_gh76$var <- full_gh76$var %>% 
    sub("EIC", "all extracted ion chromatograms", .) %>% 
    sub("FLR", "fluorescent signal", .)

full_gh76$var <- factor(full_gh76$var,
                        levels = c("fluorescent signal",
                                   "all extracted ion chromatograms"))

zoom_gh76$ion <- zoom_gh76$ion %>% 
    sub("FLR", "fluorescent signal", .)

zoom_gh76$ion <- factor(zoom_gh76$ion,
                        levels = unique(zoom_gh76$ion))

tiff("./analysis/neutral/EIC_FLR/gh76_full_v2.tiff", 
     res = 300, 
     height = 9, 
     width = 12, 
     units = "in")

gh76_p1 <- ggplot() +
    geom_line(mapping = aes(rt_min,
                            intensity,
                            color = ion),
              data = full_gh76,
              size = 1) +
    scale_colour_manual(name = "",
                        values = pal_ionshort) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 1),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom"
    ) +
    facet_grid(rows = vars(var),
               cols = vars(group_fmt),
               scales = "free_y"
    ) +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    scale_x_continuous(breaks = seq(0, 60, 1),
                       labels = major_breaks_full,
                       limits = c(0,60),
                       expand = c(0,0)) + 
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                       labels = scales::scientific ) +
    theme(strip.background = element_blank(),
          strip.text.y = element_text(angle = 360, hjust = 0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_vline(xintercept = 5, 
               linetype="longdash", 
               color = "#8c001aff", 
               size=1) + 
    geom_vline(xintercept = 25, 
               linetype="longdash", 
               color = "#8c001aff", 
               size=1)

dev.off()

tiff("./analysis/neutral/EIC_FLR/gh76_v2.tiff", 
     res = 300, 
     height = 9, 
     width = 12, 
     units = "in")

gh76_p2 <- ggplot(zoom_gh76,
                  aes(rt_min,
                      intensity,
                      color = ion)) +
    geom_line(size = 1) +
    scale_colour_manual(name = "",
                        values = pal_ionshort) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 1),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none"
    ) +
    facet_grid(rows = vars(ion),
               cols = vars(group_fmt),
               #scales = "free_y"
    ) +
    labs(x= "Retention time (min)") +
    scale_x_continuous(breaks = seq(5, 25, 1),
                       labels = major_breaks_zoom,
                       limits = c(5, 25),
                       expand = expansion(mult = c(0.02, 0.02))) + 
    theme(strip.background = element_blank(),
          strip.text.y = element_text(angle = 360, hjust = 0),
          strip.text.x = element_blank()
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(name = "Intensity (a.u.)",
                       labels = scales::scientific,
                       breaks = waiver(),
                       n.breaks = 3) +
    geom_vline(xintercept = 5, 
               linetype="longdash", 
               color = "#8c001aff", 
               size=1) + 
    geom_vline(xintercept = 25, 
               linetype="longdash", 
               color = "#8c001aff", 
               size=1)

dev.off()

gh1 <- ggplotGrob(gh76_p1)
gh2 <- ggplotGrob(gh76_p2)
gh <- rbind(gh1, gh2, size = "first")
gh$widths <- unit.pmax(gh1$widths, gh2$widths)


tiff("./analysis/neutral/EIC_FLR/gh76_v1.tiff", 
     res = 300, 
     height = 9, 
     width = 12, 
     units = "in")
grid.newpage()
grid.draw(gh)
dev.off()

#FLR/EIC FUCOIDAN ONLY -------
#1:extract and format eic -----
##extract eic 
fucoidan_mz.found.vector <- neutral_mz.found.vector[c(1:2,4)]
fucoidan_ions.found.vector <- neutral_ions.found.vector[c(1:2,4)]


data_fucoidan <- filterFile(neutral,
                            file = which(grepl("fucoidan", 
                                               neutral$name)))

fucoidan.names <- data_fucoidan$name
fucoidan.groups <- c("FITDOG", "undigested")


#get chromatograms
fucoidan_chr_list <- list()
error = 0.001

for (i in 1:length(fucoidan_mz.found.vector)){
    mzr = c(fucoidan_mz.found.vector[i] - error,
            fucoidan_mz.found.vector[i] + error)
    fucoidan_chr_list[[i]] <- chromatogram(data_fucoidan, 
                                           mz = mzr)
}

save(fucoidan_chr_list,
     file = "./analysis/fucoidan/RData/fucoidan_chr_list.RData")

#extract intensity and rt values
fucoidan_chr_int_list <- list()
for (i in 1:length(fucoidan.names)){
    fucoidan_chr_int_list[[i]] <- lapply(fucoidan_chr_list, function(x) {
        x[[i]]@intensity
    }) 
}

fucoidan_chr_rt_list <- list()
for (i in 1:length(fucoidan.names)){
    fucoidan_chr_rt_list[[i]] <- lapply(fucoidan_chr_list, function(x) {
        x[[i]]@rtime
    }) 
}

#build data frame (long format)
fucoidan.df <- data.frame(ion = as.character(),
                          sample = as.character(),
                          group = as.character(),
                          rt = as.numeric(),
                          intensity = as.numeric())

for (i in 1:length(fucoidan.names)){
    for (j in 1:length(fucoidan_mz.found.vector)){
        rt = fucoidan_chr_rt_list[[i]][[j]]
        intensity = fucoidan_chr_int_list[[i]][[j]]
        sample = rep(fucoidan.names[i], length(rt))
        group = rep(fucoidan.groups[i], length(rt))
        ion = rep(fucoidan_ions.found.vector[j], length(rt))
        temp <- data.frame(ion = ion,
                           sample = sample,
                           group = group,
                           rt = rt,
                           intensity = intensity)
        fucoidan.df <- rbind(fucoidan.df,
                             temp)
    }
}

fucoidan.df[is.na(fucoidan.df)] <- 0

#set variables as factors
fucoidan.df$ion <- factor(fucoidan.df$ion,
                          levels = fucoidan_ions.found.vector)

fucoidan.df$group <- factor(fucoidan.df$group,
                            levels = unique(fucoidan.df$group))

fucoidan.df$rt_min <- fucoidan.df$rt / 60



#2: import FLD and format-----------
#read in text files
fld_fp <- dir(path = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_mannan/triple-quad/Data_final/FLD-files", 
              all.files = FALSE, 
              full.names = TRUE)
fld_fp <- fld_fp[grep("fucoidan" ,
                      fld_fp)]
fld_list <- lapply(fld_fp, 
                   fread)

#make data frame
fld.groups <- fld_fp %>% 
    sub(".*fitdog.*", "FITDOG", .) %>% 
    sub(".*nodigest.*", "undigested", .)

fld.sample <- basename(fld_fp)%>% 
    sub("2020.*HILIC_", "", .) %>% 
    sub("-FLD.txt|.txt", "", .)

fuc_fld.df <- data.frame(sample = as.character(),
                         group = as.character(),
                         rt = as.numeric(),
                         intensity = as.numeric())

for (i in 1:length(fld.sample)){
    rt = fld_list[[i]]$V1
    intensity = fld_list[[i]]$V2
    sample = rep(fld.sample[i], length(rt))
    group = rep(fld.groups[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rt = rt,
                       intensity = intensity)
    fuc_fld.df <- rbind(fuc_fld.df,
                        temp)
}


fuc_fld.df$ion <- "FLR"
fuc_fld.df$group <- factor(fuc_fld.df$group,
                           levels = unique(fuc_fld.df$group))

#transform FLD on x axis - done manually to fit peaks
fuc_fld.df$rt_trans <- fuc_fld.df$rt - 0.75

#get FLD minimum and offset each sample so that minimum is zero
#rt range at which minimum is found based on initial plots
fuc_fld.df.notzeroed <- fuc_fld.df #keep to be safe

for (i in 1:length(fld.sample)){
    x <- fuc_fld.df.notzeroed[fuc_fld.df.notzeroed$sample==fld.sample[i],]
    int.min = x %>% 
        filter(between(rt, 5, 25)) %>% 
        select(intensity) %>% 
        min()
    x$intensity <- x$intensity + abs(int.min)
    fuc_fld.df$intensity[fuc_fld.df$sample==fld.sample[i]] <- x$intensity
}


#3: plot FLD: full chromatography, one plot per group ----
#make breaks
x <- seq(0, 60, 1)
major_breaks_full <- vector(length = length(x),
                            mode = "character")
for (i in 1:length(major_breaks_full)){
    if (x[i] %% 5 == 0){
        major_breaks_full[i] <- as.character(x[i])
    } else if (x[i] %% 5 != 0) {
        major_breaks_full[i] <- ""
    }
}


#plot
tiff("./analysis/fucoidan/analysis_plots/fucoidan_fld_full_v1.tiff", 
     res = 300, 
     height = 4.5, 
     width = 12, 
     units = "in")
svg("./analysis/fucoidan/analysis_plots/fucoidan_fld_full_v1.svg", 
    height = 4.5, 
    width = 12)

ggplot() +
    geom_line(mapping = aes(rt_trans,
                            intensity),
              data = fuc_fld.df,
              colour = "#8E8F91") +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 1),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    facet_grid(
        #rows = vars(ion),
        cols = vars(group),
        #scales = "free_y"
    ) +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)",
         title = "fucoidan LC-FLR") +
    scale_x_continuous(breaks = seq(0, 60, 1),
                       labels = major_breaks_full,
                       limits = c(0,60),
                       expand = c(0,0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_vline(xintercept = 5, 
               linetype="longdash", 
               color = "#8c001aff", 
               size=1) + 
    geom_vline(xintercept = 25, 
               linetype="longdash", 
               color = "#8c001aff", 
               size=1)

dev.off()

#4: plot EIC: full chromatography, one plot per group----

pal_fuc <- c("#E495A5",
             "#ABB065",
             "#ACA4E2")

tiff("./analysis/fucoidan/analysis_plots/fucoidan_eic_full_v1_legend.tiff", 
     res = 300, 
     height = 4.5, 
     width = 12, 
     units = "in")

svg("./analysis/fucoidan/analysis_plots/fucoidan_eic_full_v1_nolegend.svg", 
    height = 4.5, 
    width = 12)


ggplot() +
    geom_line(mapping = aes(rt_min,
                            intensity,
                            color = ion),
              data = fucoidan.df,
              size = 1) +
    scale_colour_manual(name = "Ions in EIC",
                        values = pal_fuc) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 1),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    facet_grid(#rows = vars(ion),
        cols = vars(group),
        #scales = "free_y"
    ) +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)",
         title = "fucoidan EIC") +
    scale_x_continuous(breaks = seq(0, 60, 1),
                       labels = major_breaks_full,
                       limits = c(0,60),
                       expand = c(0,0)) + 
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank(),
          legend.position = "none"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_vline(xintercept = 5, 
               linetype="longdash", 
               color = "#8c001aff", 
               size=1) + 
    geom_vline(xintercept = 25, 
               linetype="longdash", 
               color = "#8c001aff", 
               size=1)

dev.off()

#5: plot EIC and FLD: 5-25 min, facet by ion and group----
#make breaks
x <- seq(5, 25, 1)
major_breaks_zoom <- vector(length = length(x),
                            mode = "character")

for (i in 1:length(major_breaks_zoom)){
    if (x[i] %% 2 != 0){
        major_breaks_zoom[i] <- as.character(x[i])
    } else if (x[i] %% 2 == 0) {
        major_breaks_zoom[i] <- ""
    }
}

#subset FLD
fuc_fld.df.zoom <- fuc_fld.df %>% 
    filter(between(rt_trans, 5, 25))

#combine df
fuc_fld.df.zoom$rt_min <- fuc_fld.df.zoom$rt_trans

fuc_fld_eic.df <- rbind.fill(fuc_fld.df.zoom,
                             fucoidan.df)

fuc_fld_eic.df$ion <- factor(fuc_fld_eic.df$ion,
                             levels = c("FLR", fucoidan_ions.found.vector))


#make palette
pal_fuc_comb <- c("#8E8F91",
                  pal_fuc)

tiff("./analysis/fucoidan/analysis_plots/fucoidan_fld_eic_zoom_v2_nolegend.tiff", 
     res = 300, 
     height = 9, 
     width = 12, 
     units = "in")

svg("./analysis/fucoidan/analysis_plots/fucoidan_fld_eic_zoom_v2_nolegend.svg", 
    height = 9, 
    width = 12)

ggplot(fuc_fld_eic.df,
       aes(rt_min,
           intensity,
           color = ion)) +
    geom_line(size = 1) +
    scale_colour_manual(name = "Ions in EIC",
                        values = pal_fuc_comb) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 1),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none"
    ) +
    facet_grid(rows = vars(ion),
               cols = vars(group),
               scales = "free_y"
    ) +
    labs(x= "Retention time (min)",
         title = "fucoidan EIC and FLR") +
    scale_x_continuous(breaks = seq(5, 25, 1),
                       labels = major_breaks_zoom,
                       limits = c(5, 25),
                       expand = expansion(mult = c(0.02, 0.02))) + 
    theme(strip.background = element_blank(),
          strip.text.y = element_text(angle = 360, hjust = 0)
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(name = "Intensity (a.u.)",
                       labels = scales::scientific,
                       breaks = waiver(),
                       n.breaks = 3) +
    geom_vline(xintercept = 5, 
               linetype="longdash", 
               color = "#8c001aff", 
               size=1) + 
    geom_vline(xintercept = 25, 
               linetype="longdash", 
               color = "#8c001aff", 
               size=1)

dev.off()


#FLR/EIC FUCOIDAN and LAMBDA CARRAGEENAN -------
#1: import FLD and format-----------
#read in text files
fld_fp <- dir(path = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_mannan/triple-quad/Data_final/FLD-files", 
              all.files = FALSE, 
              full.names = TRUE)
fld_fp <- fld_fp[grep("fucoidan|lambda-carra" ,
                      fld_fp)]
fld_list <- lapply(fld_fp, 
                   fread)

#make data frame
fld.groups <- fld_fp %>% 
    sub(".*fitdog.*", "FITDOG", .) %>% 
    sub(".*nodigest.*|.*heat-inact.*", "undigested", .) %>% 
    sub(".*gh82.*", "GH82 enzyme digest", .) 

fld.substrate <- fld_fp %>% 
    sub(".*fucoidan.*", "fucoidan", .) %>% 
    sub(".*lambda-carra.*", "lambda carrageenan", .)

fld.sample <- basename(fld_fp) %>% 
    sub("2020-12-15_HILIC_", "", .) %>% 
    sub(".txt", "", .)

anionic_fld.df <- data.frame(sample = as.character(),
                             group = as.character(),
                             substrate = as.character(),
                             rt = as.numeric(),
                             intensity = as.numeric())

for (i in 1:length(fld.sample)){
    rt = fld_list[[i]]$V1
    intensity = fld_list[[i]]$V2
    sample = rep(fld.sample[i], length(rt))
    group = rep(fld.groups[i], length(rt))
    substrate = rep(fld.substrate[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       substrate = substrate,
                       rt = rt,
                       intensity = intensity)
    anionic_fld.df <- rbind(anionic_fld.df,
                            temp)
}


#get FLD minimum and offset each sample so that minimum is zero
#rt range at which minimum is found based on initial plots
anionic_fld.df.notzeroed <- anionic_fld.df #keep to be safe

for (i in 1:length(fld.sample)){
    x <- anionic_fld.df.notzeroed[anionic_fld.df.notzeroed$sample==fld.sample[i],]
    int.min = x %>% 
        filter(between(rt, 5, 25)) %>% 
        select(intensity) %>% 
        min()
    x$intensity <- x$intensity + abs(int.min)
    anionic_fld.df$intensity[anionic_fld.df$sample==fld.sample[i]] <- x$intensity
}


#2: plot FLD ----

anionic_fld.df$group_fmt <- anionic_fld.df$group %>% 
    sub("GH82.*", "positive control", .)


anionic_fld.df$group_fmt <- factor(anionic_fld.df$group_fmt,
                                   levels =c("positive control",
                                             "FITDOG",
                                             "undigested"))

anionic_fld.df_full <- anionic_fld.df 
anionic_fld.df <- anionic_fld.df_full %>% 
    filter(between(rt, 5, 40))


pal_fld <-  c("#182031",
              "#B02158",
              "#ACA4E2")

names(pal_fld) <- levels(anionic_fld.df$group_fmt)

#make breaks
x <- seq(5, 40, 1)
major_breaks_full <- vector(length = length(x),
                            mode = "character")
for (i in 1:length(major_breaks_full)){
    if (x[i] %% 5 == 0){
        major_breaks_full[i] <- as.character(x[i])
    } else if (x[i] %% 5 != 0) {
        major_breaks_full[i] <- ""
    }
}

svg("./analysis/neutral/analysis_plots/anionic_fitdog_flr_v2.svg", 
    height = 4.5, 
    width = 12)

ggplot() +
    geom_line(mapping = aes(rt,
                            intensity,
                            colour = group_fmt,
                            group = sample),
              data = anionic_fld.df,
              lwd = 1.2) +
    facet_grid(vars(substrate)) + 
    scale_colour_manual(name = "",
                        values = pal_fld) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 12),
          legend.position = "none",
          axis.text = element_text(size = 12),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    scale_x_continuous(breaks = seq(5, 40, 1),
                       labels = major_breaks_full,
                       limits = c(5, 40),
                       expand = c(0,0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                       labels = scales::scientific)

dev.off()


#EXTRACT AND PLOT MS2 -----
#1. load files with MS2----
neutral_fp <- dir(path = "./neutral_fitdog/mzML-files", 
                  all.files = FALSE, 
                  full.names = TRUE)
neutral_fp <- neutral_fp[grep("20201222", neutral_fp)]


pd_ms2 <- data.frame(name = basename(neutral_fp),
                 sample_type = neutral_fp %>% 
                     sub(".*omix.*|.*gh76.*", "positive control", .) %>% 
                     sub(".*fitdog.*", "FITDOG", .),
                 substrate = neutral_fp %>% 
                     sub(".*lam.*", "laminarin", .) %>% 
                     sub(".*yeastmannan.*", "yeast mannan", .),
                 stringsAsFactors = FALSE)

data_ms2 <- readMSData(files = neutral_fp, 
                   pdata = new("NAnnotatedDataFrame", 
                               pd_ms2), 
                   mode = "onDisk")

#2: pick peaks, group and fill----
data_ms2 <- findChromPeaks(data_ms2, 
                           param = cwp)

pdp <- PeakDensityParam(sampleGroups = data_ms2$substrate,
                        binSize = 0.005,
                        bw = 6) 
data_ms2 <- groupChromPeaks(data_ms2,
                            param = pdp)

data_ms2 <- fillChromPeaks(data_ms2,
                           param = fpp)

cp <- chromPeaks(data_ms2)

#format feature list----
peaksFiltered <- neutral_peaks_final

peaksFiltered$mz <- NULL
peaksFiltered$mzmin <- NULL
peaksFiltered$mzmax <- NULL
names(peaksFiltered)[names(peaksFiltered) == "i.mz"] <- "mz"
names(peaksFiltered)[names(peaksFiltered) == "i.mzmin"] <- "mzmin"
names(peaksFiltered)[names(peaksFiltered) == "i.mzmax"] <- "mzmax"

peaksFiltered<- peaksFiltered %>% 
    select(c("mz",
             "mzmin",
             "mzmax"))

#filter peaks to match identified features----

#match using foverlaps from data.table (very fast)
cp <- as.data.frame(cp)
setDT(cp)
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

cp.filtered <- unique(cp.filtered)

cp.new <- as.matrix(cp.filtered)


#reassign chromPeaks ----
data_ms2_old <- data_ms2
chromPeaks(data_ms2) <- cp.new

#extract ms2 associated with chromPeaks ----
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
                                    fcol = "peak_id",
                                    mzd = 0.005,
                                    intensityFun = mean)

#see how many spectra per sample
table(fromFile(ms2_features_comb))

#LAMINARIN FITDOG SPECTRA  ----
ms2_lamFIT <- ms2_features_comb[fromFile(ms2_features_comb) == 
                                     grep("lam_fitdog", data_ms2$name)]
ms2_lamFIT_nocomb <-  ms2_features[fromFile(ms2_features) == 
                                        grep("lam_fitdog", data_ms2$name)]

#normalise 
ms2_lamFIT_old1 <- ms2_lamFIT
ms2_lamFIT@listData <- lapply(ms2_lamFIT_old1, 
                               normalise, 
                               method = "max")

#remove low intensity peaks
ms2_lamFIT_old2 <- ms2_lamFIT
ms2_lamFIT@listData <- lapply(ms2_lamFIT_old2, 
                               removePeaks, 
                               t = 0.02)

#see precursor 
precursorMz(ms2_lamFIT)
rtime(ms2_lamFIT)/60

#extract data
ms2_lamFIT.df <- data.frame(precursorMz = as.numeric(),
                             rt = as.numeric(),
                             mz = as.numeric(),
                             intensity = as.numeric())

for (i in 1:length(ms2_lamFIT)){
    mz = sprintf("%.4f",ms2_lamFIT[[i]]@mz)
    intensity = ms2_lamFIT[[i]]@intensity * 100
    rt = rep(sprintf("%.4f", ms2_lamFIT[[i]]@rt / 60),
             length(mz))
    precursorMz = rep(sprintf("%.4f",ms2_lamFIT[[i]]@precursorMz),
                      length(mz))
    temp <- data.frame(precursorMz = precursorMz,
                       rt = rt,
                       mz = mz,
                       intensity = intensity)
    ms2_lamFIT.df <- rbind(temp,
                           ms2_lamFIT.df)
}

ms2_lamFIT.df$precursorMz <- as.numeric(ms2_lamFIT.df$precursorMz)
ms2_lamFIT.df$rt <- as.numeric(ms2_lamFIT.df$rt)
ms2_lamFIT.df$mz <- as.numeric(ms2_lamFIT.df$mz)


ms2_lamFIT_nocomb.df <- data.frame(precursorMz = as.numeric(),
                                    rt = as.numeric())

for (i in 1:length(ms2_lamFIT_nocomb)){
    rt = floor(ms2_lamFIT_nocomb[[i]]@rt / 60)
    precursorMz = round(ms2_lamFIT_nocomb[[i]]@precursorMz, 4)
    temp <- data.frame(precursorMz = precursorMz,
                       rt = rt)
    ms2_lamFIT_nocomb.df <- rbind(temp,
                                  ms2_lamFIT_nocomb.df)
}

#get how many spectra combined for each ion
ms2_lamFIT_summary <- ms2_lamFIT_nocomb.df %>% 
    group_by(precursorMz, rt) %>% 
    summarise(n = n(), .groups = "keep")

ms_lamFIT_forloop <- ms2_lamFIT.df
setDT(ms_lamFIT_forloop)
ms_lamFIT_forloop <- ms_lamFIT_forloop[order(precursorMz, rt),]
ms_lamFIT_forloop <- unique(ms_lamFIT_forloop,
                            by = c("precursorMz", "rt"))
#plot
for (i in 1:length(precursorMz(ms2_lamFIT))){
    precursorMz = ms_lamFIT_forloop$precursorMz[i]
    rt = ms_lamFIT_forloop$rt[i]
    n = ms2_lamFIT_summary$n[i]
    filename = paste0("./analysis/neutral/ms2_plots/laminarin/FITDOG/",
                      precursorMz, "_",
                      rt, "min_",
                      n, "spec")
    
    p1 <- ggplot(ms2_lamFIT.df %>% 
               filter(precursorMz == !!precursorMz &
                          rt == !!rt &
                          intensity > 0.1),
           aes(x = mz, y = intensity)) +
        geom_segment(aes(x=mz, 
                         xend=mz, 
                         y=0, 
                         yend=intensity),
                     lwd = 1.3) +
        geom_text(aes(x = mz,
                      y = intensity,
                      label = as.character(mz)),
                  angle = 90,
                  nudge_y = 7,
                  size = 5.5,
                  family = "Avenir") +
        scale_x_continuous(breaks = seq(0, 2000, by = 50),
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
    
    tiff(paste0(filename, ".tiff"),
         height = 9,
         width = 9,
         units = "in",
         res = 300)
    
    print(p1)
    
    dev.off()

    svg(paste0(filename, ".svg"),
         height = 9,
         width = 9)
    
    print(p1)
    
    dev.off()
    
}






#LAMINARIN STANDARDS SPECTRA  ----
ms2_lamPos <- ms2_features_comb[fromFile(ms2_features_comb) == 
                                    grep("lam_omix", data_ms2$name)]
ms2_lamPos_nocomb <-  ms2_features[fromFile(ms2_features) == 
                                       grep("lam_omix", data_ms2$name)]

#normalise 
ms2_lamPos_old1 <- ms2_lamPos
ms2_lamPos@listData <- lapply(ms2_lamPos_old1, 
                              normalise, 
                              method = "max")

#remove low intensity peaks
ms2_lamPos_old2 <- ms2_lamPos
ms2_lamPos@listData <- lapply(ms2_lamPos_old2, 
                              removePeaks, 
                              t = 0.02)

#see precursor 
precursorMz(ms2_lamPos)
rtime(ms2_lamPos)/60

#extract data
ms2_lamPos.df <- data.frame(precursorMz = as.numeric(),
                            rt = as.numeric(),
                            mz = as.numeric(),
                            intensity = as.numeric())

for (i in 1:length(ms2_lamPos)){
    mz = sprintf("%.4f",ms2_lamPos[[i]]@mz)
    intensity = ms2_lamPos[[i]]@intensity * 100
    rt = rep(sprintf("%.4f", ms2_lamPos[[i]]@rt / 60),
             length(mz))
    precursorMz = rep(sprintf("%.4f",ms2_lamPos[[i]]@precursorMz),
                      length(mz))
    temp <- data.frame(precursorMz = precursorMz,
                       rt = rt,
                       mz = mz,
                       intensity = intensity)
    ms2_lamPos.df <- rbind(temp,
                           ms2_lamPos.df)
}

ms2_lamPos.df$precursorMz <- as.numeric(ms2_lamPos.df$precursorMz)
ms2_lamPos.df$rt <- as.numeric(ms2_lamPos.df$rt)
ms2_lamPos.df$mz <- as.numeric(ms2_lamPos.df$mz)

breaks = seq(0, 60, 1)


ms2_lamPos_nocomb.df <- data.frame(precursorMz = as.numeric(),
                                   rt = as.numeric())

for (i in 1:length(ms2_lamPos_nocomb)){
    rt = cut(ms2_lamPos_nocomb[[i]]@rt / 60, breaks = breaks)
    precursorMz = round(ms2_lamPos_nocomb[[i]]@precursorMz, 4)
    temp <- data.frame(precursorMz = precursorMz,
                       rt = rt)
    ms2_lamPos_nocomb.df <- rbind(temp,
                                  ms2_lamPos_nocomb.df)
}

#get how many spectra combined for each ion
ms2_lamPos_summary <- ms2_lamPos_nocomb.df %>% 
    group_by(precursorMz, rt) %>% 
    summarise(n = n(), .groups = "keep")

ms_lamPos_forloop <- ms2_lamPos.df
setDT(ms_lamPos_forloop)
ms_lamPos_forloop <- ms_lamPos_forloop[order(precursorMz, rt),]
ms_lamPos_forloop <- unique(ms_lamPos_forloop,
                            by = c("precursorMz", "rt"))
#plot
for (i in 1:length(precursorMz(ms2_lamPos))){
    precursorMz = ms_lamPos_forloop$precursorMz[i]
    rt = ms_lamPos_forloop$rt[i]
    n = ms2_lamPos_summary$n[i]
    filename = paste0("./analysis/neutral/ms2_plots/laminarin/pos/",
                      precursorMz, "_",
                      rt, "min_",
                      n, "spec")
    
    p1 <- ggplot(ms2_lamPos.df %>% 
                     filter(precursorMz == !!precursorMz &
                                rt == !!rt &
                                intensity > 0.1),
                 aes(x = mz, y = intensity)) +
        geom_segment(aes(x=mz, 
                         xend=mz, 
                         y=0, 
                         yend=intensity),
                     lwd = 1.3) +
        geom_text(aes(x = mz,
                      y = intensity,
                      label = as.character(mz)),
                  angle = 90,
                  nudge_y = 7,
                  size = 5.5,
                  family = "Avenir") +
        scale_x_continuous(breaks = seq(0, 2000, by = 50),
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
    
    tiff(paste0(filename, ".tiff"),
         height = 9,
         width = 9,
         units = "in",
         res = 300)
    
    print(p1)
    
    dev.off()
    
    svg(paste0(filename, ".svg"),
        height = 9,
        width = 9)
    
    print(p1)
    
    dev.off()
    
}






#LAMINARIN MIRROR PLOTS-----
#-400 ion----

tiff("./analysis/neutral/ms2_plots/laminarin/400.2448_9.45min.tiff",
     height = 9,
     width = 9,
     units = "in",
     res = 300)

par(family = "Avenir")
plot(ms2_lamFIT[[1]], 
     ms2_lamPos[[1]])
dev.off()

svg("./analysis/neutral/ms2_plots/laminarin/400.2448_9.45min.svg",
     height = 9,
     width = 9)

par(family = "Avenir")
plot(ms2_lamFIT[[1]], 
     ms2_lamPos[[1]])
dev.off()

#-562 ion----

tiff("./analysis/neutral/ms2_plots/laminarin/562.2976_13.2min.tiff",
     height = 9,
     width = 9,
     units = "in",
     res = 300)

par(family = "Avenir")
plot(ms2_lamFIT[[3]], 
     ms2_lamPos[["CP057.F2.S02996"]])
dev.off()

svg("./analysis/neutral/ms2_plots/laminarin/562.2976_13.2min.svg",
    height = 9,
    width = 9)

par(family = "Avenir")
plot(ms2_lamFIT[[3]], 
     ms2_lamPos[["CP057.F2.S02996"]])
dev.off()





#-724 ion----

tiff("./analysis/neutral/ms2_plots/laminarin/724.3504_15.5min.tiff",
     height = 9,
     width = 9,
     units = "in",
     res = 300)

par(family = "Avenir")
plot(ms2_lamFIT[[5]], 
     ms2_lamPos[["CP074.F2.S03510"]])
dev.off()

svg("./analysis/neutral/ms2_plots/laminarin/724.3504_15.5min.svg",
    height = 9,
    width = 9)

par(family = "Avenir")
plot(ms2_lamFIT[[5]], 
     ms2_lamPos[["CP074.F2.S03510"]])
dev.off()








#-886 ion----

tiff("./analysis/neutral/ms2_plots/laminarin/886.4033_17.5min.tiff",
     height = 9,
     width = 9,
     units = "in",
     res = 300)

par(family = "Avenir")
plot(ms2_lamFIT[[6]], 
     ms2_lamPos[["CP081.F2.S03944"]])
dev.off()

svg("./analysis/neutral/ms2_plots/laminarin/886.4033_17.5min.svg",
    height = 9,
    width = 9)

par(family = "Avenir")
plot(ms2_lamFIT[[6]], 
     ms2_lamPos[["CP081.F2.S03944"]])
dev.off()








#-1048 ion----

tiff("./analysis/neutral/ms2_plots/laminarin/1048.456_19.3min.tiff",
     height = 9,
     width = 9,
     units = "in",
     res = 300)

par(family = "Avenir")
plot(ms2_lamFIT[[7]], 
     ms2_lamPos[["CP090.F2.S04314"]])
dev.off()

svg("./analysis/neutral/ms2_plots/laminarin/1048.456_19.3min.svg",
    height = 9,
    width = 9)

par(family = "Avenir")
plot(ms2_lamFIT[[7]], 
     ms2_lamPos[["CP090.F2.S04314"]])
dev.off()









#-all ----
tiff("./analysis/neutral/ms2_plots/laminarin/all.tiff",
     height = 9,
     width = 9,
     units = "in",
     res = 300)

par(family = "Avenir",
    mfrow = c(3,2),
    mar = c(4,4,1,1)) 
plot(ms2_lamFIT[[1]], 
     ms2_lamPos[[1]])
plot(ms2_lamFIT[[3]], 
     ms2_lamPos[["CP057.F2.S02996"]])
plot(ms2_lamFIT[[5]], 
     ms2_lamPos[["CP074.F2.S03510"]])
plot(ms2_lamFIT[[6]], 
     ms2_lamPos[["CP081.F2.S03944"]])
plot(ms2_lamFIT[[7]], 
     ms2_lamPos[["CP090.F2.S04314"]])
dev.off()

svg("./analysis/neutral/ms2_plots/laminarin/all.svg",
     height = 9,
     width = 9)

par(family = "Avenir",
    mfrow = c(3,2),
    mar = c(4,4,1,1)) 
plot(ms2_lamFIT[[1]], 
     ms2_lamPos[[1]])
plot(ms2_lamFIT[[3]], 
     ms2_lamPos[["CP057.F2.S02996"]])
plot(ms2_lamFIT[[5]], 
     ms2_lamPos[["CP074.F2.S03510"]])
plot(ms2_lamFIT[[6]], 
     ms2_lamPos[["CP081.F2.S03944"]])
plot(ms2_lamFIT[[7]], 
     ms2_lamPos[["CP090.F2.S04314"]])
dev.off()



#YEAST MANNAN FITDOG SPECTRA  ----
ms2_ymFIT <- ms2_features_comb[fromFile(ms2_features_comb) == 
                                    grep("yeastmannan_fitdog", data_ms2$name)]
ms2_ymFIT_nocomb <-  ms2_features[fromFile(ms2_features) == 
                                       grep("yeastmannan_fitdog", data_ms2$name)]

#normalise 
ms2_ymFIT_old1 <- ms2_ymFIT

ms2_ymFIT@listData <- lapply(ms2_ymFIT_old1, 
                              normalise, 
                              method = "max")

#remove low intensity peaks
ms2_ymFIT_old2 <- ms2_ymFIT
ms2_ymFIT@listData <- lapply(ms2_ymFIT_old2, 
                              removePeaks, 
                              t = 0.02)

#see precursor 
precursorMz(ms2_ymFIT)
rtime(ms2_ymFIT)/60

#extract data
ms2_ymFIT.df <- data.frame(precursorMz = as.numeric(),
                            rt = as.numeric(),
                            mz = as.numeric(),
                            intensity = as.numeric())

for (i in 1:length(ms2_ymFIT)){
    mz = sprintf("%.4f",ms2_ymFIT[[i]]@mz)
    intensity = ms2_ymFIT[[i]]@intensity * 100
    rt = rep(sprintf("%.4f", ms2_ymFIT[[i]]@rt / 60),
             length(mz))
    precursorMz = rep(sprintf("%.4f",ms2_ymFIT[[i]]@precursorMz),
                      length(mz))
    temp <- data.frame(precursorMz = precursorMz,
                       rt = rt,
                       mz = mz,
                       intensity = intensity)
    ms2_ymFIT.df <- rbind(temp,
                           ms2_ymFIT.df)
}

ms2_ymFIT.df$precursorMz <- as.numeric(ms2_ymFIT.df$precursorMz)
ms2_ymFIT.df$rt <- as.numeric(ms2_ymFIT.df$rt)
ms2_ymFIT.df$mz <- as.numeric(ms2_ymFIT.df$mz)


ms2_ymFIT_nocomb.df <- data.frame(precursorMz = as.numeric(),
                                   rt = as.numeric())

for (i in 1:length(ms2_ymFIT_nocomb)){
    rt = floor(ms2_ymFIT_nocomb[[i]]@rt / 60)
    precursorMz = round(ms2_ymFIT_nocomb[[i]]@precursorMz, 4)
    temp <- data.frame(precursorMz = precursorMz,
                       rt = rt)
    ms2_ymFIT_nocomb.df <- rbind(temp,
                                  ms2_ymFIT_nocomb.df)
}

#get how many spectra combined for each ion
ms2_ymFIT_summary <- ms2_ymFIT_nocomb.df %>% 
    group_by(precursorMz, rt) %>% 
    summarise(n = n(), .groups = "keep")

ms_ymFIT_forloop <- ms2_ymFIT.df
setDT(ms_ymFIT_forloop)
ms_ymFIT_forloop <- ms_ymFIT_forloop[order(precursorMz, rt),]
ms_ymFIT_forloop <- unique(ms_ymFIT_forloop,
                            by = c("precursorMz", "rt"))
#plot
for (i in 1:length(precursorMz(ms2_ymFIT))){
    precursorMz = ms_ymFIT_forloop$precursorMz[i]
    rt = ms_ymFIT_forloop$rt[i]
    n = ms2_ymFIT_summary$n[i]
    filename = paste0("./analysis/neutral/ms2_plots/yeastmannan/FITDOG/",
                      precursorMz, "_",
                      rt, "min_",
                      n, "spec")
    
    p1 <- ggplot(ms2_ymFIT.df %>% 
                     filter(precursorMz == !!precursorMz &
                                rt == !!rt &
                                intensity > 0.1),
                 aes(x = mz, y = intensity)) +
        geom_segment(aes(x=mz, 
                         xend=mz, 
                         y=0, 
                         yend=intensity),
                     lwd = 1.3) +
        geom_text(aes(x = mz,
                      y = intensity,
                      label = as.character(mz)),
                  angle = 90,
                  nudge_y = 7,
                  size = 5.5,
                  family = "Avenir") +
        scale_x_continuous(breaks = seq(0, 2000, by = 50),
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
    
    tiff(paste0(filename, ".tiff"),
         height = 9,
         width = 9,
         units = "in",
         res = 300)
    
    print(p1)
    
    dev.off()
    
    svg(paste0(filename, ".svg"),
        height = 9,
        width = 9)
    
    print(p1)
    
    dev.off()
    
}






#YEAST MANNAN GH76 SPECTRA  ----
ms2_ymPos <- ms2_features_comb[fromFile(ms2_features_comb) == 
                                   grep("yeastmannan_gh76", data_ms2$name)]
ms2_ymPos_nocomb <-  ms2_features[fromFile(ms2_features) == 
                                      grep("yeastmannan_gh76", data_ms2$name)]

#normalise 
ms2_ymPos_old1 <- ms2_ymPos

ms2_ymPos@listData <- lapply(ms2_ymPos_old1, 
                             normalise, 
                             method = "max")

#remove low intensity peaks
ms2_ymPos_old2 <- ms2_ymPos
ms2_ymPos@listData <- lapply(ms2_ymPos_old2, 
                             removePeaks, 
                             t = 0.02)

#see precursor 
precursorMz(ms2_ymPos)
rtime(ms2_ymPos)/60

#extract data
ms2_ymPos.df <- data.frame(precursorMz = as.numeric(),
                           rt = as.numeric(),
                           mz = as.numeric(),
                           intensity = as.numeric())

for (i in 1:length(ms2_ymPos)){
    mz = sprintf("%.4f",ms2_ymPos[[i]]@mz)
    intensity = ms2_ymPos[[i]]@intensity * 100
    rt = rep(sprintf("%.4f", ms2_ymPos[[i]]@rt / 60),
             length(mz))
    precursorMz = rep(sprintf("%.4f",ms2_ymPos[[i]]@precursorMz),
                      length(mz))
    temp <- data.frame(precursorMz = precursorMz,
                       rt = rt,
                       mz = mz,
                       intensity = intensity)
    ms2_ymPos.df <- rbind(temp,
                          ms2_ymPos.df)
}

ms2_ymPos.df$precursorMz <- as.numeric(ms2_ymPos.df$precursorMz)
ms2_ymPos.df$rt <- as.numeric(ms2_ymPos.df$rt)
ms2_ymPos.df$mz <- as.numeric(ms2_ymPos.df$mz)


ms2_ymPos_nocomb.df <- data.frame(precursorMz = as.numeric(),
                                  rt = as.numeric())

for (i in 1:length(ms2_ymPos_nocomb)){
    rt = floor(ms2_ymPos_nocomb[[i]]@rt / 60)
    precursorMz = round(ms2_ymPos_nocomb[[i]]@precursorMz, 4)
    temp <- data.frame(precursorMz = precursorMz,
                       rt = rt)
    ms2_ymPos_nocomb.df <- rbind(temp,
                                 ms2_ymPos_nocomb.df)
}

#get how many spectra combined for each ion
ms2_ymPos_summary <- ms2_ymPos_nocomb.df %>% 
    group_by(precursorMz, rt) %>% 
    summarise(n = n(), .groups = "keep")

ms_ymPos_forloop <- ms2_ymPos.df
setDT(ms_ymPos_forloop)
ms_ymPos_forloop <- ms_ymPos_forloop[order(precursorMz, rt),]
ms_ymPos_forloop <- unique(ms_ymPos_forloop,
                           by = c("precursorMz", "rt"))
#plot
for (i in 1:length(precursorMz(ms2_ymPos))){
    precursorMz = ms_ymPos_forloop$precursorMz[i]
    rt = ms_ymPos_forloop$rt[i]
    n = ms2_ymPos_summary$n[i]
    filename = paste0("./analysis/neutral/ms2_plots/yeastmannan/pos/",
                      precursorMz, "_",
                      rt, "min_",
                      n, "spec")
    
    p1 <- ggplot(ms2_ymPos.df %>% 
                     filter(precursorMz == !!precursorMz &
                                rt == !!rt &
                                intensity > 0.1),
                 aes(x = mz, y = intensity)) +
        geom_segment(aes(x=mz, 
                         xend=mz, 
                         y=0, 
                         yend=intensity),
                     lwd = 1.3) +
        geom_text(aes(x = mz,
                      y = intensity,
                      label = as.character(mz)),
                  angle = 90,
                  nudge_y = 7,
                  size = 5.5,
                  family = "Avenir") +
        scale_x_continuous(breaks = seq(0, 2000, by = 50),
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
    
    tiff(paste0(filename, ".tiff"),
         height = 9,
         width = 9,
         units = "in",
         res = 300)
    
    print(p1)
    
    dev.off()
    
    svg(paste0(filename, ".svg"),
        height = 9,
        width = 9)
    
    print(p1)
    
    dev.off()
    
}




#YEAST MANNAN  MIRROR PLOTS-----
#-400 ion----

tiff("./analysis/neutral/ms2_plots/yeastmannan/400.2448_9.45min.tiff",
     height = 9,
     width = 9,
     units = "in",
     res = 300)

par(family = "Avenir")
plot(ms2_ymFIT[[2]], 
     ms2_ymPos[[2]])
dev.off()

svg("./analysis/neutral/ms2_plots/yeastmannan/400.2448_9.45min.svg",
    height = 9,
    width = 9)

par(family = "Avenir")
plot(ms2_ymFIT[[2]], 
     ms2_ymPos[[2]])
dev.off()

#-562 ion----

tiff("./analysis/neutral/ms2_plots/yeastmannan/562.2976_14.2min.tiff",
     height = 9,
     width = 9,
     units = "in",
     res = 300)

par(family = "Avenir")
plot(ms2_ymFIT[[4]], 
     ms2_ymPos[[5]])
dev.off()

svg("./analysis/neutral/ms2_plots/yeastmannan/562.2976_14.2min.svg",
    height = 9,
    width = 9)

par(family = "Avenir")
plot(ms2_ymFIT[[4]], 
     ms2_ymPos[[5]])
dev.off()





#-724 ion----

tiff("./analysis/neutral/ms2_plots/yeastmannan/724.3504_17.7min.tiff",
     height = 9,
     width = 9,
     units = "in",
     res = 300)

par(family = "Avenir")
plot(ms2_ymFIT[[5]], 
     ms2_ymPos[[6]])
dev.off()

svg("./analysis/neutral/ms2_plots/yeastmannan/724.3504_17.7min.svg",
    height = 9,
    width = 9)

par(family = "Avenir")
plot(ms2_ymFIT[[5]], 
     ms2_ymPos[[6]])
dev.off()







#-886 ion----

tiff("./analysis/neutral/ms2_plots/yeastmannan/886.4033_20.5min.tiff",
     height = 9,
     width = 9,
     units = "in",
     res = 300)

par(family = "Avenir")
plot(ms2_ymFIT[[6]], 
     ms2_ymPos[[7]])
dev.off()

svg("./analysis/neutral/ms2_plots/yeastmannan/886.4033_20.5min.svg",
    height = 9,
    width = 9)

par(family = "Avenir")
plot(ms2_ymFIT[[6]], 
     ms2_ymPos[[7]])
dev.off()









#-all ----

tiff("./analysis/neutral/ms2_plots/yeastmannan/all.tiff",
     height = 9,
     width = 9,
     units = "in",
     res = 300)

par(family = "Avenir",
    mfrow = c(2,2),
    mar = c(4,4,1,1)) 
plot(ms2_ymFIT[[2]], 
     ms2_ymPos[[2]])
plot(ms2_ymFIT[[4]], 
     ms2_ymPos[[5]])
plot(ms2_ymFIT[[5]], 
     ms2_ymPos[[6]])
plot(ms2_ymFIT[[6]], 
     ms2_ymPos[[7]])
dev.off()

svg("./analysis/neutral/ms2_plots/yeastmannan/all.svg",
     height = 9,
     width = 9)

par(family = "Avenir",
    mfrow = c(2,2),
    mar = c(4,4,1,1)) 
plot(ms2_ymFIT[[2]], 
     ms2_ymPos[[2]])
plot(ms2_ymFIT[[4]], 
     ms2_ymPos[[5]])
plot(ms2_ymFIT[[5]], 
     ms2_ymPos[[6]])
plot(ms2_ymFIT[[6]], 
     ms2_ymPos[[7]])
dev.off()






