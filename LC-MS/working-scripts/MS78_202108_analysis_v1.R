setwd("/Users/margotbligh/Google_Drive/MPI_PhD/Lab-things/FITDOG/methods-dev/MS78_202108")
load("./analysis/RData/RData_20210905.RData")

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

#2. Import and inspect MS data --------------------------------------------------------
#get file paths to mzML files
fp <- dir(path = "mzML-files", 
          all.files = FALSE, 
          full.names = TRUE,
          recursive = TRUE,
          include.dirs = FALSE)



#create phenodata data.frame
#each sample must have a unique name!
pd <- data.frame(name = basename(fp) %>%
                     sub("MS78_", "", .) %>% 
                     sub(".mzML", "", .),
                 sample_type = basename(fp) %>% 
                     sub(".*_0min_.*|.*_h2o_.*", "neg", .) %>% 
                     sub(".*_60min_.*", "FITDOG", .) %>% 
                     sub(".*blank.*", "solvent blank", ., ignore.case = T) %>% 
                     sub(".*HILIC.*", "HILIC standard", .) %>% 
                     sub(".*QC.*", "QC", .),
                 substrate = basename(fp) %>% 
                     sub(".*_h2o_.*", "water", .) %>% 
                     sub(".*_ves_.*", "fucoidan", .) %>% 
                     sub(".*bic.*", "laminarin", .) %>% 
                     sub("MS78.*", NA, .),
                 reduced = basename(fp) %>% 
                     sub(".*NaBH4_.*", "y", .) %>% 
                     sub(".*\\dmin_.*", "n", .) %>% 
                     sub("MS78.*", NA, .),
                 reactiontime = basename(fp) %>% 
                     sub(".*_60min_.*", 60, .) %>% 
                     sub(".*0min_.*", 0, .) %>% 
                     sub("MS78.*", NA, .),
                 group = basename(fp) %>% 
                     sub(".*blank.*", "solvent blank", ., ignore.case = T) %>% 
                     sub(".*QC.*", "QC", .) %>% 
                     sub(".*HILIC.*", "HILIC standard", .) %>% 
                     sub("\\+NaBH4.*", "reduced", .) %>% 
                     sub("60min", "t=60", .) %>% 
                     sub("0min", "t=0", .) %>% 
                     sub(".*h2o", "water", .) %>% 
                     sub(".*ves", "fucoidan", .) %>% 
                     sub(".*bic", "laminarin", .) %>% 
                     sub(".mzML|_\\d{1,2}.mzML", "", .) %>% 
                     gsub("_", " ", .),
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
dir.create("./analysis/processing_plots/peakpicking",
           showWarnings = FALSE)

#set parameters
# cwp<-CentWaveParam()
# cwp@ppm<-1.6
# cwp@peakwidth<-c(10,100)
# cwp@snthresh<-5

cwp<-CentWaveParam()
cwp@ppm<-5
cwp@peakwidth<-c(10,60)
cwp@snthresh<-5
cwp@noise <- 5000
cwp@prefilter <- c(5, 1000)
cwp@mzdiff <- 0.001

#pick peaks
data_peaks<-findChromPeaks(data, 
                           param=cwp)

error = 0.0025
pal1 <- hcl.colors(n = length(pd$name),
                   palette = "Dark3")
groups <- data_peaks$group %>% unique()

#hex-1 m+cl
chr1 <- chromatogram(data_peaks,
                     mz = c(215.03279157990	 - error,
                            215.03279157990	 + error))


cairo_pdf("./analysis/processing_plots/peakpicking/peakpicking_5ppm_5sn_pw10to60_noise5000_prefilter5-1000_mzdiff0.001_hex-1_cl-adduct.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mfrow=c(5,3))
for (i in 1:length(groups)){
    plot(chr1[,chr1$group == groups[i]],
         lwd = 2,
         cex.main = 1,
         peakCol = "black",
         peakType = "rectangle",
         xlim = c(50, 200),
         ylim = c(0,1e5),
         main = groups[i]
    )
}
dev.off()

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

names(pal1) <- data$name

cairo_pdf("./analysis/processing_plots/peakgrouping/peakgrouping_binsize0.005_bw6_hex-1_cl-adduct.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mar=c(4,4,3,10))
plotChromPeakDensity(chr1, 
                     col = pal1, 
                     param = pdp,
                     peakBg = pal1[chromPeaks(chr1)[, "sample"]],
                     peakCol = pal1[chromPeaks(chr1)[, "sample"]],
                     peakPch = 16)
legend("topright",
       legend = paste0(seq(1,length(pal1),1),
                       "=",
                       names(pal1)),
       inset=c(0,0),
       fill = pal1,
       pt.cex = 0.3,
       cex = 0.5,
       bty = "n",
       horiz = FALSE,
       xpd=TRUE,
       ncol = 3)
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

#8: PCA plot: filled vs not filled features ----
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

#essentially no separation into groups, only 20% of variation explained...


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
tiff("./analysis/processing_plots/pca/pca-filled-featureIntensities.tiff",
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
ft_ints.fl <- ft_ints.fl[,grep("HILIC", names(ft_ints.fl), invert = TRUE)]
#perform PCA with the intensities mean centered
pc2 <- prcomp(t(ft_ints.fl), center = TRUE)
#plot
tiff("./analysis/processing_plots/pca/pca-filled-featureIntensities-noHILICstd.tiff",
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
     bg = pal_group[pd$group[pd$group!="HILIC standard"]], 
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
text(x = pc2$x[, 1],
     y = pc2$x[, 2],
     labels = names(ft_ints.fl),
     col = "black",
     pos = 1,
     cex = 0.5)
dev.off()

#some very separated, all others in one clump
#remove very separated and run again

x <- c("20210824_32_bic_60min_39",
       "20210824_QC_3",
       "20210826_QC_3",
       "20210824_Blank_1",
       "20210826_blank_1")


#subset data
ft_ints.fl2 <- ft_ints.fl %>% dplyr::select(-x)
#make new palette
pd2 <- pd[!pd$name %in% x & pd$group != "HILIC standard",]
pal_group2 <- hcl.colors(n = length(unique(pd2$group)),
                         palette = "Dark3")
names(pal_group2) <- unique(pd2$group)
#perform PCA with the intensities mean centered
pc2 <- prcomp(t(ft_ints.fl2), center = TRUE)
#plot
tiff("./analysis/processing_plots/pca/pca-filled-featureIntensities-noHILICstd-subset.tiff",
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
     bg = pal_group[pd2$group], 
     cex = 2)
grid()
legend("topleft",
       legend = names(pal_group2) %>% unique(),
       pt.bg = pal_group2 %>% unique(),
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
text(x = pc2$x[, 1],
     y = pc2$x[, 2],
     labels = names(ft_ints.fl2),
     col = "black",
     pos = 1,
     cex = 0.5)
dev.off()

#clustering appears quite random... no separation into groups
#continue with filled data for now?

res <- data_peaks_grouped_filled
rm(data, data_peaks_grouped, data_peaks, data_peaks_grouped_filled)

#9: Save diffreport of xdata -----
xset <- as(res, "xcmsSet")
sampnames(xset) <- pData(res)$name
sampclass(xset) <- pData(res)$group

#10. Isotope picking----
##create xsannotate object
#extracts the peaktable from a provided xcmsSet,
#which is used for all further analysis

an <- xsAnnotate(xset)

##Group peaks of a xsAnnotate object according to their retention time 
#Perfwhm = parameter defines the window width, which is used for matching
an <- groupFWHM(an, 
                 perfwhm = 0.6)

##Annotate isotope peaks
#Mzabs = the allowed m/z error
an <- findIsotopes(an, 
                    mzabs=0.01)

##Peak grouping after correlation information into pseudospectrum groups 
#cor_eic_th = correlation threshold for EIC correlation
an <- groupCorr(an, 
                  cor_eic_th=0.75)

##Find adducts
an <- findAdducts(an, 
                    polarity="negative")


#16: Compare samples with heatmap ----
##ALL PEAKS - laminarin
#add unique rownames
setDF(pl)
x <- paste0("rt_",
            round(pl$rt/60,2),
            "mz",
            round(pl$mz,4))
rownames(pl) <- x

#subset to only have intensity counts
sampleColNames <- pd$name %>% sub("\\+", ".", .) %>% sub("^2", "X2", .)
sampleColNames <- sampleColNames[grep("_bic_", sampleColNames)]
sampleColNames <- sampleColNames[-grep("NaBH4", sampleColNames)]

counts_all <- pl %>% 
    dplyr::select(sampleColNames)

#change NA to 0
counts_all[is.na(counts_all)] <- 0

#set factor level
group <- pd$group[pd$substrate=="laminarin" & pd$reduced == "n"]
group <- group[!is.na(group)]
group<-factor(group)

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
selY_n_all <- dge_n_all[rownames(result_n_all$table)[result_n_all$table$PValue<0.05 & 
                                                         result_n_all$table$logFC < -2 |
                                                         result_n_all$table$PValue<0.05 & 
                                                         result_n_all$table$logFC > 2],]


#make heatmap
cimColour <- rev(viridis(1000))

par(mar= c(20, 15, 15, 40))
cim(selY_n_all, 
    color = cimColour,
    symkey=FALSE,
    mar=c(10,20),
    col.cex = 0.8,
    #row.sideColors = cimColourRows,
    row.names = TRUE,
    #row.cex = 0.8,
    #keysize = c(0.1, 0.1)
    #save = "tiff",
    #name.save = "./analysis/heatmap/heatmap_matched_unmatched.tiff"
)
dev.off()


##ALL PEAKS - laminarin - reduced
#subset to only have intensity counts
sampleColNames <- pd$name %>% sub("\\+", ".", .) %>% sub("^2", "X2", .)
sampleColNames <- sampleColNames[grep("_bic_", sampleColNames)]
sampleColNames <- sampleColNames[grep("NaBH4", sampleColNames)]

counts_all <- pl %>% 
    dplyr::select(sampleColNames)

#change NA to 0
counts_all[is.na(counts_all)] <- 0

#set factor level
group <- pd$group[pd$substrate=="laminarin" & pd$reduced == "y"]
group <- group[!is.na(group)]
group<-factor(group)

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
selY_n_all <- dge_n_all[rownames(result_n_all$table)[result_n_all$table$PValue<0.05 & 
                                                         result_n_all$table$logFC < -2 |
                                                         result_n_all$table$PValue<0.05 & 
                                                         result_n_all$table$logFC > 2],]


#make heatmap
cimColour <- rev(viridis(1000))

par(mar= c(20, 15, 15, 40))
cim(selY_n_all, 
    color = cimColour,
    symkey=FALSE,
    mar=c(10,20),
    col.cex = 0.8,
    #row.sideColors = cimColourRows,
    row.names = TRUE,
    #row.cex = 0.8,
    #keysize = c(0.1, 0.1)
    #save = "tiff",
    #name.save = "./analysis/heatmap/heatmap_matched_unmatched.tiff"
)
dev.off()


lam_reduced <- rownames(result_n_all$table)[result_n_all$table$PValue<0.05 & 
                                                result_n_all$table$logFC > 2]
lam_reduced <- data.frame(mz = as.numeric(sub(".*mz","",lam_reduced)),
                          rt = as.numeric(gsub("rt_|mz.*","",lam_reduced)))
mz_predicted <- fread("neutral-predicted-sugars_dp1to10_alditol_unsaturated_pentose.txt")

#remove "extra" columns
extraCol <- c('mass',
              'formula')

mz_predicted <- mz_predicted %>% 
    dplyr::select(-all_of(extraCol))

#make long format
predicted <- mz_predicted %>% 
    gather(key = "ion",
           value = "mz",
           -name,
           -dp)

#make data.table
setDT(predicted)
setDT(lam_reduced)
#create interval to overlap with (same width as for peak grouping)
predicted$mz <- as.numeric(predicted$mz)
predicted$mzmin <- predicted$mz-0.005
predicted$mzmax <- predicted$mz+0.005
lam_reduced$mzmin <- lam_reduced$mz-0.005
lam_reduced$mzmax <- lam_reduced$mz+0.005

#remove NA rows
predicted <- na.omit(predicted)

#match using foverlaps from data.table (very fast)
setkey(predicted, mzmin, mzmax)
lam_reduced <- foverlaps(lam_reduced,
                         predicted)

#change NA values created during matching (features with no match) to be blank
#remove extra columns
lam_reduced$mzmin <- NULL
lam_reduced$mzmax <- NULL

lam_reduced <- lam_reduced %>% 
    replace_na(list("name"="",
                    "ion"= "", 
                    "mz" = "",
                    "dp" = ""))
#only keep matched features
peaks_matched <- peaks[!peaks$name=="",]

#order by retention time
peaks_matched <- peaks_matched[order(rt_min),]

#make id and ion column
peaks_matched$id_ion <- paste0(peaks_matched$name, 
                               ": ",
                               peaks_matched$ion)



##ALL PEAKS - fucoidanm
#add unique rownames
setDF(pl)
x <- paste0("rt_",
            round(pl$rt/60,2),
            "mz",
            round(pl$mz,4))
rownames(pl) <- x

#subset to only have intensity counts
sampleColNames <- pd$name %>% sub("\\+", ".", .) %>% sub("^2", "X2", .)
sampleColNames <- sampleColNames[grep("_ves_", sampleColNames)]
sampleColNames <- sampleColNames[-grep("NaBH4", sampleColNames)]

counts_all <- pl %>% 
    dplyr::select(sampleColNames)

#change NA to 0
counts_all[is.na(counts_all)] <- 0

#set factor level
group <- pd$group[pd$substrate=="fucoidan" & pd$reduced == "n"]
group <- group[!is.na(group)]
group<-factor(group)

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

par(mar= c(20, 15, 15, 40))
cim(selY_n_all, 
    color = cimColour,
    symkey=FALSE,
    mar=c(10,20),
    col.cex = 0.8,
    #row.sideColors = cimColourRows,
    row.names = TRUE,
    #row.cex = 0.8,
    #keysize = c(0.1, 0.1)
    #save = "tiff",
    #name.save = "./analysis/heatmap/heatmap_matched_unmatched.tiff"
)
dev.off()

#positive fold change means higher at 60 min than 0 min

ves.df <- 





#11. Peak list filtering and formatting----
#get peak list
pl <-getPeaklist(an)

#filter by blank exclusion (detected peaks)
pl_be <-pl[pl$negative.control==0 & 
               pl$solvent.blank==0,]

#make rownames from rt and mz of features
rownames(pl_be)<-paste(round(pl_be$rt,1),
                       round(pl_be$mz,3),
                       sep="_")
#change NA to 0
pl_be[is.na(pl_be)] <- 0

#change name
peaks <- pl_be

#add rounded retention time as first colum
peaks <- cbind(rt_min = round(peaks$rt/60, 
                              1),
               peaks)

#12: Collapse features with multiple isotopes -----
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

#13: Annotate features based on predictions----
#import prediction table
#see https://github.com/margotbligh/sugarMassesPredict
#created with: sugarMassesPredict.py -dp 1 7 -p 0 -m sulphate unsaturated -i neg -s 175 1400
mz_predicted <- fread("neutral-predicted-sugars_dp1to10_alditol_unsaturated_pentose.txt")

#remove "extra" columns
extraCol <- c('mass',
              'formula')

mz_predicted <- mz_predicted %>% 
    dplyr::select(-all_of(extraCol))

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
predicted$mzmin <- predicted$mz-0.005
predicted$mzmax <- predicted$mz+0.005

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

#order by retention time
peaks_matched <- peaks_matched[order(rt_min),]

#make id and ion column
peaks_matched$id_ion <- paste0(peaks_matched$name, 
                               ": ",
                               peaks_matched$ion)


#write to table
fwrite(peaks_matched,
       "./analysis/analysis_tables/matched-peaks-v1.txt",
       sep = "\t")

fwrite(peaks,
       "./analysis/analysis_tables/peaks-v1.txt",
       sep = "\t")

#14: Extract and format eic to check peaks of identified sugars -----
#get vectors
mz.found.vector <- peaks_matched$i.mz %>% 
    round(., 3)

ions.found.vector  <- peaks_matched$id_ion

ions.rt.found.vector <- paste0(ions.found.vector,
                               "_rt",
                               peaks_matched$rt_min)

    
#remove solvent blanks and HILIC standard
res2 <- filterFile(res,
                    file = which(!grepl("standard|Solvent",
                                       res$name)))

#get phenodata vectors
res2.names <- res2$name
res2.groups <- res2$sample_type


#extract chromatograms
chr_list <- list()
error = 0.001

for (i in 1:length(mz.found.vector)){
    mzr = c(mz.found.vector[i] - error,
            mz.found.vector[i] + error)
    chr_list[[i]] <- chromatogram(res2, 
                                  mz = mzr)
}

#save list
save(chr_list,
     file = "./analysis/RData/chr_list.RData")

#extract intensity and rt values
chr_int_list <- list()
for (i in 1:length(res2.names)){
    chr_int_list[[i]] <- lapply(chr_list, function(x) {
        x[[i]]@intensity
    }) 
}

chr_rt_list <- list()
for (i in 1:length(res2.names)){
    chr_rt_list[[i]] <- lapply(chr_list, function(x) {
        x[[i]]@rtime
    }) 
}

#build data frame (long format)
res2.df <- data.frame(ion = as.character(),
                      sample = as.character(),
                      group = as.character(),
                      rt = as.numeric(),
                      intensity = as.numeric())

for (i in 1:length(res2.names)){
    for (j in 1:length(mz.found.vector)){
        rt = chr_rt_list[[i]][[j]]/60
        intensity = chr_int_list[[i]][[j]]
        sample = rep(res2.names[i], length(rt))
        group = rep(res2.groups[i], length(rt))
        ion = rep(ions.rt.found.vector[j], length(rt))
        temp <- data.frame(ion = ion,
                           sample = sample,
                           group = group,
                           rt = rt,
                           intensity = intensity)
        res2.df <- rbind(res2.df,
                         temp)
    }
}



#set NA to 0
res2.df[is.na(res2.df)] <- 0

#set variables as factors
res2.df$ion <- factor(res2.df$ion,
                      levels = unique(ions.rt.found.vector))

res2.df$group <- factor(res2.df$group,
                           levels = unique(res2.groups))

#write to file
fwrite(res2.df,
       file = "./analysis/analysis_tables/eic-table_v1.txt",
       sep = "\t")

#15: Plot EIC ----
#make directory
dir.create("./analysis/analysis_plots/eic_checking",
           showWarnings = FALSE)

#make palette
pal2 <- hcl.colors(n = length(unique(res2.groups)),
                   palette = "Dark3")

#plot for each ion the full and restricted rt chromatograms by group
for (i in 1:length(unique(ions.rt.found.vector))){
    ion = ions.rt.found.vector[i]
    rt = ions.rt.found.vector[i] %>% 
        sub(".*_rt", "", .) %>% 
        as.numeric()
    rtmin = rt - 2.5
    if (rtmin < 0) {rtmin <- 0}
    rtmax = rt + 2.5
    
    df1 <- res2.df %>% 
        filter(between(rt, rtmin, rtmax)) %>% 
        filter(ion == !!ion)
    
    df2 <- res2.df %>% 
        filter(ion == !!ion)
    
    #plot restricted retention time
    p1 <- ggplot() +
        geom_line(mapping = aes(rt,
                                intensity,
                                colour = group,
                                group = sample),
                  data = df1,
                  lwd = 1.2) +
        scale_colour_manual(values = pal2) +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(group)) +
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
                                colour = group,
                                group = sample),
                  data = df2,
                  lwd = 1.2) +
        scale_colour_manual(values = pal2) +
        labs(x= "Retention time (min)",
             y = "Intensity (a.u.)") +
        facet_grid(rows = vars(group)) +
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
    combined <- p1 + p2 & theme(legend.position = "bottom")
    
    #save as pdf
    cairo_pdf(paste0("./analysis/analysis_plots/eic_checking/",
                     ion,
                     ".pdf"),
              width = 12,
              height = 9)
    print(combined +
              plot_layout(ncol=2, guides = "collect") + 
              plot_annotation(title = ion))
    
    dev.off()
}

#plot just the EIC for hex-1-sulphate [M-H]- in standards
ion = "hex-1-sulphate-1: [M-H]-_rt2.3" #picked one rt randomly

df1 <- res2.df %>% 
    filter(between(rt, 0, 5)) %>% 
    filter(ion == !!ion) %>% 
    filter(group == "standard")

pal2 <- hcl.colors(n = length(pd$name[pd$sample_type == "standard"]),
                   palette = "Dark3")
names(pal2) <- pd$name[pd$sample_type == "standard"]

cairo_pdf("./analysis/analysis_plots/mannose-sulphate_standards_eic.pdf",
          width = 12,
          height = 9)
ggplot() +
    geom_line(mapping = aes(rt,
                            intensity,
                            colour = sample),
              data = df1,
              lwd = 1.2) +
    scale_colour_manual(values = pal2) +
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







