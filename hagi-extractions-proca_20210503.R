dir <- "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/CarbX/GCC-SPE_ProcA_tests"

setwd(dir)
#load("analysis/RData/RData_20210218.RData")

#1: Install packages --------------------------------------------------------
library(BiocStyle)
library(xcms)
library(faahKO)
library(pander)
library(RColorBrewer)
library(magrittr)
library(pheatmap)
library(MSnbase)
library(msdata)
library(png)
library(IPO)
library(tidyr)
library(detect)
library(devtools)
library(MetMSLine)
library(pcaMethods)
library(statTarget)
library(randomForest)
library(rlist)
library(purrr)
library(wesanderson)
library(ggplot2)
library(reshape2)
library(extrafont)
library(Rmisc)
library(edgeR)
library(limma)
library(mixOmics)
library(HTSFilter)
library(rstatix)
library(reshape2)
library(scales)
library(data.table)
library(remotes)
library(ggridges)
library(gridExtra)
library(tidyverse)
library(ggpubr)
library(viridis)
library(lemon)
library(cowplot)
library(ggsci)
library(ggfortify)
library(ropls)
library(gplots)
library(grid)
library(dplyr)


#2. Import and inspect MS data --------------------------------------------------------
#get file paths to mzML files
fp <- dir(path = "./mzML-files", 
                      all.files = FALSE, 
                      full.names = TRUE)
fp <- fp[-grepl("Icon", fp)]

#create phenodata data.frame
pd <- data.frame(name = basename(fp) %>%
                     sub(".*Blank", "SolventBlank", .) %>% 
                     sub(".*Pos531_2_11", "Poseidon_14h", .) %>% 
                     sub(".*Pos531_2_33", "Poseidon_02h", .) %>% 
                     sub(".*Std_1000ng_oligo_mix", "StandardMixGCCSPE", .) %>% 
                     sub(".*Std_POS_1000ng_oligo_mix", "StandardMix", .) %>% 
                     sub("Std.*","ExtractionStandard", .) %>% 
                     sub("blk.*", "ExtractionBlank", .) %>% 
                     sub("a_\\d{2}", "a", .) %>% 
                     sub("b_\\d{2}", "b", .) %>% 
                     sub("c_\\d{2}", "c", .) %>% 
                     sub(".mzML", "", .) %>% 
                     sub(".*procA_", "", .),
                 sample_type = basename(fp) %>%
                     sub(".*Blank", "SolventBlank", .) %>% 
                     sub(".*Pos531_2_11.*", "Poseidon_14h", .) %>% 
                     sub(".*Pos531_2_33.*", "Poseidon_02h", .) %>% 
                     sub(".*Std_1000ng_oligo_mix", "StandardMixGCCSPE", .) %>% 
                     sub(".*Std_POS_1000ng_oligo_mix", "StandardMix", .) %>% 
                     sub(".*col[123]_Std","ExtractionStandard", .) %>% 
                     sub(".*col[123]_blk.*", "ExtractionBlank", .) %>% 
                     sub("_\\d{2}.mzML", "", .),
                 stringsAsFactors = FALSE)

#read in data
data <- readMSData(files = fp, 
                   pdata = new("NAnnotatedDataFrame", 
                               pd), 
                   mode = "onDisk")

#split into only ms1, and ms1 and ms2
#xset only works with ms1 data
#do all steps until then for both, 
#then reassign chrompeaks for ms2 after filtering

data_ms2 <- data
data_ms1 <- data[data@featureData@data$msLevel == 1]

#3: Create initial output directories -------------------------------------
dir.create("./analysis",
           showWarnings = FALSE)
dir.create("./analysis/RData",
           showWarnings = FALSE)
dir.create("./analysis/processing_plots",
           showWarnings = FALSE)
dir.create("./analysis/analysis_plots",
           showWarnings = FALSE)
dir.create("./analysis/processing_tables",
           showWarnings = FALSE)
dir.create("./analysis/analysis_tables",
           showWarnings = FALSE)

#EXTRACTION STANDARDS AND BLANKS ONLY----

stds <- filterFile(data_ms1,
                   grep("Extraction|SolventBlank",
                        data_ms1$name)) 

#4: Peak picking ---------------------------
#MatchedFilter (testing)----
#define parameters
# mfp <- MatchedFilterParam()
# mfp@binSize <- 0.005
# mfp@max <- 15
# mfp@snthresh <- 6
# mfp@steps <- 2
# mfp@mzdiff <- 0.6
# mfp@impute <- "linbase"
# mfp@distance <- 3
# mfp@sigma <- 10
# 
# 
# #test peak finding:
# 
# #chromatogram test 1: disaccharide
# chr_raw1 <- chromatogram(data_ms1,
#                          rt = c(1000, 1300),
#                          mz = c(562, 563))
# chr_raw1_mfp <- findChromPeaks(chr_raw1,
#                                param = mfp)
# plot(chr_raw1_mfp)
# 
# 
# #chromatogram test 2: trisaccharide
# chr_raw2 <- chromatogram(data_ms1,
#                          rt = c(900, 1200),
#                          mz = c(724, 725))
# chr_raw2_mfp <- findChromPeaks(chr_raw2,
#                                param = mfp)
# plot(chr_raw2_mfp)
# 
# chr_raw2_cwp <- findChromPeaks(chr_raw2,
#                                param = cwp)
# plot(chr_raw2_cwp)
# 
# #chromatogram test 3: tetrasaccharide
# chr_raw3 <- chromatogram(data_ms1,
#                          rt = c(1000, 1600),
#                          mz = c(886, 887))
# chr_raw3_mfp <- findChromPeaks(chr_raw3,
#                                param = mfp)
# plot(chr_raw3_mfp)
# 
# chr_raw3_cwp <- findChromPeaks(chr_raw3,
#                                param = cwp)
# plot(chr_raw3_cwp)
# 
# rm(chr_raw1,
#    chr_raw1_mfp,
#    chr_raw2,
#    chr_raw2_cwp,
#    chr_raw2_mfp,
#    chr_raw3,
#    chr_raw3_cwp,
#    chr_raw3_mfp)

#CentWave (use for now) ----
cwp<-CentWaveParam()
cwp@ppm<-1.6
cwp@peakwidth<-c(10,50)
cwp@snthresh<-3

chr_raw1 <- chromatogram(stds,
                         rt = c(1000, 1300),
                         mz = c(562, 563))
chr_raw1_cwp <- findChromPeaks(chr_raw1,
                               param = cwp)
plot(chr_raw1_cwp)


chr_raw2 <- chromatogram(stds,
                     mz = c(724.3, 724.4),
                     rt = c(1000, 1200))
chr_raw2_cwp <- findChromPeaks(chr_raw2,
                               param = cwp)
plot(chr_raw2_cwp)


stds_pks <-findChromPeaks(stds, 
                          param=cwp)


#Refine peaks by merging neighbouring peaks-----
mpp <- MergeNeighboringPeaksParam()
mpp@expandRt <- 2
mpp@expandMz <- 0
mpp@ppm <- 2

chr_raw2_cwp_mpp <- refineChromPeaks(chr_raw2_cwp,
                                     param = mpp)

plot(chr_raw2_cwp_mpp)


stds_pks_mp <-refineChromPeaks(stds_pks, 
                               param=mpp)


#5: Group peaks to create "features"---------
#parameters
pdp <- PeakDensityParam(sampleGroups = stds_pks$sample_type,
                        binSize = 0.005,
                        bw = 4,
                        minSamples = 1,
                        minFraction = 0.25) 

#dry run
chr1 <- chromatogram(stds_pks_mp,
                     mz = c(724.3, 724.4),
                     rt = c(1000, 1200))
plotChromPeakDensity(chr1, 
                     param = pdp,
                     peakPch = 16)

chr2 <- chromatogram(stds_pks_mp,
                     mz = c(562.28, 562.31),
                     rt = c(750, 950))
plotChromPeakDensity(chr2, 
                     param = pdp,
                     peakPch = 16)

chr3 <- chromatogram(stds_pks_mp,
                     mz = c(887.35, 887.45),
                     rt = c(1110, 1230))
plotChromPeakDensity(chr3, 
                     param = pdp,
                     peakPch = 16)

stds_pks_mp_grp <-groupChromPeaks(stds_pks_mp, 
                                  param=pdp)

#6: Fill in missing peaks----------
fpp <- FillChromPeaksParam()
stds_pks_mp_grp_fld <- fillChromPeaks(stds_pks_mp_grp)

#7: Set object to xcmsSet -----
xset <- as(stds_pks_mp_grp_fld, "xcmsSet")
sampnames(xset) <- pData(stds_pks_mp_grp_fld)$name
sampclass(xset) <- pData(stds_pks_mp_grp_fld)$sample_type

#8. Isotope/adduct detection ----
#NOTES
#all operations (grouping, isotope detection, adduct detection etc) are done on
#the FEATURES from xcms correspondence analysis - not the detected chromatographic
#peaks!!! 

##create xsannotate object
an <- xsAnnotate(xset)

##Group peaks 
an <- groupFWHM(an)

##Annotate isotope peaks
an <- findIsotopes(an, 
                   mzabs=0.04,
                   minfrac = 0.25)

##Peak grouping after correlation information 
an <- groupCorr(an, 
                  cor_eic_th=0.75)

##Find adducts
an <- findAdducts(an, 
                    polarity="positive")

#9. Peak list filtering ----
##get peak list
pl <-getPeaklist(an)

##make rownames from rt and mz of features
rownames(pl)<-paste(round(pl$rt,1),
                    round(pl$mz,3),
                    sep="_")

#blank exclusion (detected peaks (not filled))
pl_be <- pl[pl$SolventBlank == 0,]

##filter for isotopes or adducts
pl_be_isoadd <- pl_be[pl_be$isotopes!=""| 
                          pl_be$adduct!="",]

#set NA to be 0
pl_be_isoadd[is.na(pl_be_isoadd)] <- 0

#blank exclusion (peak intensity, includes integration on filled peaks)
sample_peaks <- pl_be_isoadd %>% 
    filter_at(vars(contains("SolventBlank")),
                   all_vars(.<1000)) 

sampleColNames <-  pd$name[pd$sample_type!="SolventBlank" & 
                               pd$sample_type!= "ExtractionBlank"]

sample_peaks <- sample_peaks %>% 
    filter_at(vars(any_of(sampleColNames)),
              any_vars(.>2e4))

sample_peaks <- cbind(rt_round = round_any(sample_peaks$rt, 
                                           5),
                      sample_peaks)

#10. Collapse features with multiple isotopes ----
##collapse features with multiple isotopes
setDT(sample_peaks)
#split out features without an isotope detected
sample_peaks_noiso <- sample_peaks[sample_peaks$isotopes=="",]
sample_peaks_iso <- sample_peaks[!sample_peaks$isotopes=="",]
#make column for the isotope group
sample_peaks_iso$isotope_group <- sample_peaks_iso$isotopes %>% 
    sub("\\[M.*", "", .)
#order isotopes within each group correctly
sample_peaks_iso$isotope_number <- sample_peaks_iso$isotopes %>% 
    sub(".*\\[M\\].*", "0", .) %>% 
    sub(".*\\[M\\+", "", .) %>% 
    sub("\\].*", "", .) %>% 
    as.numeric()
sample_peaks_iso <- sample_peaks_iso[order(isotope_group, 
                                           isotope_number),]
#get concatenated list of isotopes per group
iso_concat <- sample_peaks_iso[, 
                               list(isotopes = paste(isotopes, 
                                                     collapse = ', ')), 
                               by = isotope_group]
#remove duplicates within each isotope group (will keep [M] isotope)
#because of ordering
sample_peaks_iso <- unique(sample_peaks_iso, 
                           by = "isotope_group")
#merge to get concatenated isotope lists
sample_peaks_iso <- merge(sample_peaks_iso,
                          iso_concat,
                          by = "isotope_group")
#clean up df
sample_peaks_iso <- sample_peaks_iso %>% 
    select(-c("isotope_group",
              "isotope_number",
              "isotopes.x"))
names(sample_peaks_iso)[names(sample_peaks_iso) == 'isotopes.y'] <- 'isotopes'

#replace features that don't contain [M] isotope with [M] isotope
temp <- sample_peaks_iso %>% 
    filter(!grepl("\\[M\\]", isotopes))

sample_peaks_iso <- sample_peaks_iso %>% 
    filter(grepl("\\[M\\]", isotopes))

temp.vec <- temp$isotopes %>% 
    sub("\\[M.*", "", .)

pl_be_isoadd$isotope_group <- pl_be_isoadd$isotopes %>% 
    sub("\\[M.*", "", .)
pl_be_isoadd$isotope_number <- pl_be_isoadd$isotopes %>% 
    sub(".*\\[M\\].*", "0", .) %>% 
    sub(".*\\[M\\+", "", .) %>% 
    sub("\\].*", "", .) %>% 
    as.numeric()

temp <- pl_be_isoadd[pl_be_isoadd$isotope_group %in% temp.vec,] %>% 
    filter(isotope_number == 0)

temp$isotope_group <- NULL
temp$isotope_number <- NULL

temp <- cbind(rt_round = round_any(temp$rt, 
                                   5),
              temp)

sample_peaks_iso <- rbind(sample_peaks_iso, 
                          temp)

rm(temp,
   temp.vec)

#merge features with and without isotopes
sample_peaks <- rbind.fill(sample_peaks_noiso,
                           sample_peaks_iso)

setDT(sample_peaks)
sample_peaks <- sample_peaks[order(rt_round,
                                   mz),]

fwrite(sample_peaks, 
       file = "./analysis/analysis_tables/sample_peaks.txt",
       sep = "\t")

#11. Annotate features with predicted ions ----
#import table of predicted ions
mz_predicted <- fread("dp1to8-procainamide-allmod-posmode.txt")

#remove "extra" columns
extraCol <- c('index',
              'formula', #formulas are wrong in current version of script!
              'mass',
              'isomers' #also not confident in this
              )

mz_predicted <- mz_predicted %>% 
    select(-all_of(extraCol))
#make long format
mz_predicted_wide <- mz_predicted
mz_predicted <- gather(mz_predicted_wide, 
                       key = "ion", 
                       value = "mz", 
                       -dp,
                       -name)

#remove ions with m/z value of NA (i.e. ions with m/z values outside scan range)
mz_predicted <- na.omit(mz_predicted)


#add glucosamine 
pred_glucosamine_procA <- data.frame(dp = 1,
                                     name = "hex-1-amine-1-procA",
                                     ion = "[M+H]",
                                     mz = 399.260197)
mz_predicted <- rbind(mz_predicted,
                      pred_glucosamine_procA)


#make data.table
setDT(mz_predicted)
setDT(sample_peaks)
#create interval to overlap with (same width as for peak grouping)
mz_predicted$mz <- as.numeric(mz_predicted$mz)
mz_predicted$mzmin <- mz_predicted$mz-0.0025
mz_predicted$mzmax <- mz_predicted$mz+0.0025
#match using foverlaps from data.table (very fast)
setkey(mz_predicted, mzmin, mzmax)
sample_peaks_annot <- foverlaps(sample_peaks,
                                mz_predicted)

#change NA values created during matching (features with no match) to be blank
sample_peaks_annot <- sample_peaks_annot %>% 
    replace_na(list("dp"="",
                    "name" = "", 
                    "ion"= "", 
                    "mz" = "", 
                    "mzmin" = "",
                    "mzmax"= ""))

sample_peaks_matched <- sample_peaks_annot[sample_peaks_annot$dp!="",]
sample_peaks_unmatched <- sample_peaks_annot[sample_peaks_annot$dp=="",]



#aggregate so that if there are multiple predictions for one feature
#they are shown in the same row. delete all of the other extra columns added 
#during matching

theoretical <- paste0(sample_peaks_matched$name, 
                      ":", 
                      sample_peaks_matched$ion)

sample_peaks_matched <- cbind(theoretical,
                              sample_peaks_matched)

names <- setdiff(names(sample_peaks_matched), 
                 names(mz_predicted))
names <- names[!names == "theoretical"]
setDT(sample_peaks_matched)
sample_peaks_matched <- sample_peaks_matched[, 
                       list(theoretical = paste(theoretical, 
                                                collapse = ', ')), 
                       by = names]
names(sample_peaks_matched)[names(sample_peaks_matched) == 'i.mz'] <- 'mz'
names(sample_peaks_matched)[names(sample_peaks_matched) == 'i.mzmin'] <- 'mzmin'
names(sample_peaks_matched)[names(sample_peaks_matched) == 'i.mzmax'] <- 'mzmax'

sample_peaks_matched_old <- sample_peaks_matched

colOrder <- names(sample_peaks_matched)
colOrder <- colOrder[1:length(colOrder)-1]
colOrder <- c("theoretical", colOrder)

setDF(sample_peaks_matched)
sample_peaks_matched <- sample_peaks_matched[,colOrder]


fwrite(sample_peaks_matched, 
       file = "./analysis/analysis_tables/sample_peaks_matched.txt",
       sep = "\t")

#12: heatmap and volcano plot-------
dir.create("./analysis/heatmap",
           showWarnings = FALSE)
#all sample peaks ----

sample_peaks_unmatched$theoretical <- paste0("m/z_",
                                             round(sample_peaks_unmatched$i.mz,
                                                   3),
                                             "_unknown")

sample_peaks_all <- rbind.fill(sample_peaks_unmatched,
                               sample_peaks_matched)

#subset to only have intensity counts
mask <- pd$name %>% 
    gsub("\\+|\\s", ".", .)
mask <- mask[grep("Extraction", mask)]

setDF(sample_peaks_all)
rownames(sample_peaks_all) <- paste0("rt",
                                 round(sample_peaks_all$rt, 1),
                                 "_",
                                 sample_peaks_all$theoretical)


counts <- sample_peaks_all %>% 
    select(mask)

#set factor level 
group<-factor(mask %>% sub("col[123]_", "", .))

#DGEList:Creates a DGEList object from a table of counts 
#(rows=features, columns=samples), 
#group indicator for each column, 
#library size (optional) and a table of feature annotation (optional).

y_n <- DGEList(counts=counts,
               group=group)
#filterByExpr {edgeR}
#determine which features have sufficiently large counts to be retained for stats
#output is a logical vector
keep_n <- filterByExpr(y_n)
y_n <- y_n[keep_n,,keep.lib.sizes=FALSE]

#Calculate normalisation factors to scale the raw library sizes
y_n <- calcNormFactors(y_n)

#creates a design (or model) matrix, e.g., by expanding factors to a 
#set of dummy variables (depending on the contrasts) and 
#expanding interactions similarly.
design <- model.matrix(~group)

#estimate disparity
y_n <- estimateDisp(y_n,design)
y_n <- estimateCommonDisp(y_n)

#test difference
#output:
#log2-fold-change (logFC), 
#the average log2-counts-per-million (logCPM), 
#and the two-sided p-value (PValue).
tested_n <-exactTest(y_n)
cairo_pdf("./analysis/heatmap/pvalue_hist_sample-peaks-all.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
hist(tested_n$table[,"PValue"], breaks=50)
dev.off()

#extract most different
result_n <- topTags(tested_n, 
                    n=nrow(tested_n$table))

#set up data for volcano plot
volcanoData_n <- cbind(result_n$table$logFC,
                       result_n$table$FDR,
                       -log10(result_n$table$FDR))

volcanoData_n <- as.data.frame(volcanoData_n)
rownames(volcanoData_n) <- rownames(result_n[["table"]])
colnames(volcanoData_n) <- c("logFC", 
                             "FDR",
                             "negLogFDR")


volcanoData_n$diff <- NA
volcanoData_n$diff[volcanoData_n$logFC < -2 &
                       volcanoData_n$FDR < 0.05] <- "DOWN"
volcanoData_n$diff[volcanoData_n$logFC > 2 &
                       volcanoData_n$FDR < 0.05] <- "UP"

volcanoData_n$diff[is.na(volcanoData_n$diff)] <- "stable"

volcanoData_n$label <- rownames(volcanoData_n)
volcanoData_n$label[volcanoData_n$diff=="stable"] <- ""
volcanoData_n$label <- volcanoData_n$label %>% 
    sub(".*unknown.*", "unknown", .) %>% 
    sub(".*procA.*", "matched", .)
volcanoData_n$label_diff <- paste0(volcanoData_n$label,
                                   volcanoData_n$diff)

volcanoData_n$label_diff <- factor(volcanoData_n$label_diff,
                                   levels = c("unknownUP",
                                              "matchedUP",
                                              "stable",
                                              "unknownDOWN",
                                              "matchedDOWN"))

pal <- c("#98B56A",
         "#086D70",
         "grey",
         "#DA95CD",
         "#8D0141")
names(pal) <- c("unknownDOWN",
                "matchedDOWN",
                "stable",
                "unknownUP",
                "matchedUP")

volcanoData_n$annot <- rownames(volcanoData_n)



# volcanoData_n$annot <- volcanoData_n$annot %>% 
#     sub(".*unknown", NA, .)
volcanoData_n$annot[volcanoData_n$label_diff == "stable"] <- NA

tiff("./analysis/heatmap/volcano_plot_sample-peaks-all_nolabels.tiff",
     res = 300,
     units = "in",
     width = 12,
     height = 4)

ggplot(data = volcanoData_n, 
       aes(x = logFC, 
           y = negLogFDR, 
           fill=label_diff,
           group = label_diff,
           label = annot
       )) +
    geom_point(alpha=1, 
               size=4,
               shape = 21) +
    # geom_label(size = 2,
    #            family = "Avenir",
    #            colour = "white",
    #            #nudge_y = 1,
    #            #nudge_x = 5
    #            ) +
    #geom_text(aes(label = label)) +
    scale_fill_manual(values=pal) +
    #xlim(c(-4.5, 4.5)) +
    geom_vline(xintercept=c(-2,2),lty=2,col="black",lwd=0.5) +
    geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.5) +
    labs(x="log2(fold change)",
         y="-log10(false discovery rate)")  +
    theme_classic() +
    theme(text = element_text(family = "Avenir", size = 16),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 12),
          legend.position = "none") 

dev.off()

tiff("./analysis/heatmap/volcano_plot_sample-peaks-all_labels.tiff",
     res = 300,
     units = "in",
     width = 12,
     height = 4)

ggplot(data = volcanoData_n, 
       aes(x = logFC, 
           y = negLogFDR, 
           fill=label_diff,
           group = label_diff,
           label = annot
       )) +
    geom_point(alpha=1, 
               size=4,
               shape = 21) +
    geom_label(size = 2,
               family = "Avenir",
               colour = "white",
               #nudge_y = 1,
               #nudge_x = 5
               ) +
    #geom_text(aes(label = label)) +
    scale_fill_manual(values=pal) +
    #xlim(c(-4.5, 4.5)) +
    geom_vline(xintercept=c(-2,2),lty=2,col="black",lwd=0.5) +
    geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.5) +
    labs(x="log2(fold change)",
         y="-log10(false discovery rate)")  +
    theme_classic() +
    theme(text = element_text(family = "Avenir", size = 16),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 12),
          legend.position = "none") 

dev.off()



#counts per million and log2 (normalise and transform)
# dge_n <- cpm(y_n, 
#              log=TRUE, 
#              prior.count = 1)

dge_n <- log2(y_n$counts + 1)

#subset
selY_n <- dge_n[rownames(result_n$table)[result_n$table$FDR<0.05 & 
                                             result_n$table$logFC < -2 |
                                             result_n$table$FDR<0.05 & 
                                             result_n$table$logFC > 2],]



#make heatmap
cimColour <- viridis(1000)
cimColurCols <- c(rep("#DA95CD", 3),
                  rep("#ADA4E2", 3))


pal <- c("#98B56A",
         "#086D70",
         "grey",
         "#DA95CD",
         "#8D0141")

cimColourRows <- rownames(volcanoData_n)[volcanoData_n$diff!="stable"]
cimColourRows[cimColourRows %in% 
                  rownames(volcanoData_n)[
                      volcanoData_n$label_diff=="unknownDOWN"]] <- "#98B56A"
cimColourRows[cimColourRows %in% 
                  rownames(volcanoData_n)[
                      volcanoData_n$label_diff=="matchedDOWN"]] <- "#086D70"
cimColourRows[cimColourRows %in% 
                  rownames(volcanoData_n)[
                      volcanoData_n$label_diff=="unknownUP"]] <- "#DA95CD"
cimColourRows[cimColourRows %in% 
                  rownames(volcanoData_n)[
                      volcanoData_n$label_diff=="matchedUP"]] <- "#8D0141"    


tiff("./analysis/heatmap/heatmap_matched_unmatched.tiff",
     units = "in",
     res = 300,
     width = 12,
     height = 16)

# svg("./analysis/heatmap/heatmap_matched_unmatched.svg",
#     width = 12,
#     height = 9)

par(mar= c(10, 15, 15, 40))
cim(selY_n, 
    color = cimColour,
    symkey=FALSE,
    #mar=c(5,20),
    row.sideColors = cimColourRows,
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

tiff("./analysis/heatmap/heatmap_matched_unmatched_3.tiff",
     units = "in",
     res = 300,
     width = 12,
     height = 16)

# svg("./analysis/heatmap/heatmap_matched_unmatched_3.svg",
#     width = 12,
#     height = 16)


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
          key.xlab = expression("log"[2]*"(integrated intensity)"),
          key.ylab = NA,

          #layout param
          margins=c(10,58), # height, width margins around plot
          lwid= c(0.08,0.005, 0.3),
          lhei = c(1.3,15),
          lmat = rbind(c(5,5,4), c(3,1,2)),
          cexRow = 1.5,# text size rows
          cexCol = 2, # text size cols
          srtCol = 45
          
)

dev.off()



#sample peaks annotated----

#subset to only have intensity counts
mask2 <- pd$name %>% 
    gsub("\\+|\\s", ".", .)
mask2 <- mask2[grep("Extraction", mask2)]

rownames(sample_peaks_matched) <- paste0("rt",
                                         round(sample_peaks_matched$rt, 1),
                                         "_",
                                         sample_peaks_matched$theoretical)


counts2 <- sample_peaks_matched %>% 
    select(mask2)

#set factor level 
group2<-factor(mask2 %>% sub("col[123]_", "", .))

#DGEList:Creates a DGEList object from a table of counts 
#(rows=features, columns=samples), 
#group indicator for each column, 
#library size (optional) and a table of feature annotation (optional).

y_n2 <- DGEList(counts=counts2,
               group=group2)
#filterByExpr {edgeR}
#determine which features have sufficiently large counts to be retained for stats
#output is a logical vectir
keep_n2 <- filterByExpr(y_n2)
y_n2 <- y_n2[keep_n2,,keep.lib.sizes=FALSE]

#Calculate normalisation factors to scale the raw library sizes
y_n2 <- calcNormFactors(y_n2)

#creates a design (or model) matrix, e.g., by expanding factors to a 
#set of dummy variables (depending on the contrasts) and 
#expanding interactions similarly.
design2 <- model.matrix(~group2)

#estimate disparity
y_n2 <- estimateDisp(y_n2,design2)
y_n2 <-estimateCommonDisp(y_n2)

#test difference
#output:
#log2-fold-change (logFC), 
#the average log2-counts-per-million (logCPM), 
#and the two-sided p-value (PValue).
tested_n2 <-exactTest(y_n2)
cairo_pdf("./analysis/heatmap/pvalue_hist_sample-peaks-matched.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
hist(tested_n2$table[,"PValue"], breaks=50)
dev.off()

#extract most different
result_n2 <- topTags(tested_n2, 
                    n=nrow(tested_n2$table))

#set up data for volcano plot
volcanoData_n2 <- cbind(result_n2$table$logFC, 
                        result_n2$table$FDR,
                       -log10(result_n2$table$FDR))


volcanoData_n2 <- as.data.frame(volcanoData_n2)
rownames(volcanoData_n2) <- rownames(result_n2[["table"]])
colnames(volcanoData_n2) <- c("logFC", "FDE", "negLogFDR")

volcanoData_n2$diff <- NA
volcanoData_n2$diff[volcanoData_n2$logFC < -2 &
                        volcanoData_n2$negLogFDR > -log10(0.05)] <- "DOWN"
volcanoData_n2$diff[volcanoData_n2$logFC > 2 &
                       volcanoData_n2$negLogFDR > -log10(0.05)] <- "UP"

volcanoData_n2$diff[is.na(volcanoData_n2$diff)] <- "stable"

volcanoData_n2$label <- rownames(volcanoData_n2)
volcanoData_n2$label[volcanoData_n2$diff == "stable"] <- NA

tiff("./analysis/heatmap/volcano_plot_sample-peaks-matched_nolabels.tiff",
     res = 300,
     units = "in",
     width = 12,
     height = 9)

ggplot(data = volcanoData_n2, 
       aes(x = logFC, 
           y = negLogFDR, 
           fill=diff,
           )) +
    geom_point(alpha=1, 
               size=2.5,
               shape = 21) +
    #geom_text(aes(label = label)) +
    scale_fill_manual(values=c("grey", "#DA95CD")) +
    #xlim(c(-4.5, 4.5)) +
    geom_vline(xintercept=c(-2,2),lty=2,col="black",lwd=0.5) +
    geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.5) +
    labs(x="log2(fold change)",
         y="-log10(false discovery rate)")  +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 12),
          legend.position = "none") 

dev.off()

tiff("./analysis/heatmap/volcano_plot_sample-peaks-matched_nolabels.tiff",
     res = 300,
     units = "in",
     width = 12,
     height = 9)

ggplot(data = volcanoData_n2, 
       aes(x = logFC, 
           y = negLogFDR, 
           fill=diff,
           label = label
       )) +
    geom_point(alpha=1, 
               size=2.5,
               shape = 21) +
    geom_label(size = 2,
               family = "Avenir",
               colour = "white",
               #nudge_y = 1,
               #nudge_x = 5
    ) +
    scale_fill_manual(values=c("grey", "#DA95CD")) +
    #xlim(c(-4.5, 4.5)) +
    geom_vline(xintercept=c(-2,2),lty=2,col="black",lwd=0.5) +
    geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.5) +
    labs(x="log2(fold change)",
         y="-log10(false discovery rate)")  +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 12),
          legend.position = "none") 

dev.off()


chr_tri <- chromatogram(stds_pks_mp_grp_fld,
                        mz = c(724.33, 724.36),
                        rt = c(1050, 1200))
chr_tri.ex <- chr_tri[,chr_tri$name %in% mask]
par(mfrow=c(3,2))
for (i in 1:6){
    plot(chr_tri.ex[,i],
         main = mask[i])
}


#counts per million and log2 (normalise and transform)
# dge_n <- cpm(y_n, 
#              log=TRUE, 
#              prior.count = 1)

dge_n <- log2(y_n$counts + 1)

#subset
selY_n <- dge_n[rownames(result_n$table)[result_n$table$FDR<0.05 & 
                                             result_n$table$logFC < -2 |
                                             result_n$table$FDR<0.05 & 
                                             result_n$table$logFC > 2],]


#make heatmap
cimColour <- rev(viridis(1000))
cimColurCols <- c(rep("#DA95CD", 3),
                  rep("#ADA4E2", 3))

cairo_pdf("./analysis/heatmap/heatmap_matched_v2.pdf",
          family = "Avenir",
          width = 12,
          height = 9)

svg("./analysis/heatmap/heatmap_matched_v2.svg",
          width = 12,
          height = 9)
cim(selY_n, 
    color = cimColour,
    symkey=FALSE,
    mar=c(12,40),
    col.sideColors = cimColurCols
)
dev.off()

#sample peaks annotated and extraction blanks----

#subset to only have intensity counts
mask <- pd$name %>% 
    gsub("\\+|\\s", ".", .)
mask <- mask[grep("Poseidon|ExtractionBlank", mask)]

rownames(sample_peaks_matched) <- paste0("rt",
                                         round(sample_peaks_matched$rt, 1),
                                         "_",
                                         sample_peaks_matched$theoretical)


counts <- sample_peaks_matched %>% 
    select(mask)

#set factor level 
group<-factor(mask 
              %>% sub("_\\d{2}h[abc]$", "", .) %>% 
                  sub("col[123]_", "", .))

#DGEList:Creates a DGEList object from a table of counts 
#(rows=features, columns=samples), 
#group indicator for each column, 
#library size (optional) and a table of feature annotation (optional).

y_n <- DGEList(counts=counts,
               group=group)
#filterByExpr {edgeR}
#determine which features have sufficiently large counts to be retained for stats
#output is a logical vectir
keep_n <- filterByExpr(y_n)
y_n <- y_n[keep_n,,keep.lib.sizes=FALSE]

#Calculate normalisation factors to scale the raw library sizes
y_n <- calcNormFactors(y_n)

#creates a design (or model) matrix, e.g., by expanding factors to a 
#set of dummy variables (depending on the contrasts) and 
#expanding interactions similarly.
design <- model.matrix(~group)

#estimate disparity
y_n <- estimateDisp(y_n,design)
y_n<-estimateCommonDisp(y_n)

#test difference
#output:
#log2-fold-change (logFC), 
#the average log2-counts-per-million (logCPM), 
#and the two-sided p-value (PValue).
tested_n <-exactTest(y_n)
cairo_pdf("./analysis/heatmap/pvalue_hist_matched_exblk.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
hist(tested_n$table[,"PValue"], breaks=50)
dev.off()

#extract most different
result_n <- topTags(tested_n, 
                    n=nrow(tested_n$table))

#plot volcano plot
volcanoData_n <- cbind(result_n$table$logFC, 
                       -log10(result_n$table$FDR))
colnames(volcanoData_n) <- c("logFC", "negLogFDR")

cairo_pdf("./analysis/heatmap/volcano_plot_matched_exblk.pdf",
          family = "Avenir",
          width = 12,
          height = 9)

plot(volcanoData_n, pch=19)
abline(h=-log10(0.05), col="red") #add line for cutoff to include in plot
text(-17,-log10(0.001), "-log10(0.05)", col = "red")
dev.off()

#counts per million and log2 (normalise and transform)
dge_n <- cpm(y_n, 
             log=TRUE, 
             prior.count = 1)

#subset
selY_n <- dge_n[rownames(result_n$table)[result_n$table$FDR<=0.05],]


#make heatmap
cimColour <- rev(viridis(1000))

cairo_pdf("./analysis/heatmap/heatmap_matched_exblk.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
cim(t(selY_n), 
    color = cimColour,
    symkey=FALSE,
    mar=c(30,12),
    #row.sideColors = cimColurRows.all
)
dev.off()

#13: Extract and format eic -----
#get vectors
setDF(sample_peaks_matched)
sample_peaks_matched <- sample_peaks_matched[
    order(sample_peaks_matched$rt_round),]

ions.found.vector <- sample_peaks_matched$theoretical
mz.found.vector <- sample_peaks_matched$mz
rt.found.vector <- sample_peaks_matched$rt_round
rt.found.vector.min <- rt.found.vector - 30
rt.found.vector.max <- rt.found.vector + 30

#remove solvent blanks
stds_noblanks <- filterFile(stds,
                    file = which(grepl("Extraction",
                                       stds$name)))

#get phenodata vectors
stds.names <- stds_noblanks$name
stds.groups <- stds_noblanks$sample_type

#extract chromatograms
#each peak is extracted with a retention time range of 1 min around peak
chr_list <- list()

error = 0.001

for (i in 1:length(mz.found.vector)){
    mzr = c(mz.found.vector[i] - error,
            mz.found.vector[i] + error)
    rtr = c(rt.found.vector.min[i],
            rt.found.vector.max[i])
    chr_list[[i]] <- chromatogram(stds_noblanks, 
                                  mz = mzr,
                                  rt = rtr)
}

#extract intensity and rt values
chr_int_list <- list()
for (i in 1:length(stds.names)){
    chr_int_list[[i]] <- lapply(chr_list, function(x) {
        x[[i]]@intensity
    }) 
}

chr_rt_list <- list()
for (i in 1:length(stds.names)){
    chr_rt_list[[i]] <- lapply(chr_list, function(x) {
        x[[i]]@rtime
    }) 
}

#build data frame (long format)
eic.df <- data.frame(ion = as.character(),
                     sample = as.character(),
                     group = as.character(),
                     rt = as.numeric(),
                     feature_number = as.numeric(),
                     intensity = as.numeric())

for (i in 1:length(stds.names)){
    for (j in 1:length(mz.found.vector)){
        rt = chr_rt_list[[i]][[j]]/60
        intensity = chr_int_list[[i]][[j]]
        sample = rep(stds.names[i], length(rt))
        group = rep(stds.groups[i], length(rt))
        ion = rep(ions.found.vector[j], length(rt))
        feature_number = rep(j, length(rt))
        temp <- data.frame(ion = ion,
                           sample = sample,
                           group = group,
                           rt = rt,
                           feature_number = feature_number,
                           intensity = intensity)
        eic.df <- rbind(eic.df,
                        temp)
    }
}

eic.df[is.na(eic.df)] <- 0

#set variables as factors
eic.df$group <- factor(eic.df$group,
                       levels = c("ExtractionStandard",
                                  "ExtractionBlank"))

#create variable for line grouping
eic.df$feature_number_sample <- paste0(eic.df$feature_number, 
                                       "_",
                                       eic.df$sample)

#artificially add vertical space to each chromatogram
# eic.df2 <- eic.df[eic.df$ion == ions.found.vector[1],]
# 
# i <- 2
# while (i < length(ions.found.vector)+1){
#     temp <- eic.df[eic.df$ion == ions.found.vector[i],]
#     temp$intensity <- temp$intensity + ((i-1)*1e5)
#     eic.df2 <- rbind(eic.df2,
#                      temp)
#     
#     rm(temp)
#     i = i+1
# }

#make dataframe for label at each peak apex
#note that the ions are seen exactly -2H below main saccharides at same retention tim
#do not label these - believe that they are an artefact of ionisation...

#hex-6-procA:[M+H] = mz 1210
#mz 1208 = "hex-5-pent-1-carboxyl-1-omethyl-1-procA:[M+H], hex-6-carboxyl-1-deoxy-1-procA:[M+H], hex-4-pent-1-phosphate-1-nacetyl-2-omethyl-2-procA:[M+H], hex-5-phosphate-1-deoxy-1-nacetyl-2-omethyl-1-procA:[M+H]" 

#hex-5-procA:[M+H] = mz 1048
#mz 1046 = "hex-4-pent-1-carboxyl-1-omethyl-1-procA:[M+H], hex-5-carboxyl-1-deoxy-1-procA:[M+H]"                                                                                        

#hex-3-procA:[M+H] = mz 724
#mz 722 = "hex-2-pent-1-carboxyl-1-omethyl-1-procA:[M+H], hex-3-carboxyl-1-deoxy-1-procA:[M+H]"    

#hex-2-procA:[M+H] = mz 562
#mz 560  = "hex-1-pent-1-carboxyl-1-omethyl-1-procA:[M+H], hex-2-carboxyl-1-deoxy-1-procA:[M+H]"       

eic.apex.labels <- eic.df[order(-eic.df$intensity),]

setDT(eic.apex.labels)
eic.apex.labels <- eic.apex.labels %>% 
    unique(.,
           by = c("feature_number", "group"))

eic.apex.labels <- eic.apex.labels[eic.apex.labels$intensity > 5e4,]
eic.apex.labels$label <- eic.apex.labels$ion %>% 
    sub("pent-1-procA:\\[M\\+H\\]", "xylose", .) %>% 
    sub("pent-1-omethyl-1.*", "fucose", .) %>% 
    sub("hex-1-procA:\\[M\\+H\\]", "glucose", .) %>% 
    sub("hex-2-procA:\\[M\\+H\\]", "laminaribiose", .) %>% 
    sub("hex-3-procA:\\[M\\+H\\]", "mannotriose", .) %>% 
    sub("hex-4-procA:\\[M\\+H\\]", "cellotetraose", .) %>% 
    sub("hex-5-procA:\\[M\\+H\\]", "pentasaccharide", .) %>% 
    sub("hex-6-procA:\\[M\\+H\\]", "laminarihexaose", .) %>% 
    sub("hex-5-pent-1-carboxyl-1-omethyl-1.*|hex-4-pent-1-carboxyl-1-omethyl-1.*", NA, .) %>% 
    sub("hex-2-pent-1-carboxyl-1-omethyl-1.*", NA, .) %>% 
    sub(".*-procA.*", "", .) 

#set extraction blank 3 as its own thing

eic.df$var <- ""
eic.df$var[eic.df$sample == "col1_ExtractionBlank" |
               eic.df$sample == "col2_ExtractionBlank"] <- "Blanks"

eic.df$var[eic.df$sample == "col3_ExtractionBlank" ] <- "Carryover"
eic.df$var[eic.df$group == "ExtractionStandard" ] <- "Standards"

eic.df$var <- factor(eic.df$var,
                     levels = c("Blanks",
                                "Standards",
                                "Carryover"))


#14: Plot EIC ----

# svg("./analysis/analysis_plots/standards-eic-v1.svg",
#     height = 6,
#     width = 12)
# cairo_pdf("./analysis/analysis_plots/standards-eic-v1.pdf",
#     height = 6,
#     width = 12)
# ggplot() +
#     geom_line(mapping = aes(rt,
#                             intensity,
#                             group = feature_number_sample,
#                             #colour = ion
#                             ),
#               data = eic.df,
#               lwd = 1.2,
#               colour = "black"
#               ) +
#     geom_text(data = eic.apex.labels,
#                 mapping = aes(rt,
#                               intensity + 5e5,
#                               label = label),
#               family = "Avenir") +
#     facet_grid(rows = vars(group)) +
#     labs(x= "Retention time (min)",
#          y = "Intensity (a.u.)") +
#     scale_x_continuous(breaks = seq(5, 23, 1),
#                        limits = c(5,23),
#                        expand = c(0,0)) +
#     scale_y_continuous(expand = expansion(mult = c(0.02, 0.1)),
#                        labels = scales::scientific) +
#     theme_classic() +
#     theme(strip.background = element_blank(),
#           strip.text.y = element_blank(),
#           strip.text.x =element_blank(),
#           panel.border = element_rect(colour = "#848587",
#                                       size = 0.5,
#                                       fill = NA),
#           axis.line = element_blank(),
#           axis.text = element_text(size = 12, family = "Avenir"),
#           axis.title = element_text(size = 14, family = "Avenir LT 65 Medium"),
#           legend.position = "bottom") 
# dev.off()


#facet by ion, three columns
ions.for.facet <- c("hex-2-procA:[M+H]",
                    "hex-3-procA:[M+H]",
                    "hex-4-procA:[M+H]",
                    "hex-6-procA:[M+H]")

eic.df.ionfacet <- eic.df[eic.df$ion %in% ions.for.facet,]
eic.df.ionfacet$ion.label <- eic.df.ionfacet$ion %>% 
    sub("hex-2-procA.*", "DP2", .) %>% 
    sub("hex-3-procA.*", "DP3", .) %>% 
    sub("hex-4-procA.*", "DP4", .) %>% 
    sub("hex-6-procA.*", "DP6", .)


ggplot() +
    geom_line(mapping = aes(rt,
                            intensity,
                            group = feature_number_sample),
    data = eic.df.ionfacet,
    lwd = 0.7,
    colour = "black") +
    facet_grid(rows = vars(ion.label),
               cols = vars(var),
               scales = "free_y") +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    scale_x_continuous(breaks = seq(5, 25, 5),
                       limits = c(5,25),
                       expand = c(0,0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                       labels = scales::scientific) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.y = element_text(angle = 360),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          text = element_text(family = "Helvetica"),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 8),
          panel.spacing = unit(1, "lines"),
          legend.position = "none") 








