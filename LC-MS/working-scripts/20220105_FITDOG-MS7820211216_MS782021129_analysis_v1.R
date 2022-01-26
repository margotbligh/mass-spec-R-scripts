#analyse full MS run and MS/MS run data together
#running on linux server of MPI

setwd("/home/mbligh/ownCloud2/marglyco/Data/LC-MS_data/FITDOG/20211129_20211216_analysis")

#1: Install packages 
library(xcms)
library(ggplot2)
library(tidyverse)
library(scales)
library(data.table)
library(MSnbase)
library(CAMERA)
library(plyr)
library(viridis)
library(reticulate)

register(SerialParam())

scientific_function <- function(x) {
    text <- gsub("E0", "", gsub("e\\+0", "E", scales::scientific_format()(x)))
    text
}
#2. Import and inspect MS data --------------------------------------------------------
#get file paths to mzML files
fp1 <- dir(path = "../20211129/mzML", all.files = FALSE, full.names = TRUE)
fp2 <- dir(path = "../20211216/mzML", all.files = FALSE, full.names = TRUE)
fp <- c(fp1, fp2)

#create phenodata data.frame
#each sample must have a unique name!
pd <- data.frame(name = basename(fp) %>%
                     sub("MS78_", "", .) %>% 
                     sub("_\\d\\d.mzML", "", .),
                 sample_type = basename(fp) %>% 
                     sub(".*solvent.*", "solvent blank", .) %>% 
                     sub(".*MS78.*", "FITDOG digest", .),
                 substrate = basename(fp) %>% 
                     sub(".*bic.*|.*dig.*", "laminarin", .) %>% 
                     sub(".*hyp.*|.*ves.*", "fucoidan", .) %>% 
                     sub(".*twa.*", "sulphated mannan", .) %>% 
                     sub(".*yea.*", "mannan", .) %>% 
                     sub(".*MS78.*", NA, .),
                 source = basename(fp) %>%
                     sub(".*bic.*", "E. bicyclis", .) %>% 
                     sub(".*dig.*", "L. digitata", .) %>% 
                     sub(".*ves.*", "F. vesiculosus", .) %>%
                     sub(".*hyp.*", "L. hyperborea", .) %>% 
                     sub(".*twa.*", "T. weissflogii", .) %>% 
                     sub(".*yea.*", "S. cerevisisae", .) %>% 
                     sub(".*MS78.*", NA, .),
                 runtype = basename(fp) %>%
                     sub(".*20211129.*", "full MS", .) %>% 
                     sub(".*20211216.*", "full MS to DDA-MS2", .),
                 stringsAsFactors = FALSE)

pd$replicate <- NA
pd$replicate[pd$sample_type == "FITDOG digest"] <- pd$name[
    pd$sample_type == "FITDOG digest"] %>% 
    sub("^\\D{3}", "", .)

pd$source[is.na(pd$source)] <- "solvent blank"

#read in data
all_data <- readMSData(files = fp, pdata = new("NAnnotatedDataFrame", pd), 
                       mode = "onDisk")

#split MS1 and MS2
data <- all_data[all_data@featureData@data$msLevel == 1]
data_ms2 <- all_data

3: Create initial output directories -------------------------------------
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
save(data, file = "./analysis/RData/data.RData")
save(data_ms2, file = "./analysis/RData/data_ms2.RData")

#4: Plot TIC ----
dir.create("./analysis/processing_plots/tic", showWarnings = FALSE)

pal_group <- hcl.colors(n = length(unique(pd$source)),
                        palette = "Dark3")
names(pal_group) <- unique(pd$source)

#plot the tic as boxplot
tc <- split(tic(all_data), f = fromFile(all_data))
cairo_pdf("./analysis/processing_plots/tic/tic_boxplot.pdf",
          family = "Arial",width = 12, height = 9)
par(mar=c(9,5,1,1))
boxplot(tc, col = pal_group[all_data$source], ylab = "intensity", 
        main = "total ion current", names = all_data$name, las=2,cex.axis = 0.8)
dev.off()

#5: Peak picking (CentWave) ---------------------------
#set parameters
cwp<-CentWaveParam()
cwp@ppm<-5.5
cwp@peakwidth<-c(10,100)
cwp@snthresh<-5
cwp@noise <- 5000
cwp@prefilter <- c(3, 1000)

#pick peaks
data_peaks<-findChromPeaks(data, param=cwp)
data_ms2<-findChromPeaks(data_ms2, param=cwp)

#save RData object
save(data_peaks, file = "./analysis/RData/data_peaks.RData")
save(data_ms2, file = "./analysis/RData/data_ms2.RData")


#6: Group peaks to create "features"---------
#parameters
pdp <- PeakDensityParam(sampleGroups = data$source,
                        binSize = 0.005, bw = 6, minFraction = 1/3) 

#group peaks
data_peaks_grouped <- groupChromPeaks(data_peaks, param = pdp)

save(data_peaks_grouped, file = "./analysis/RData/data_peaks_grouped.RData")

#7: Fill in missing peaks----------
fpp <- FillChromPeaksParam()
data_peaks_grouped_filled <- fillChromPeaks(data_peaks_grouped)

save(data_peaks_grouped_filled, file = "./analysis/RData/data_peaks_grouped_filled.RData")

res <- data_peaks_grouped_filled

#8: Remove large data files from environment----------
rm(data, data_peaks, data_peaks_grouped, data_peaks_grouped_filled)    

#9: Save diffreport of xdata -----
xset <- as(res, "xcmsSet")
sampnames(xset) <- pData(res)$name
sampclass(xset) <- pData(res)$source

#10. Isotope picking----
##create xsannotate object
an <- xsAnnotate(xset)
##Group peaks of a xsAnnotate object according to their retention time 
an <- groupFWHM(an, perfwhm = 0.6)
##Annotate isotope peaks
an <- findIsotopes(an, mzabs=0.01)
##Peak grouping after correlation information into pseudospectrum groups 
an <- groupCorr(an, cor_eic_th=0.75)
##Find adducts
an <- findAdducts(an, polarity="negative")

#11. Peak list filtering and formatting----
#get peak list
pl <-getPeaklist(an)

#filter by blank exclusion (detected peaks)
pl_be <-pl[pl$solvent.blank==0,]

#add rounded retention time as first colum
pl_be <- cbind(rt_min = round(pl_be$rt/60, 1), pl_be)

#filter by retention time
pl_be_rt <- pl_be %>% filter(between(rt_min, 1, 22))

#make rownames from rt and mz of features
rownames(pl_be_rt)<-paste(round(pl_be_rt$rt,1), round(pl_be_rt$mz,3), sep="_")

#change NA to 0
pl_be_rt[is.na(pl_be_rt)] <- 0

#change name
peaks <- pl_be_rt

save.image("./analysis/RData/RData_20220126.RData")










