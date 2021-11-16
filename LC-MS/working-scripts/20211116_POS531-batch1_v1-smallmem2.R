setwd("/scratch/mbligh/POS531")

#1: Install packages 
library(xcms)
library(ggplot2)
library(tidyverse)
library(data.table)
library(CAMERA)

#2. Import and inspect MS data --------------------------------------------------------
fp <- dir(path = "mzML", all.files = FALSE, full.names = TRUE)

#create phenodata data.frame
pd <- data.frame(name = basename(fp) %>%gsub("MS31_20211020_|.mzML", "", .),
                 sampletype = basename(fp) %>% sub(".*POS531.*", "extract", .) %>% sub(".*procedure.*", "procedureblank", .) %>% sub(".*QCmix.*", "QCmix", .) %>% 
                     sub(".*Qcspike.*", "QCspikemix", .),
                 station = basename(fp) %>% gsub(".*POS531_|_Ex_.*", "", .) %>% sub("MS31_.*", "NA", .),
                 timepoint = basename(fp) %>% gsub(".*Ex_|[ABCDEF]_\\d{2}.mzML", "", .) %>% sub("MS31_.*", "NA", .),
                 group = basename(fp) %>% gsub("MS31_20211020_|[ABCDEF]_\\d{2}.mzML", "", .) %>% sub(".*procedure.*", "procedureblank", .) %>% 
                     sub(".*QCmix.*", "QCmix", .) %>% sub(".*Qcspike.*", "QCspikemix", .),
                 replicate = basename(fp) %>% gsub(".*Ex_\\d{1,2}|_\\d{2}.mzML", "", .) %>% sub(".*mix|.*spike-|.*blank", "", .),
                 run = basename(fp) %>%  gsub(".*Ex_\\d{1,2}[ABCDEF]_|.mzML", "",.) %>% sub(".*[1-9]_|.*spike-\\d{2,3}_", "", .) %>% sub(".*[[:alpha:]]_", "", .) %>% as.numeric(),
                 stringsAsFactors = FALSE)
#read in data
all_data <- readMSData(files = fp, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")

#split MS1 and MS2
data <- all_data[all_data@featureData@data$msLevel == 1]
data_ms2 <- all_data

#parallelization - get errors in BiocParallel otherwise
register(SerialParam())

#3. Pick peaks ------------------------------------------------------
cwp<-CentWaveParam()
cwp@ppm<-1.2
cwp@peakwidth<-c(10,60)
cwp@snthresh<-20
cwp@noise <- 5000
cwp@prefilter <- c(3, 1000)

data.pks <- findChromPeaks(data, param = cwp)
save(data.pks, file = "data.pks.RData")

#4: Group peaks to create "features"---------------------------------
#parameters
pdp <- PeakDensityParam(sampleGroups = data$group, binSize = 0.005, bw = 6, minFraction = 2/3) 
data.pks.grp <- groupChromPeaks(data.pks, param = pdp)
save(data.pks.grp, file = "data.pks.grp.RData")

#5: Fill missing peaks-----------------------------------
fpp <- FillChromPeaksParam()
data.pks.grp.fill <- fillChromPeaks(data.pks.grp)

save(data.pks.grp.fill, file = "data.pks.grp.fill.RData")
res <- data.pks.grp.fill

#6: Save diffreport of xdata ---------------------
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
an <- groupFWHM(an, perfwhm = 0.6)

##Annotate isotope peaks
#Mzabs = the allowed m/z error
an <- findIsotopes(an, mzabs=0.01)

##Peak grouping after correlation information into pseudospectrum groups 
#cor_eic_th = correlation threshold for EIC correlation
an <- groupCorr(an, cor_eic_th=0.75)

##Find adducts
an <- findAdducts(an, polarity="positive")

save(an, file = "an.RData")

save.image("smallmem2-RData.RData")




