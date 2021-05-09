setwd("/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_mannan/orbitrap/data/a-mannan")
load("./analysis/RData/RData_20210224.RData")




#1: Install packages --------------------------------------------------------
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
library(cowplot)
library(ggsci)


#2. Import and inspect MS data --------------------------------------------------------
#set working directory
setwd("/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_mannan/orbitrap/data/a-mannan")

#get file paths to mzML files
#keep only ones run only in positive mode
raw_files_path <- dir(path = "./mzML-files", 
                      all.files = FALSE, 
                      full.names = TRUE)
raw_files_path <- raw_files_path[grep("20201222", 
                                      raw_files_path)]
#also get negative and positive controls
#get neg acetone precipitation and procainamide labelling
#also get two neg solvent blanks, one run with method for controls 
#and one run with method for samples
control_files_path <- dir(path = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_mannan/orbitrap/data/controls/mzML-files/20201223", 
                          all.files = FALSE, 
                          full.names = TRUE)
control_files_path <- control_files_path[grep("neg|P1000|_07|_13",
                                              control_files_path)]

fp <- c(raw_files_path, 
        control_files_path)

rm(raw_files_path, 
   control_files_path)

#create phenodata data.frame
#each sample must have a unique name!
pd <- data.frame(name = basename(fp) %>%
                     sub("MS31_20201222_", "", .) %>% 
                     sub("-procA-50ng", "", .) %>% 
                     sub("-rep", "", .) %>% 
                     sub("pos_", "", .) %>% 
                     sub("SolventBlank", "solvent blank",.) %>% 
                     sub("P1000.*", "standard mix", .) %>% 
                     sub("neg_", "blank ", .) %>% 
                     sub("-", " ", .) %>% 
                     sub("_07", " 1", .) %>% 
                     sub("_13", " 2", .) %>% 
                     sub("_\\d\\d.mzML|.mzML", "", .) %>% 
                     gsub("_", "+", .),
                 sample_type = basename(fp) %>% 
                     sub("_\\d\\d.mzML", "", .) %>% 
                     sub("MS31_20201222_", "", .) %>% 
                     sub("pos_.*", "spiked pos", .) %>% 
                     sub("amannan_gh92_spike.*", "sample", .) %>% 
                     sub("amannan_spike.*|^gh92_spike.*", "spiked neg", .) %>% 
                     sub("P1000.*", "pos", .) %>% 
                     sub("neg_.*|SolventBlank", "neg", .),
                 replicate = basename(fp) %>% 
                     sub("_\\d\\d.mzML", "", .) %>% 
                     sub("MS31_20201222_", "", .) %>% 
                     sub(".*-rep", "", .) %>% 
                     sub("-procA-50ng", "", .) %>% 
                     sub(".*\\D+.*", "NA", .),
                 method = basename(fp) %>% 
                     sub(".mzML", "", .) %>% 
                     sub("MS31_20201222_.*_", "", .) %>% 
                     as.numeric() %>% 
                     sub("[8-9]|10|11|12|13", "pos mode top 5 MSMS", .) %>% 
                     sub("^[3-7]", "pos mode inclusion list MSMS", .),
                 stringsAsFactors = FALSE)

#read in data
data <- readMSData(files = fp, 
                   pdata = new("NAnnotatedDataFrame", 
                               pd), 
                   mode = "onDisk")

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

#4: Peak picking (CentWave) ---------------------------
cwp<-CentWaveParam()
cwp@ppm<-1.6
cwp@peakwidth<-c(10,100)
cwp@snthresh<-5

data_peaks<-findChromPeaks(data, 
                           param=cwp)

#5: Group peaks to create "features"---------
#parameters
pdp <- PeakDensityParam(sampleGroups = data$sample_type,
                        binSize = 0.005,
                        bw = 6) 

data_peaks_grouped <- groupChromPeaks(data_peaks, param = pdp)

#6: Fill in missing peaks----------
fpp <- FillChromPeaksParam()
data_peaks_grouped_filled <- fillChromPeaks(data_peaks_grouped)

#7: Save diffreport of xdata -----
#get ms2 data info
j <- grep("neg_", 
          fileNames(data), 
          invert = TRUE) #file indices
MS2.file.paths <- fileNames(data)[j] # file paths of files
MS2.file.names <- gsub(".*/", "", MS2.file.paths) # names of files

#split data objects
data_ms2 <- filterFile(data_peaks_grouped_filled, 
                       MS2.file.names)
data_ms1 <- data_peaks_grouped_filled[
    data_peaks_grouped_filled@featureData@data$msLevel == 1
    ]

xset <- as(data, "xcmsSet")
sampnames(xset) <- pData(data)$name
sampclass(xset) <- pData(data)$sample_type

#8. Isotope picking and filtering ----
##create xsannotate object
#extracts the peaktable from a provided xcmsSet,
#which is used for all further analysis

an <- xsAnnotate(xset)

##Group peaks of a xsAnnotate object according to their retention time 
#Perfwhm = parameter defines the window width, which is used for matching
anF <- groupFWHM(an, 
                 perfwhm = 0.6)

##Annotate isotope peaks
#Mzabs = the allowed m/z error
anI <- findIsotopes(anF, 
                    mzabs=0.01)

##Peak grouping after correlation information into pseudospectrum groups 
#cor_eic_th = correlation threshold for EIC correlation
anIC <- groupCorr(anI, 
                  cor_eic_th=0.75)

##Find adducts
anFA <- findAdducts(anIC, 
                    polarity="positive")

write.csv(getPeaklist(anFA), 
          file="./analysis/processing_tables/peaklist_xsannotate.csv")
pl <-getPeaklist(anFA)

pl_red <- getReducedPeaklist(anFA)

#filter by retention time 
#remove everything before 5 min and everything after 35 min
pl_rt <- pl %>% 
    filter(between(rt,
                   300,
                   2100))

#Filter for isotopes
pl_rt_iso <-pl_rt[pl_rt$isotopes!="",]

#Filter by blank exclusion
pl_rt_iso_be <-pl_rt_iso[pl_rt_iso$neg==0,]
pl_rt_be <-pl_rt[pl_rt$neg==0,]

#make rownames from rt and mz of features
rownames(pl_rt_iso_be)<-paste(round(pl_rt_iso_be$rt,1),
                              round(pl_rt_iso_be$mz,3),
                              sep="_")

rownames(pl_rt_be)<-paste(round(pl_rt_be$rt,1),
                          round(pl_rt_be$mz,3),
                          sep="_")

#change NA to 0
pl_rt_iso_be[is.na(pl_rt_iso_be)] <- 0
pl_rt_be[is.na(pl_rt_be)] <- 0

#subset to see peaks in samples
pl_rt_iso_be_samp <- pl_rt_iso_be[pl_rt_iso_be$sample >= 1,]
pl_rt_be_samp <- pl_rt_be[pl_rt_be$sample >= 1,]

#filter for isotopes or adducts
pl_rt_be_samp_isoadd <- pl_rt_be_samp[pl_rt_be_samp$isotopes!=""| 
                                          pl_rt_be_samp$adduct!="",]

pl_rt_be_isoadd <- pl_rt_be[pl_rt_be$isotopes!=""| 
                                pl_rt_be$adduct!="",]


#filter for rt of interest in samples
pl_rt_be_samp_isoadd_550to600 <- pl_rt_be_samp_isoadd %>% 
    filter(between(rt, 550, 600))

pl_rt_be_isoadd_470to520 <- pl_rt_be_isoadd %>% 
    filter(between(rt, 470, 520))

write.csv(pl_rt_be_samp_isoadd_550to600, 
          file="./analysis/processing_tables/peaklist_samples_rt550to600.csv")

write.csv(pl_rt_be_isoadd_470to520, 
          file="./analysis/processing_tables/peaklist_rt470to520.csv")

#PLOTTING FLR EIC MANNAN SAMPLES ONLY ---------

#1:get features----
pl_rt_be_isoadd_man <- pl_rt_be_isoadd[pl_rt_be_isoadd$sample >= 1,]
mannan_peaks <- pl_rt_be_isoadd_man %>% 
    filter(solvent.blank.1 == 0 & 
               solvent.blank.2 == 0 & 
               blank.acetone.precipitation == 0 & 
               blank.procainamide.reaction == 0 &
               amannan.gh92.spike2 > 1000 &
               amannan.gh92.spike3 > 1000 )
mannan_peaks <- mannan_peaks %>% 
    filter(amannan.spike3 > 1000 & 
               gh92.spike > 10000 & 
               spike > 1000 |
               amannan.spike3 == 0 & 
               gh92.spike == 0  & 
               spike == 0 )
mannan_peaks <- cbind(rt_round = round_any(mannan_peaks$rt, 
                                           5),
                      mannan_peaks)

##collapse features with multiple isotopes
setDT(mannan_peaks)
#split out features without an isotope detected
mannan_peaks_noiso <- mannan_peaks[mannan_peaks$isotopes=="",]
mannan_peaks_iso <- mannan_peaks[!mannan_peaks$isotopes=="",]
#make column for the isotope group
mannan_peaks_iso$isotope_group <- mannan_peaks_iso$isotopes %>% 
    sub("\\[M.*", "", .)
#order isotopes within each group correctly
mannan_peaks_iso$isotope_number <- mannan_peaks_iso$isotopes %>% 
    sub(".*\\[M\\].*", "0", .) %>% 
    sub(".*\\[M\\+", "", .) %>% 
    sub("\\].*", "", .) %>% 
    as.numeric()
mannan_peaks_iso <- mannan_peaks_iso[order(isotope_group, 
                                           isotope_number),]
#get concatenated list of isotopes per group
iso_concat <- mannan_peaks_iso[, 
                               list(isotopes = paste(isotopes, 
                                                     collapse = ', ')), 
                               by = isotope_group]
#remove duplicates within each isotope group (will keep [M] isotope)
#because of ordering
mannan_peaks_iso <- unique(mannan_peaks_iso, 
                           by = "isotope_group")
#merge to get concatenated isotope lists
mannan_peaks_iso <- merge(mannan_peaks_iso,
                          iso_concat,
                          by = "isotope_group")
#clean up df
mannan_peaks_iso <- mannan_peaks_iso %>% 
    select(-c("isotope_group",
              "isotope_number",
              "isotopes.x"))
names(mannan_peaks_iso)[names(mannan_peaks_iso) == 'isotopes.y'] <- 'isotopes'

#replace features that don't contain [M] isotope with [M] isotope
temp <- mannan_peaks_iso %>% 
    filter(!grepl("\\[M\\]", isotopes))

mannan_peaks_iso <- mannan_peaks_iso %>% 
    filter(grepl("\\[M\\]", isotopes))

temp.vec <- temp$isotopes %>% 
    sub("\\[M.*", "", .)

pl_rt_be_isoadd_man$isotope_group <- pl_rt_be_isoadd_man$isotopes %>% 
    sub("\\[M.*", "", .)
pl_rt_be_isoadd_man$isotope_number <- pl_rt_be_isoadd_man$isotopes %>% 
    sub(".*\\[M\\].*", "0", .) %>% 
    sub(".*\\[M\\+", "", .) %>% 
    sub("\\].*", "", .) %>% 
    as.numeric()

temp <- pl_rt_be_isoadd_man[pl_rt_be_isoadd_man$isotope_group %in% temp.vec,] %>% 
    filter(isotope_number == 0)

temp$isotope_group <- NULL
temp$isotope_number <- NULL

temp <- cbind(rt_round = round_any(temp$rt, 
                                   5),
              temp)

mannan_peaks_iso <- rbind(mannan_peaks_iso, 
                          temp)

rm(temp,
   temp.vec)

#merge features with and without isotopes
mannan_peaks <- rbind.fill(mannan_peaks_noiso,
                           mannan_peaks_iso)
rm(mannan_peaks_noiso,
   mannan_peaks_iso,
   iso_concat)

setDT(mannan_peaks)
mannan_peaks <- mannan_peaks[order(rt_round,
                                   mz),]

write.csv(mannan_peaks,
          file = "./analysis/analysis_tables/mannan_peaklist.csv") 

##annotate features with predicted ions
#import table
man_mz_predicted <- fread("mannan_mass-list.txt",
                          blank.lines.skip = TRUE)

#remove "extra" columns
extraCol <- c('formula',
              'mass',
              'charge')

man_predicted <- man_mz_predicted %>% 
    select(-all_of(extraCol))
#make data.table
setDT(man_predicted)
setDT(mannan_peaks)
#create interval to overlap with (same width as for peak grouping)
man_predicted$mz <- as.numeric(man_predicted$mz)
man_predicted$mzmin <- man_predicted$mz-0.005
man_predicted$mzmax <- man_predicted$mz+0.005
#match using foverlaps from data.table (very fast)
setkey(man_predicted, mzmin, mzmax)
mannan_peaks <- foverlaps(mannan_peaks,
                          man_predicted)

#change NA values created during matching (features with no match) to be blank
mannan_peaks <- mannan_peaks %>% 
    replace_na(list("sugar"="unknown",
                    "dp"="unknown",
                    "id" = "", 
                    "ion"= "unknown", 
                    "mz" = "", 
                    "mzmin" = "",
                    "mzmax"= ""))

#format ion names for plot
mannan_peaks$ion <- mannan_peaks$ion %>% 
    sub("\\+1", "+", .)

mannan_peaks$id_ion <- paste0(mannan_peaks$id,
                              ":",
                              mannan_peaks$ion,
                              " mz=",
                              round(mannan_peaks$i.mz,
                                    3))

mannan_peaks$id_ion <- mannan_peaks$id_ion %>% 
    sub("^:", "", .)

mannan_peaks <- mannan_peaks[order(rt_round),] 

fwrite(mannan_peaks, 
       file = "./analysis/analysis_tables/mannan_peaklist_matched.txt",
       sep = "\t")


#2:extract and format eic -----
##extract eic 
man_mz.found.vector <- mannan_peaks$i.mz %>% 
    round(., 3) %>% 
    unique()

man_ions.found.vector  <- mannan_peaks$id_ion %>% 
    unique()

man_ions.found.vector <- man_ions.found.vector %>% 
    sub("unknown", "unknown:", .) %>% 
    sub("K2", "k-carrageenan DP2", .)

data_man <- filterFile(data,
                       file = which(grepl("spike", data$name)))

man.names <- grep("spike", data$name, value = TRUE)
man.groups <- pd$sample_type[which(grepl("spike", data$name))]

#get chromatograms
man_chr_list <- list()
error = 0.001

for (i in 1:length(man_mz.found.vector)){
    mzr = c(man_mz.found.vector[i] - error,
            man_mz.found.vector[i] + error)
    man_chr_list[[i]] <- chromatogram(data_man, 
                                      mz = mzr)
}

#extract intensity and rt values
man_chr_int_list <- list()
for (i in 1:length(man.names)){
    man_chr_int_list[[i]] <- lapply(man_chr_list, function(x) {
        x[[i]]@intensity
    }) 
}

man_chr_rt_list <- list()
for (i in 1:length(man.names)){
    man_chr_rt_list[[i]] <- lapply(man_chr_list, function(x) {
        x[[i]]@rtime
    }) 
}

#build data frame (long format)
man.df <- data.frame(ion = as.character(),
                     sample = as.character(),
                     group = as.character(),
                     rt = as.numeric(),
                     intensity = as.numeric())

for (i in 1:length(man.names)){
    for (j in 1:length(man_mz.found.vector)){
        rt = man_chr_rt_list[[i]][[j]]
        intensity = man_chr_int_list[[i]][[j]]
        sample = rep(man.names[i], length(rt))
        group = rep(man.groups[i], length(rt))
        ion = rep(man_ions.found.vector[j], length(rt))
        temp <- data.frame(ion = ion,
                           sample = sample,
                           group = group,
                           rt = rt,
                           intensity = intensity)
        man.df <- rbind(man.df,
                        temp)
    }
}

man.df[is.na(man.df)] <- 0

#set variables as factors
man.df$ion <- factor(man.df$ion,
                     levels = man_ions.found.vector)
man.df$group_fmt <- man.df$group %>% 
    sub("spiked neg", 
        "k-carrageenan spiked negative controls", 
        .) %>% 
    sub("spiked pos", 
        "k-carrageenan spike only control", 
        .) %>% 
    sub("sample",
        "k-carrageenan spiked a-mannan GH92 digests", 
        .)
man.df$group_fmt <- factor(man.df$group_fmt,
                           levels = unique(man.df$group_fmt))
man.df$rt_min <- man.df$rt / 60

man.df$rep <- man.df$sample %>% 
    sub("amannan\\+gh92\\+spike", "", .) %>% 
    sub("^2", "1", .) %>% 
    sub("^3", "2", .) %>% 
    sub("^\\D.*", "1", .)
man.df$rep <- factor(man.df$rep,
                     levels = c(1,2))


#identify and remove ions with zero intensity in samples
mannan_zeroIntensityIons <- man.df %>% 
    group_by(sample, ion) %>% 
    summarise(sum = sum(intensity), .groups="keep") %>% 
    filter(sum == 0)

mannan_zeroIntensityIons$sample_ion <- paste0(mannan_zeroIntensityIons$sample,
                                              "_",
                                              mannan_zeroIntensityIons$ion)

mannan_zeroIntensityIons <- mannan_zeroIntensityIons$sample_ion


man.df$sample_ion <- paste0(man.df$sample,
                            "_",
                            man.df$ion)

man.df <- subset(man.df, 
                 !(sample_ion  %in% mannan_zeroIntensityIons))



#3:import FLD and format-----------
#read in text files
fld_fp <- dir(path = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_mannan/triple-quad/Data_final/FLD-files", 
              all.files = FALSE, 
              full.names = TRUE)
fld_fp <- fld_fp[grep("spike" ,
                      fld_fp)]
fld_man_list <- lapply(fld_fp, 
                       fread)

#make data frame
fld_man.groups <- fld_fp %>% 
    sub(".*mannan_gh82_spiked.*", "sample", .) %>% 
    sub(".*mannan_spiked.*", "spike neg", .) %>% 
    sub(".*spike_pos.*", "spiked pos", .)


fld_man.rep <- fld_fp%>% 
    sub(".*mannan_gh82_spiked_proca_rep1.*", "1", .) %>% 
    sub(".*mannan_gh82_spiked_proca_rep2.*", "2", .) %>% 
    sub(".*.txt", "1", .)

fld_man.sample <- basename(fld_fp)%>% 
    sub("2020.*HILIC_", "", .) %>% 
    sub("-FLD.txt", "", .)

fld_man.df <- data.frame(sample = as.character(),
                         group = as.character(),
                         rep = as.numeric(),
                         rt = as.numeric(),
                         intensity = as.numeric())

for (i in 1:length(fld_man.sample)){
    rt = fld_man_list[[i]]$V1
    intensity = fld_man_list[[i]]$V2
    sample = rep(fld_man.sample[i], length(rt))
    group = rep(fld_man.groups[i], length(rt))
    rep = rep(fld_man.rep[i], length(rt))
    temp <- data.frame(sample = sample,
                       group = group,
                       rep = rep,
                       rt = rt,
                       intensity = intensity)
    fld_man.df <- rbind(fld_man.df,
                        temp)
}



#transform FLD on x axis - done manually to fit peaks
fld_man.df$rt_trans <- fld_man.df$rt - 0.15

#get FLD minimum and offset each sample so that minimum is zero
#rt range at which minimum is found based on initial plots
fld_man.df.notzeroed <- fld_man.df #keep to be safe

for (i in 1:length(fld_man.sample)){
    x <- fld_man.df.notzeroed[fld_man.df.notzeroed$sample==fld_man.sample[i],]
    int.min = x %>% 
        filter(between(rt, 5, 25)) %>% 
        select(intensity) %>% 
        min()
    x$intensity <- x$intensity + abs(int.min)
    fld_man.df$intensity[fld_man.df$sample==fld_man.sample[i]] <- x$intensity
}

#set variables as factors
fld_man.df$rep <- factor(fld_man.df$rep,
                         levels = c(1,2))
fld_man.df$group_fmt <- fld_man.df$group %>% 
    sub("spike neg", 
        "k-carrageenan spiked negative controls", 
        .) %>% 
    sub("spiked pos", 
        "k-carrageenan spike only control", 
        .) %>% 
    sub("sample",
        "k-carrageenan spiked a-mannan GH92 digests", 
        .)
fld_man.df$group_fmt <- factor(fld_man.df$group_fmt,
                               levels = unique(fld_man.df$group_fmt))

#4: plot FLR ----

fld_man.df$var <- "FLR"

fld_man.df$group_fmt <- factor(fld_man.df$group_fmt,
                               levels = c("k-carrageenan spiked a-mannan GH92 digests",
                                          "k-carrageenan spiked negative controls",
                                          "k-carrageenan spike only control"))

fld_man.df$group_fmt <- factor(fld_man.df$group_fmt,
                               levels = rev(levels(fld_man.df$group_fmt)))


#make breaks
x <- seq(5, 15, 0.5)
major_breaks_zoom <- vector(length = length(x),
                            mode = "character")

for (i in 1:length(major_breaks_zoom)){
    if (x[i] %% 1 == 0){
        major_breaks_zoom[i] <- as.character(x[i])
    } else if (x[i] %% 1 != 0) {
        major_breaks_zoom[i] <- ""
    }
}


#filter by retention time

fld_man.df.zoom <- fld_man.df %>% 
    filter(between(rt_trans, 5, 25))

fld_man.df.zoom$group2 <- fld_man.df.zoom$group %>% 
    sub("spiked.*", "neg", .)
fld_man.df.zoom$group2 <- factor(fld_man.df.zoom$group2,
                                 levels = c("sample",
                                            "neg"))


#make palette
pal <- c("#182031",
         "#21a4b0ff")

tiff("./analysis/test1.tiff", 
     res = 300, 
     height = 3, 
     width = 12, 
     units = "in")

#plot
man_p1 <- ggplot() +
    geom_line(mapping = aes(rt_trans,
                            intensity,
                            colour = group2,
                            group = sample),
              data = fld_man.df.zoom,
              lwd = 1.2) +
    scale_colour_manual(values = pal) +
    facet_grid(rows = vars(var)) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 12),
          legend.position = "none",
          strip.text = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    scale_x_continuous(breaks = seq(5, 15, 0.5),
                       labels = major_breaks_zoom,
                       limits = c(5,15),
                       expand = c(0,0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                       labels = scales::scientific)
dev.off()


#5: plot  eic ----
man.df$var <- man.df$group %>% 
    sub("spiked.*", "neg", .)
man.df$var <- factor(man.df$var,
                     levels = c("sample",
                                "neg"))

pal <- c("#182031",
         "#21a4b0ff")

tiff("./analysis/test2.tiff", 
     res = 300, 
     height = 3, 
     width = 12, 
     units = "in")
man_p2 <- ggplot() +
    geom_line(mapping = aes(rt_min,
                            intensity,
                            group = sample_ion,
                            colour = var),
              data = man.df %>% 
                  filter(ion == "k-carrageenan DP2:[M+H]+ m/z=624.242" |
                             ion == "sulphated monosaccharide:[M+H]+ m/z=480.201"),
              lwd = 1.5) +
    facet_grid(rows = vars(var)) +
    scale_colour_manual(values = pal) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 12),
          legend.position = "none",
          strip.text.x = element_blank(),
          strip.text.y = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    scale_x_continuous(breaks = seq(5, 15, 0.5),
                       labels = major_breaks_zoom,
                       limits = c(5,15),
                       expand = c(0,0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                       labels = scales::scientific,
                       n.breaks = 4)
dev.off()


#7:plot togther -----
g1 <- ggplotGrob(man_p1)
g2 <- ggplotGrob(man_p2)

g <- rbind(g1, 
           g2, 
           size = "first")

g$widths <- unit.pmax(g1$widths, 
                      g2$widths)

tiff("./analysis/summary_plot_v3.tiff", 
     res = 300, 
     height = 6, 
     width = 8, 
     units = "in")
grid.newpage()
grid.draw(g)
dev.off()

svg("./analysis/summary_plot_v3.svg", 
    height = 6, 
    width = 6)
grid.newpage()
grid.draw(g)
dev.off()


#PLOTTING FLR EIC CONTROLS ONLY ---------
#1: get features----
pl_rt_be_isoadd_pos <- pl_rt_be_isoadd[pl_rt_be_isoadd$standard.mix >= 1,]
pos_control_peaks <- pl_rt_be_isoadd_pos %>% 
    filter(solvent.blank.1 == 0 & 
               solvent.blank.2 == 0 & 
               blank.acetone.precipitation == 0 & 
               blank.procainamide.reaction == 0)

pos_control_peaks <- cbind(rt_round = round_any(pos_control_peaks$rt, 
                                                5),
                           pos_control_peaks)

pl_rt_be <- cbind(rt_round = round_any(pl_rt_be$rt, 
                                       5),
                  pl_rt_be)

##collapse features with multiple isotopes
setDT(pos_control_peaks)
#split out features without an isotope detected
pos_control_peaks_noiso <- pos_control_peaks[pos_control_peaks$isotopes=="",]
pos_control_peaks_iso <- pos_control_peaks[!pos_control_peaks$isotopes=="",]
#make column for the isotope group
pos_control_peaks_iso$isotope_group <- pos_control_peaks_iso$isotopes %>% 
    sub("\\[M.*", "", .)
#order isotopes within each group correctly
pos_control_peaks_iso$isotope_number <- pos_control_peaks_iso$isotopes %>% 
    sub(".*\\[M\\].*", "0", .) %>% 
    sub(".*\\[M\\+", "", .) %>% 
    sub("\\].*", "", .) %>% 
    as.numeric()
pos_control_peaks_iso <- pos_control_peaks_iso[order(isotope_group, 
                                                     isotope_number),]
#get concatenated list of isotopes per group
iso_concat <- pos_control_peaks_iso[, 
                                    list(isotopes = paste(isotopes, 
                                                          collapse = ', ')), 
                                    by = isotope_group]
#remove duplicates within each isotope group (will keep [M] isotope)
#because of ordering
pos_control_peaks_iso <- unique(pos_control_peaks_iso, 
                                by = "isotope_group")
#merge to get concatenated isotope lists
pos_control_peaks_iso <- merge(pos_control_peaks_iso,
                               iso_concat,
                               by = "isotope_group")
#clean up df
pos_control_peaks_iso <- pos_control_peaks_iso %>% 
    select(-c("isotope_group",
              "isotope_number",
              "isotopes.x"))
names(pos_control_peaks_iso)[names(pos_control_peaks_iso) == 'isotopes.y'] <- 'isotopes'

#merge features with and without isotopes
pos_control_peaks <- rbind.fill(pos_control_peaks_noiso,
                                pos_control_peaks_iso)
rm(pos_control_peaks_noiso,
   pos_control_peaks_iso,
   iso_concat)

write.csv(pos_control_peaks,
          file = "./analysis/analysis_tables/pos_control_peaklist.csv") 

##identify features based on predictions
#import table
pos_mz_predicted <- fread("poscontrol_mass-list.txt",
                          blank.lines.skip = TRUE)

#remove "extra" columns
extraCol <- c('formula',
              'mass',
              'charge')

predicted <- pos_mz_predicted %>% 
    select(-all_of(extraCol))
#make data.table
setDT(predicted)
setDT(pos_control_peaks)
#create interval to overlap with (same width as for peak grouping)
predicted$mz <- as.numeric(predicted$mz)
predicted$mzmin <- predicted$mz-0.005
predicted$mzmax <- predicted$mz+0.005
#match using foverlaps from data.table (very fast)
setkey(predicted, mzmin, mzmax)
pos_control_peaks_nopred <- pos_control_peaks
pos_control_peaks <- foverlaps(pos_control_peaks,
                               predicted)

#change NA values created during matching (features with no match) to be blank
pos_control_peaks <- pos_control_peaks %>% 
    replace_na(list("sugar"="",
                    "dp"="",
                    "id" = "", 
                    "ion"= "", 
                    "mz" = "", 
                    "mzmin" = "",
                    "mzmax"= ""))
#only keep matched features
pos_control_peaks_matched <- pos_control_peaks[!pos_control_peaks$sugar=="",]

#format ion names for plot
pos_control_peaks_matched$ion <- pos_control_peaks_matched$ion %>% 
    sub("\\+ProcA", "", .) %>% 
    sub("\\+1", "+", .)

pos_control_peaks_matched$id_ion <- paste0(pos_control_peaks_matched$id,
                                           ":",
                                           pos_control_peaks_matched$ion,
                                           " mz=",
                                           round(pos_control_peaks_matched$i.mz,
                                                 3))

pos_control_peaks_matched <- pos_control_peaks_matched[order(rt_round),] 

fwrite(pos_control_peaks_matched, 
       file = "./analysis/analysis_tables/pos_control_peaklist_matched.txt",
       sep = "\t")

#2:extract and format eic -----
##extract eic 
mz.found.vector <- pos_control_peaks_matched$i.mz %>% 
    round(., 3) %>% 
    unique()

ions.found.vector  <- pos_control_peaks_matched$id_ion %>% 
    unique()

ions.found.vector <- ions.found.vector %>% 
    sub("K2", "k-carrageenan DP2", .) %>% 
    sub("G1", "glucose", .) %>%
    sub("L2", "laminaribiose", .) %>%
    sub("L4", "laminaritetraose", .) %>%
    sub("K4", "k-carrageenan DP4", .) %>%
    sub("BM3", "b-mannotriose", .)
    

data_controls <- filterFile(data,
                            file = which(grepl("standard|acetone|procainamide",
                                               data$name)))


control.names <- data_controls$name
control.groups <- data_controls$sample_type
control.groups <- control.groups %>% 
    sub("neg","negative controls", .) %>% 
    sub("pos","oligosaccharide standards mix",.)

#get chromatograms
control_chr_list <- list()
error = 0.001

for (i in 1:length(mz.found.vector)){
    mzr = c(mz.found.vector[i] - error,
            mz.found.vector[i] + error)
    control_chr_list[[i]] <- chromatogram(data_controls, 
                                          mz = mzr)
}

#extract intensity and rt values
control_chr_int_list <- list()
for (i in 1:length(control.names)){
    control_chr_int_list[[i]] <- lapply(control_chr_list, function(x) {
        x[[i]]@intensity
    }) 
}

control_chr_rt_list <- list()
for (i in 1:length(control.names)){
    control_chr_rt_list[[i]] <- lapply(control_chr_list, function(x) {
        x[[i]]@rtime
    }) 
}

#build data frame (long format)
control.df <- data.frame(ion = as.character(),
                         sample = as.character(),
                         group = as.character(),
                         rt = as.numeric(),
                         intensity = as.numeric())

for (i in 1:length(control.names)){
    for (j in 1:length(mz.found.vector)){
        rt = control_chr_rt_list[[i]][[j]]
        intensity = control_chr_int_list[[i]][[j]]
        sample = rep(control.names[i], length(rt))
        group = rep(control.groups[i], length(rt))
        ion = rep(ions.found.vector[j], length(rt))
        temp <- data.frame(ion = ion,
                           sample = sample,
                           group = group,
                           rt = rt,
                           intensity = intensity)
        control.df <- rbind(control.df,
                            temp)
    }
}

control.df[is.na(control.df)] <- 0

#set variables as factors
control.df$ion <- factor(control.df$ion,
                         levels = ions.found.vector)

control.df$ion <- relevel(control.df$ion, 
                          "k-carrageenan DP2:[M-SO3+H]+ mz=544.286")
control.df$ion <- relevel(control.df$ion, 
                          "k-carrageenan DP2:[M+H]+ mz=624.242")
control.df$ion <- relevel(control.df$ion, 
                          "k-carrageenan DP2:[M+Na]+ mz=646.225")
control.df$ion <- relevel(control.df$ion, 
                          "k-carrageenan DP2:[M+K]+ mz=662.199") 
#come before G1 in chromatography


control.df$group <- factor(control.df$group,
                           levels = c("oligosaccharide standards mix",
                                      "negative controls"))

control.df$rt_min <- control.df$rt / 60

#3: import FLD and format-----------
#read in text files
fld_fp <- dir(path = "/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_mannan/triple-quad/Data_final/FLD-files", 
              all.files = FALSE, 
              full.names = TRUE)
fld_fp <- fld_fp[grep("pos1000" ,
                      fld_fp)]
fld_pos <- fread(fld_fp)
fld_pos$V3 <- NULL
names(fld_pos) <- c("rt",
                    "intensity")

fld_pos$group <- "oligosaccharide standards mix"
fld_pos$group <- factor(fld_pos$group,
                        levels = "oligosaccharide standards mix")

#get FLD minimum and offset so that minimum is zero
#rt range at which minimum is found based on initial plots
fld_pos.notzeroed <- fld_pos #keep to be safe
int.min = fld_pos.notzeroed %>% 
    filter(between(rt, 5, 20)) %>% 
    select(intensity) %>% 
    min()
fld_pos$intensity <- fld_pos$intensity + abs(int.min)

#transform FLD on x axis - done manually to fit peaks
fld_pos$rt_trans <- fld_pos$rt - 1.15

#4:plot flr only: pos control ----

fld_pos$var <- "FLR"

#ZOOMED
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
fld_pos.zoom <- fld_pos %>% 
    filter(between(rt_trans, 5, 25))

control_flr <- ggplot() +
    geom_line(mapping = aes(rt_trans,
                            intensity),
              data = fld_pos.zoom,
              colour = "black",
              lwd = 1.2) +
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
    theme(strip.background = element_blank(),
          strip.text.y = element_text(angle = 360, hjust = 0),
          strip.text.x = element_blank(),
          axis.text = element_text(size = 12)) + 
    scale_x_continuous(breaks = seq(5, 25, 1),
                       labels = major_breaks_zoom,
                       limits = c(5, 25),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)))

#5:plot all eic ----
#get the most abundant ion for each sugar
pos_control_peaks_matched <- pos_control_peaks_matched[order(id, 
                                                             standard.mix,
                                                             decreasing = TRUE),]

pos_control_peaks_matched$id_ion_fmt <- pos_control_peaks_matched$id_ion  %>% 
    sub("K2", "k-carrageenan DP2", .) %>% 
    sub("G1", "glucose", .) %>%
    sub("L2", "laminaribiose", .) %>%
    sub("L4", "laminaritetraose", .) %>%
    sub("K4", "k-carrageenan DP4", .) %>%
    sub("BM3", "b-mannotriose", .)

pos_control_abundant_ions <- unique(pos_control_peaks_matched, 
                                    by = "id")
    
pos_control_abundant_ions <- pos_control_abundant_ions$id_ion_fmt

control.df_abundant_ions <- control.df[control.df$ion %in% 
                                           pos_control_abundant_ions,]

#change mz to m/z
control.df_abundant_ions$ion <- control.df_abundant_ions$ion %>% 
    sub("mz", "m/z", .)

#set factor levels
control.df_abundant_ions$ion <- factor(control.df_abundant_ions$ion,
                                       levels = unique(control.df_abundant_ions$ion))
control.df_abundant_ions$ion <- relevel(control.df_abundant_ions$ion, 
                                        "k-carrageenan DP2:[M-SO3+H]+ m/z=544.286")

control.df_abundant_ions$var <- "EIC"

#set palette
pal_cold<- hcl.colors(n = 5, "Cold")
pal_control <- c("#552046",
                 pal_cold[1:3],
                 "#339D84",
                 "#174844")

names(pal_control) <- levels(control.df_abundant_ions$ion)

#plot
pcontrol_eic <- ggplot() +
    geom_line(mapping = aes(rt_min,
                            intensity,
                            colour = ion),
              data = control.df_abundant_ions[control.df_abundant_ions$group== 
                                                  "oligosaccharide standards mix",],
              lwd = 1.2) +
    scale_color_viridis(discrete = TRUE, option = "D")+
    theme(strip.background = element_blank(),
          strip.text.y = element_text(angle = 360, hjust = 0),
          strip.text.x = element_blank()) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 12),
          legend.position = "none",
          axis.text = element_text(size = 12)) +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    scale_x_continuous(breaks = seq(5, 25, 1),
                       labels = major_breaks_zoom,
                       limits = c(5,25),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                       labels = scales::scientific)
#NEGATIVE CONTROLS

#get y axis limits of positive control plot
pcontrol_ylim <- ggplot_build(control_p3)$layout$panel_params[[1]]$y.range

ncontrol_eic <- ggplot() +
    geom_line(mapping = aes(rt_min,
                            intensity,
                            colour = ion),
              data = control.df_abundant_ions[control.df_abundant_ions$group== 
                                                  "negative controls",],
              lwd = 1.2) +
    scale_color_viridis(discrete = TRUE, option = "D") +
    theme(strip.background = element_blank(),
          strip.text.y = element_text(angle = 360, hjust = 0),
          strip.text.x = element_blank()) +
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
    scale_x_continuous(breaks = seq(5, 25, 1),
                       labels = major_breaks_zoom,
                       limits = c(5,25),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = c(0,0),
                       labels = scales::scientific,
                       limits = pcontrol_ylim)
#7:plot all together -----
cg1 <- ggplotGrob(control_flr)
cg2 <- ggplotGrob(pcontrol_eic)
#cg3 <- ggplotGrob(ncontrol_eic)


cg <- rbind(cg1, 
            cg2, 
            size = "first")
cg$widths <- unit.pmax(cg1$widths, 
                       cg2$widths)

tiff("./analysis/analysis_plots/control_plot_v4.tiff", 
     res = 300, 
     height = 5, 
     width = 12, 
     units = "in")
grid.newpage()
grid.draw(cg)
dev.off()

svg("./analysis/analysis_plots/control_plot_v4.svg", 
     height = 6, 
     width = 12)
grid.newpage()
grid.draw(cg)
dev.off()


#EXTRACT AND PLOT MS2 -----
#get peaks to extract ms2 spectra for ----
##combine features for control and mannan
peaksFiltered <- rbind.fill(pos_control_peaks_matched,
                            mannan_peaks)
setDT(peaksFiltered)
peaksFiltered <- peaksFiltered[order(i.mz, rt),]
peaksFiltered <- unique(peaksFiltered,
                                by = c("i.mz", "rt"))

##get ms2 file names and indices
j <- which(fileNames(data) %in% fileNames(data_ms2)) #index of ms2 files
MS2.file.paths <- fileNames(data)[j] # file paths 
MS2.file.names <- basename(MS2.file.paths) # names 

##get feature definitions and values
fd <- featureDefinitions(data) # extract feature defs as a new object
fv <- featureValues(data) # extract feature values as a new object
fv.filtered <- fv[, 
                  colnames(fv) %in% MS2.file.names] #filter for MS2 samples

##get and filter chromPeak data
cp <- chromPeaks(data) #extract peaks
cp <- cbind(rowid=seq(nrow(cp)), 
            cp) # add a col for matching with 'peakidx' in feature definitions 
cp.filtered <- cp[which(cp[,which(colnames(cp) == "sample")] %in% j),] # filter cp to include only DDA samples









