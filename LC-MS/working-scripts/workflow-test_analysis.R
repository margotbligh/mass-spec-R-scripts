setwd("/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/1_standards/1E_full-workflow")
load("./analysis/RData/RData_20210114.RData")

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



library(faahKO)
library(pander)
library(RColorBrewer)
library(magrittr)
library(pheatmap)
library(msdata)
library(png)
library(tidyr)
library(detect)
library(devtools)
library(MetMSLine)
library(pcaMethods)
library(statTarget)
library(randomForest)
library(rlist)
library(purrr)
library(reshape2)
library(extrafont)
library(Rmisc)
library(edgeR)
library(limma)
library(mixOmics)
library(HTSFilter)
library(rstatix)
library(reshape2)
library(remotes)
library(ggridges)
library(gridExtra)
library(ggpubr)
library(lemon)
library(cowplot)
library(ggsci)


#2. Import and inspect MS data --------------------------------------------------------
#get file paths to mzML files
#keep only ones run only in positive mode
#also get two solvent blanks
fp <- dir(path = "mzML-files/20201223", 
          all.files = FALSE, 
          full.names = TRUE)
fp <- fp[grep("neg|P1000|_07|_13",
              fp)]


#create phenodata data.frame
#each sample must have a unique name!
pd <- data.frame(name = basename(fp) %>%
                     sub("MS31_20201222_", "", .) %>% 
                     sub("P1000.*", "standard mix", .) %>% 
                     sub("SolventBlank", "solvent blank",.) %>% 
                     sub("neg_", "blank ", .) %>% 
                     sub("-", " ", .) %>% 
                     sub("_07", " 1", .) %>% 
                     sub("_13", " 2", .) %>% 
                     sub("_\\d\\d.mzML|.mzML", "", .),
                 sample_type = basename(fp) %>% 
                     sub(".*P1000.*", "pos", .) %>% 
                     sub(".*neg_.*|.*SolventBlank.*", "neg", .),
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

#4: Peak picking (CentWave) ---------------------------
cwp<-CentWaveParam()
cwp@ppm<-1.6
cwp@peakwidth<-c(10,100)
cwp@snthresh<-5

data_peaks<-findChromPeaks(data, 
                           param=cwp)

data_ms2<-findChromPeaks(data_ms2, 
                         param=cwp)

save(data_peaks, 
     file = "./analysis/RData/data_peaks.RData")

#5: Group peaks to create "features"---------
#parameters
pdp <- PeakDensityParam(sampleGroups = data$sample_type,
                        binSize = 0.005,
                        bw = 6) 

data_peaks_grouped <- groupChromPeaks(data_peaks, param = pdp)

save(data_peaks_grouped, 
     file = "./analysis/RData/data_peaks_grouped.RData")

#6: Fill in missing peaks----------
fpp <- FillChromPeaksParam()
data_peaks_grouped_filled <- fillChromPeaks(data_peaks_grouped)

save(data_peaks_grouped_filled, 
     file = "./analysis/RData/data_peaks_grouped_filled.RData")

#7: Remove large data files from environment----------
rm(data, data_peaks, data_peaks_grouped)    

#8: Save diffreport of xdata -----
xset <- as(data_peaks_grouped_filled, "xcmsSet")
sampnames(xset) <- pData(data)$name
sampclass(xset) <- pData(data)$sample_type

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
                    polarity="positive")

#11. Peak list filtering----
#get peak list
pl <-getPeaklist(an)

#filter by retention time (5-35 min)
pl_rt <- pl %>% 
    filter(between(rt,
                   300,
                   2100))

#filter by blank exclusion (detected peaks)
pl_rt_be <-pl_rt[pl_rt$neg==0,]

#make rownames from rt and mz of features
rownames(pl_rt_be)<-paste(round(pl_rt_be$rt,1),
                          round(pl_rt_be$mz,3),
                          sep="_")
#change NA to 0
pl_rt_be[is.na(pl_rt_be)] <- 0

#filter by blank exclusion (intensity, accounts for filled peaks)
pos_control_peaks <- pl_rt_be %>% 
    filter(solvent.blank.1 <1e4 & 
               solvent.blank.2 < 1e4 & 
               blank.acetone.precipitation < 1e4 & 
               blank.procainamide.reaction < 1e4)

#add rounded retention time as first colum
pos_control_peaks <- cbind(rt_round = round_any(pos_control_peaks$rt, 
                                                5),
                           pos_control_peaks)

#12: Collapse features with multiple isotopes -----
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

#13: Annotate features based on predictions----
#import prediction table
#created with: sugarMassesPredict.py -dp 1 6 -p 0 -m sulphate carboxyl anhydrobridge -i pos -s 300 2000 -l procainamide
mz_predicted <- fread("predicted_sugars.txt")

#remove "extra" columns
extraCol <- c('dp',
              'mass',
              'formula')

mz_predicted <- mz_predicted %>% 
    select(-all_of(extraCol))

#make long format
predicted <- mz_predicted %>% 
    gather(key = "ion",
           value = "mz",
           -name)

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
#remove extra columns
pos_control_peaks$mzmin <- NULL
pos_control_peaks$mzmax <- NULL

pos_control_peaks <- pos_control_peaks %>% 
    replace_na(list("name"="",
                    "ion"= "", 
                    "mz" = ""))
#only keep matched features
pos_control_peaks_matched <- pos_control_peaks[!pos_control_peaks$name=="",]

#format ion names
pos_control_peaks_matched$id <- pos_control_peaks_matched$name %>% 
    sub("hex-1-procA", "glucose", .) %>% 
    sub("hex-2-procA", "laminaribiose", .) %>%
    sub("hex-3-procA", "b-mannotriose", .) %>%
    sub("hex-4-procA", "laminaritetraose", .) %>%
    sub("hex-2-sulphate-1-anhydrobridge-1-procA", "k-carrageenan DP2", .) %>%
    sub("hex-4-sulphate-2-anhydrobridge-2-procA", "k-carrageenan DP4", .)
    
#desulphated k-carrageenan DP2
pos_control_peaks_matched$ion[pos_control_peaks_matched$id 
                              == "hex-2-anhydrobridge-1-procA"] <- "[M-SO3+H]+"
pos_control_peaks_matched$id[pos_control_peaks_matched$id 
                              == "hex-2-anhydrobridge-1-procA"] <- "k-carrageenan DP2"    

#make name + ion column
pos_control_peaks_matched$id_ion <- paste0(pos_control_peaks_matched$id,
                                           ":",
                                           pos_control_peaks_matched$ion,
                                           " mz=",
                                           sprintf("%.3f",
                                                   pos_control_peaks_matched$i.mz))

#order by retention time
pos_control_peaks_matched <- pos_control_peaks_matched[order(rt_round),] 

#write to file
fwrite(pos_control_peaks_matched, 
       file = "./analysis/analysis_tables/pos_control_peaklist_matched.txt",
       sep = "\t")

#14: Extract and format eic -----
#get vectors
mz.found.vector <- pos_control_peaks_matched$i.mz %>% 
    round(., 3)

ions.found.vector  <- pos_control_peaks_matched$id_ion

    
#remove solvent blanks
data2 <- filterFile(data_peaks_grouped_filled,
                    file = which(grepl("standard|acetone|procainamide",
                                       data_peaks_grouped_filled$name)))

#get phenodata vectors
control.names <- data2$name
control.groups <- data2$sample_type
control.groups <- control.groups %>% 
    sub("neg","negative controls", .) %>% 
    sub("pos","oligosaccharide standards mix",.)

#extract chromatograms
chr_list <- list()
error = 0.001

for (i in 1:length(mz.found.vector)){
    mzr = c(mz.found.vector[i] - error,
            mz.found.vector[i] + error)
    chr_list[[i]] <- chromatogram(data2, 
                                  mz = mzr)
}

#extract intensity and rt values
chr_int_list <- list()
for (i in 1:length(control.names)){
    chr_int_list[[i]] <- lapply(chr_list, function(x) {
        x[[i]]@intensity
    }) 
}

chr_rt_list <- list()
for (i in 1:length(control.names)){
    chr_rt_list[[i]] <- lapply(chr_list, function(x) {
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
        rt = chr_rt_list[[i]][[j]]/60
        intensity = chr_int_list[[i]][[j]]
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

control.df$group <- factor(control.df$group,
                           levels = c("oligosaccharide standards mix",
                                      "negative controls"))

#15: Import FLR and format-----------
#read in text files
fld_fp <- dir(path = "./fld-files", 
              all.files = FALSE, 
              full.names = TRUE)
fld.df <- fread(fld_fp)
fld.df$V3 <- NULL
names(fld.df) <- c("rt",
                    "intensity")

fld.df$group <- "oligosaccharide standards mix"

#get FLD minimum and offset so that minimum is zero
#rt range at which minimum is found based on initial plots
fld.df.notzeroed <- fld.df #keep to be safe
int.min = fld.df.notzeroed %>% 
    filter(between(rt, 5, 20)) %>% 
    select(intensity) %>% 
    min()
fld.df$intensity <- fld.df$intensity + abs(int.min)

#transform FLD on x axis - done manually to fit peaks
fld.df$rt_trans <- fld.df$rt - 1.15

#16: Plot FLR----

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
#if you just change plot limits it scales the axis to the procaiamide peak
fld.df.zoom <- fld.df %>% 
    filter(between(rt_trans, 5, 25))

#plot (positive control only)
flr <- ggplot() +
    geom_line(mapping = aes(rt_trans,
                            intensity),
              data = fld.df.zoom,
              colour = "black",
              lwd = 1.2) +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    scale_x_continuous(breaks = seq(5, 25, 1),
                       labels = major_breaks_zoom,
                       limits = c(5, 25),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    theme_classic() +
    theme(axis.text = element_text(size = 12,
                                   family = "Avenir"),
          axis.title = element_text(size = 14, 
                                    family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          legend.position = "none")

#17: Plot EIC ----
#a: set up data ----
#get the most abundant ion for each sugar
pos_control_peaks_matched <- pos_control_peaks_matched[order(id, 
                                                             standard.mix,
                                                             decreasing = TRUE),]
pos_control_abundant_ions <- unique(pos_control_peaks_matched, 
                                    by = "id")
pos_control_abundant_ions <- pos_control_abundant_ions$id_ion
control.df_abundant_ions <- control.df[control.df$ion %in% 
                                           pos_control_abundant_ions,]

#change mz to m/z
control.df_abundant_ions$ion <- control.df_abundant_ions$ion %>% 
    sub("mz", "m/z", .)

#set factor levels
control.df_abundant_ions$ion <- factor(control.df_abundant_ions$ion,
                                       levels = unique(control.df_abundant_ions$ion))

#create palette
pal <- viridis(n = 6)

#make legend labels
ions.for.labels <- levels(control.df_abundant_ions$ion)
ions.for.labels <- ions.for.labels %>% 
    sub(":\\[", ": [", .)
ions.for.labels <- ions.for.labels %>% 
    sub("k-car", "-car", .) %>% 
    sub("b-man", "-man", .)
labels.split <- strsplit(ions.for.labels, "=|\\+ ")
labels.split1 <- unlist(lapply(labels.split, `[[`, 1))
labels.split3 <- unlist(lapply(labels.split, `[[`, 3))

#b: plot positive control ----
pcontrol_eic <- ggplot() +
    geom_line(mapping = aes(rt,
                            intensity,
                            colour = ion),
              data = control.df_abundant_ions[control.df_abundant_ions$group== 
                                                  "oligosaccharide standards mix",],
              lwd = 1.2) +
    scale_colour_manual(values = pal,
                        labels = c(bquote(paste(kappa,
                                                .(labels.split1[1]) ^"+", 
                                                " ", 
                                                italic(m/z),
                                                "=",
                                                .(labels.split3[1]))),
                                   bquote(paste(.(labels.split1[2]) ^"+", 
                                                " ", 
                                                italic(m/z),
                                                "=",
                                                .(labels.split3[2]))),
                                   bquote(paste(kappa,
                                                .(labels.split1[3]) ^"+", 
                                                " ", 
                                                italic(m/z),
                                                "=",
                                                .(labels.split3[3]))),
                                   bquote(paste(.(labels.split1[4]) ^"+", 
                                                " ", 
                                                italic(m/z),
                                                "=",
                                                .(labels.split3[4]))),
                                   bquote(paste(beta,
                                                .(labels.split1[5]) ^"+", 
                                                " ", 
                                                italic(m/z),
                                                "=",
                                                .(labels.split3[5]))),
                                   bquote(paste(.(labels.split1[6]) ^"+", 
                                                " ", 
                                                italic(m/z),
                                                "=",
                                                .(labels.split3[6])))),
                        name = "") +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    scale_x_continuous(breaks = seq(5, 25, 1),
                       labels = major_breaks_zoom,
                       limits = c(5,25),
                       expand = c(0,0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                       labels = scales::scientific,
                       n.breaks = 3) +
    guides(colour=guide_legend(ncol=2)) +
    theme_classic() +
    theme(axis.text = element_text(size = 12,
                                   family = "Avenir"),
          axis.title = element_text(size = 14, 
                                    family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          legend.text = element_text(size = 12, 
                                     family = "Avenir"),
          axis.line = element_blank(),
          legend.position = "bottom")

#c: plot negative controls ----
#get y axis limits of positive control plot
pcontrol_ylim <- ggplot_build(pcontrol_eic)$layout$panel_params[[1]]$y.range

#plot negative controls
ncontrol_eic <- ggplot() +
    geom_line(mapping = aes(rt,
                            intensity,
                            colour = ion),
              data = control.df_abundant_ions[control.df_abundant_ions$group== 
                                                  "negative controls",],
              lwd = 1.2) +
    scale_colour_manual(values = pal,
                        name = "") +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    scale_x_continuous(breaks = seq(5, 25, 1),
                       labels = major_breaks_zoom,
                       limits = c(5,25),
                       expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0),
                       labels = scales::scientific,
                       limits = pcontrol_ylim) +
    guides(colour=guide_legend(ncol=2)) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          plot.background = element_blank(),
          legend.position = "none")

#18: Plot all plots together -----
#add neg control as inset of pos control
eic <- pcontrol_eic + inset_element(ncontrol_eic,
                             0.7, 0.7, 1, 1,
                             align_to = 'full')

tiff("./analysis/analysis_plots/flr_eic_plot_v1.tiff", 
     res = 300, 
     height = 5, 
     width = 12, 
     units = "in")
flr + eic + 
    plot_layout(ncol=1) + 
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 16,
                                  family = "Avenir Next",
                                  face = "bold"))
dev.off()

svg("./analysis/analysis_plots/flr_eic_plot_v1.svg", 
     height = 5, 
     width = 12)
flr + eic + 
    plot_layout(ncol=1) + 
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 16,
                                  family = "Avenir Next",
                                  face = "bold"))
dev.off()


#19: Extract MS2 associated with annotated features in final peak list-----
#get peak info
cp <- chromPeaks(data_ms2)

#get final feature list
peaksFiltered <- pos_control_peaks_matched
peaksFiltered$mz <- NULL
names(peaksFiltered)[names(peaksFiltered) == "i.mz"] <- "mz"
names(peaksFiltered)[names(peaksFiltered) == "i.mzmin"] <- "mzmin"
names(peaksFiltered)[names(peaksFiltered) == "i.mzmax"] <- "mzmax"
peaksFiltered <- peaksFiltered[,c("mz", "mzmin", "mzmax")]


#filter peaks to match identified features
cp <- as.data.frame(cp)
setDT(cp)
setDT(peaksFiltered)
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
cp.filtered$sample <- 1 #only one sample has ms2
cp.new <- as.matrix(cp.filtered)

#assign chromPeaks tp ms2 file
data_ms2_old <- data_ms2 #keep to be safe
data_ms2 <- filterFile(data_ms2, 3)
chromPeaks(data_ms2) <- cp.new

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

precursors.df <- data.frame(peakId = names(precursors),
                            precursorMz = precursors,
                            rtime = rtime,
                            i.mzmin = precursors - 0.0025,
                            i.mzmax = precursors + 0.0025)

setDT(precursors.df)
setDT(pos_control_peaks_matched)
setkey(precursors.df, i.mzmin, i.mzmax)
sample_peaks_matched_ms2 <- foverlaps(pos_control_peaks_matched,
                                      precursors.df)
sample_peaks_matched_ms2 <- sample_peaks_matched_ms2 %>% 
    replace_na(list("peakId"="",
                    "precursorMz" = "", 
                    "rtime"= "", 
                    "mzmin" = "", 
                    "mzmax" = ""))
sample_peaks_matched_ms2 <- sample_peaks_matched_ms2[
    sample_peaks_matched_ms2$peakId!=""]


#extract data for plotting
ms2.df <- data.frame(precursorMz = as.numeric(),
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
    temp <- data.frame(precursorMz = precursorMz,
                       rt = rt,
                       mz = mz,
                       intensity = intensity)
    ms2.df <- rbind(temp,
                    ms2.df)
}

ms2.df$precursorMz <- as.numeric(ms2.df$precursorMz)
ms2.df$rt <- as.numeric(ms2.df$rt)
ms2.df$mz <- as.numeric(ms2.df$mz)


#20: Plot MS2 -----

precursors.vec <- as.numeric(sprintf("%.4f", precursors))
rtime.vec <- as.numeric(sprintf("%.4f", rtime/60))



for (i in 1:length(ms2_features_comb)){
    #get precursor and retention time
    precursorMz = precursors.vec[i]
    rt = rtime.vec[i]
    
    #get x-axis limits
    xlim.max <- round_any(precursorMz, 50, ceiling)
    xlim.min <-  ms2.df %>% 
        filter(precursorMz == !! precursorMz &
                   rt == !! rt &
                   intensity > 1) %>% 
        select(mz) %>% 
        min() %>% 
        round_any(., 50, floor)
    
    #make filename
    filename = paste0("./analysis/ms2_plots/",
                      precursors.vec[i],
                      "mz_",
                      round(rtime.vec[i], 1),
                      "min_ms2")
    
    #plot
    p <- ggplot(ms2.df %>% 
               filter(precursorMz == !! precursorMz &
                          rt == !! rt &
                          intensity > 1),
           aes(x = mz, 
               y = intensity)) +
        geom_segment(aes(x=mz, 
                         xend=mz, 
                         y=0, 
                         yend=intensity),
                     lwd = 1.3) +
        geom_text(aes(x = mz,
                      y = intensity+13,
                      label = sprintf("%.4f",mz)),
                  angle = 90,
                  size = 5,
                  family = "Avenir") +
        geom_point(aes(x=precursorMz, 
                       y = 2),
                   shape = 25,
                   size = 3,
                   fill = "black") +
        scale_x_continuous(name = expression(italic(m/z)),
                           limits = c(xlim.min,
                                      xlim.max),
                           expand = c(0, 0)) +
        scale_y_continuous(name = "Relative intensity (%)",
                           expand = c(0, 0),
                           breaks = seq(0, 100, by = 20),
                           limits = c(0, 130)) +
        theme_classic() +
        theme(axis.text = element_text(size = 20,
                                       family = "Avenir"),
              axis.title = element_text(size = 22,
                                        family = "Avenir LT 65 Medium"),
              panel.border = element_rect(colour = "black",
                                          size = 0.5,
                                          fill = NA),
              axis.line = element_blank())
    
    #save to tiff
    ggsave(filename = paste0(filename, ".tiff"),
           plot = p,
           width = 6,
           height = 6)
    
    #save to svg
    ggsave(filename = paste0(filename, ".svg"),
           plot = p,
           width = 6,
           height = 6)
}


#21: Compute differences within fragmentation spectra associated with features ----
ms2_differences.df <- data.frame(precursorMz = as.numeric(),
                                 differences = as.numeric())

for (i in 1:length(ms2_features_comb)){
    precursorMz.var = as.numeric(sprintf("%.4f",
                              ms2_features_comb[[i]]@precursorMz))
    x <- as.numeric(sprintf("%.4f",ms2_features_comb[[i]]@mz))
    y <- data.frame(difference = c(dist(x)))
    y <- cbind(precursorMz = rep(precursorMz.var,
                                 nrow(y)),
               y)
    ms2_differences.df <- rbind(y,
                                ms2_differences.df) 
}

#significant differences:
    #73.089 -> R1 procainamide fragmentation
    #116.131 -> R2 procainamide fragmentation
    #43.0422 -> difference between R1 and R2 fragments
    #162.053 -> loss of hexose monomer
    #79.957 -> loss of sulphate


sig.ms2_differences.df <- ms2_differences.df %>% 
    filter(between(difference, 73.088, 73.090) |
               between(difference, 116.130, 116.132) |
               between(difference, 43.041, 43.043) |
               between(difference, 162.052, 162.054) |
               between(difference, 79.956, 79.958)) 

#format table
sig.ms2_differences.df$difference_rounded <- round(sig.ms2_differences.df$difference)
sig.ms2_differences.df_wide <- sig.ms2_differences.df %>% 
    pivot_wider(names_from = difference_rounded,
                values_from = difference,
                values_fn = list)
names(sig.ms2_differences.df_wide)[names(sig.ms2_differences.df_wide) == "43"] <- "R1.R2.diff"
names(sig.ms2_differences.df_wide)[names(sig.ms2_differences.df_wide) == "116"] <- "R2"
names(sig.ms2_differences.df_wide)[names(sig.ms2_differences.df_wide) == "162"] <- "hexose.monomer"
names(sig.ms2_differences.df_wide)[names(sig.ms2_differences.df_wide) == "73"] <- "R1"


#21: Compute differences within fragmentation spectra associated with features ----
ms2_differences.df <- data.frame(precursorMz = as.numeric(),
                                 differences = as.numeric())

for (i in 1:length(ms2_features_comb)){
    precursorMz.var = as.numeric(sprintf("%.4f",
                                         ms2_features_comb[[i]]@precursorMz))
    x <- as.numeric(sprintf("%.4f",ms2_features_comb[[i]]@mz))
    y <- data.frame(difference = c(dist(x)))
    y <- cbind(precursorMz = rep(precursorMz.var,
                                 nrow(y)),
               y)
    ms2_differences.df <- rbind(y,
                                ms2_differences.df) 
}

#significant differences:
#73.089 -> R1 procainamide fragmentation
#116.131 -> R2 procainamide fragmentation
#43.0422 -> difference between R1 and R2 fragments
#162.053 -> loss of hexose monomer
#79.957 -> loss of sulphate


sig.ms2_differences.df <- ms2_differences.df %>% 
    filter(between(difference, 73.088, 73.090) |
               between(difference, 116.130, 116.132) |
               between(difference, 43.041, 43.043) |
               between(difference, 162.052, 162.054) |
               between(difference, 79.956, 79.958)) 

#format table
sig.ms2_differences.df$difference_rounded <- round(sig.ms2_differences.df$difference)
sig.ms2_differences.df_wide <- sig.ms2_differences.df %>% 
    pivot_wider(names_from = difference_rounded,
                values_from = difference,
                values_fn = list)
names(sig.ms2_differences.df_wide)[names(sig.ms2_differences.df_wide) == "43"] <- "R1.R2.diff"
names(sig.ms2_differences.df_wide)[names(sig.ms2_differences.df_wide) == "116"] <- "R2"
names(sig.ms2_differences.df_wide)[names(sig.ms2_differences.df_wide) == "162"] <- "hexose.monomer"
names(sig.ms2_differences.df_wide)[names(sig.ms2_differences.df_wide) == "73"] <- "R1"


#22: Screen all MS2 for significant differences ----

#ideas from CAMERA manual and tutorials:
####
#bin-methods
    #aggregates individual spectra (Spectrum instances) or whole experiments 
    #(MSnExp instances) into discrete bins
#calculateFragments-methods
    #These method calculates a-, b-, c-, x-, y- and z-ions produced by fragmentation.

#compareSpectra-methods
    #This method compares spectra (Spectrum instances) pairwise or all spectra of an experiment (MSnExp
# instances). Currently the comparison is based on the number of common peaks fun = "common",
# the Pearson correlation fun = "cor", the dot product fun = "dotproduct" or a user-defined
# function.

#extractPrecSpectra-methods
    #Extracts the MSMS spectra that originate from the precursor(s) having the same MZ value as defined
    #in theprec argument.


#DIA ANALYSIS
#table(isolationWindowTargetMz(data_ms2))
    #https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms-lcms-ms.html#2_Analysis_of_DDA_data
    #reconstructChromPeakSpectra
    #findChromPeaksIsolationWindow




