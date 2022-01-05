setwd("/Users/margotbligh/Google_Drive/MPI_PhD/Lab-things/alpha-mannan/Polaribacter_Hel1-33-78_enzymes/202108_GH99_PL29")
load("./analysis/RData/RData_20210818.RData")
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
fp <- dir(path = "data/mzML-files/MS31_20210813", 
          all.files = FALSE, 
          full.names = TRUE)

#create phenodata data.frame
#each sample must have a unique name!
pd <- data.frame(name = basename(fp) %>%
                     sub("MS31_20210813_", "", .) %>% 
                     sub("Std_", "", .) %>% 
                     sub("_\\d\\d.mzML", "", .),
                 sample_type = basename(fp) %>% 
                     sub(".*PL29_peak.*", "PL29 purified oligo", .) %>% 
                     sub(".*GH99_peak.*", "GH99 purified oligo", .) %>%
                     sub(".*neg_.*", "negative control", .) %>% 
                     sub(".*GH99_[ABC]_.*", "GH99 digest", .) %>% 
                     sub(".*PL29_[ABC]_.*", "PL29 digest", .) %>% 
                     sub(".*Std.*", "standard", .) %>% 
                     sub(".*HILIC.*", "HILIC standard", .) %>% 
                     sub(".*Solvent.*", "solvent blank", .),
                 enzyme = basename(fp) %>% 
                     sub(".*PL29.*", "PL29", .) %>% 
                     sub(".*GH99.*", "GH99", .) %>% 
                     sub(".*Std.*|.*HILIC.*|.*Solvent.*", "", .),
                 injection = basename(fp) %>% 
                     sub(".*[ABC]_|.*g_|.*d_|.*k\\d_|.*e_", "", .) %>% 
                     sub(".mzML", "", .) %>% as.numeric(),
                 stringsAsFactors = FALSE)

pd$name[pd$sample_type == "standard"] <- paste0(pd$name[pd$sample_type == "standard"],
                                                "_",
                                                pd$injection[pd$sample_type == "standard"])

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

pal_group <- hcl.colors(n = length(unique(pd$sample_type)),
                        palette = "Dark3")
names(pal_group) <- unique(pd$sample_type)

#plot the tic as boxplot
tc <- split(tic(all_data), 
            f = fromFile(all_data))
cairo_pdf("./analysis/processing_plots/tic/tic_boxplot.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mar=c(9,5,1,1))
boxplot(tc, 
        col = pal_group[all_data$sample_type],
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
    group = rep(all_data$sample_type[i], length(rt))
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
cwp@ppm<-1.2
cwp@peakwidth<-c(7,45)
cwp@snthresh<-20
cwp@noise <- 5000
cwp@prefilter <- c(3, 1000)


#check with chromatograms for standard (sulphated mannoase)
dir.create("./analysis/processing_plots/peakpicking",
           showWarnings = FALSE)

error = 0.0025
chr1 <- chromatogram(data,
                     rt = c(50, 400),
                     mz = c(259.01293 - error,
                            259.01293 + error))
chr_cwp1 <- findChromPeaks(chr1, 
                           param=cwp)
pal1 <- hcl.colors(n = length(pd$name),
                   palette = "Dark3")
cairo_pdf("./analysis/processing_plots/peakpicking/peakpicking_1.3ppm_20sn_pw7to45_noise5000_prefilter3-1000_sulphatedmannose.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
plot(chr_cwp1,
     lwd = 2,
     cex.main = 1,
     peakCol = pal1[chromPeaks(chr_cwp1)[, "column"]],
     peakType = "rectangle")
dev.off()

#pick peaks
data_peaks<-findChromPeaks(data, 
                           param=cwp)
data_ms2<-findChromPeaks(data_ms2, 
                         param=cwp)
#save as RData objects
save(data_peaks, 
     file = "./analysis/RData/data_peaks.RData")
save(data_ms2, 
     file = "./analysis/RData/data_ms2.RData")

#6: Group peaks to create "features"---------
#parameters
pdp <- PeakDensityParam(sampleGroups = data$sample_type,
                        binSize = 0.005,
                        bw = 6,
                        minFraction = 0.2) 

#check parameters
#extract and plot chromatograms to test settings
dir.create("./analysis/processing_plots/peakgrouping",
           showWarnings = FALSE)

chr_pdp1 <- chromatogram(data_peaks,
                     rt = c(50, 400),
                     mz = c(259.01293 - error,
                            259.01293 + error))
names(pal1) <- data$name #name palette

cairo_pdf("./analysis/processing_plots/peakgrouping/peakgrouping_binsize0.005_bw6_sulphatedmannose.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mar=c(4,4,3,10))
plotChromPeakDensity(chr_pdp1, 
                     col = pal1, 
                     param = pdp,
                     peakBg = pal1[chromPeaks(chr_pdp1)[, "sample"]],
                     peakCol = pal1[chromPeaks(chr_pdp1)[, "sample"]],
                     peakPch = 16)
legend("topright",
       legend = paste0(seq(1,length(pal1),1),
                       "=",
                       names(pal1)),
       inset=c(-0.6,0),
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

res <- data_peaks_grouped_filled

#8: Remove large data files from environment----------
rm(data, data_peaks, data_peaks_grouped, data_peaks_grouped_filled)    

#9: Save diffreport of xdata -----
xset <- as(res, "xcmsSet")
sampnames(xset) <- pData(res)$name
sampclass(xset) <- pData(res)$sample_type

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
mz_predicted <- fread("predicted-masses-dp1to7-sulphatedhexose-unsaturated.txt")

mz_predicted <- fread("predicted-masses-dp1to7-sulphatedhexose-unsaturated-carboxyl-multimod.txt")


#remove "extra" columns
extraCol <- c('mass',
              'formula')

mz_predicted <- mz_predicted %>% 
    select(-all_of(extraCol))

#make long format
predicted <- mz_predicted %>% 
    gather(key = "ion",
           value = "mz",
           -name,
           -dp)

predicted <- predicted[grep("unsaturated-[234567]" , predicted$name, invert = TRUE),]

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







