#setwd("C:Users/mbligh/ownCloud/marglyco/Data/LC-MS_data/FITDOG/20211216")

setwd("~/ownCloud/marglyco/Data/LC-MS_data/FITDOG/20211216")
load("analysis/RData/RData_20220114.RData")
save.image("analysis/RData/RData_20220114.RData")

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
fp <- dir(path = "mzML", all.files = FALSE, full.names = TRUE)

#create phenodata data.frame
#each sample must have a unique name!
pd <- data.frame(name = basename(fp) %>%
                     sub("MS78_20211216_", "", .) %>% 
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
          family = "Avenir",width = 12, height = 9)
par(mar=c(9,5,1,1))
boxplot(tc, col = pal_group[all_data$source], ylab = "intensity", 
        main = "total ion current", names = all_data$name, las=2,cex.axis = 0.8)
dev.off()

#plot as a chromatogram
tic <- chromatogram(all_data, aggregationFun = "sum")
tic.df <- data.frame(sample = as.character(), group = as.character(),
                     rt = as.numeric(),intensity = as.numeric())
for (i in 1:length(pd$name)){
    rt = tic[[i]]@rtime / 60
    intensity = tic[[i]]@intensity
    sample = rep(all_data$name[i], length(rt))
    group = rep(all_data$source[i], length(rt))
    temp <- data.frame(sample = sample, group = group, rt = rt,
                       intensity = intensity)
    tic.df <- rbind(temp, tic.df)
}

cairo_pdf("./analysis/processing_plots/tic/tic_chromatogram.pdf",
          family = "Avenir", width = 12, height = 9)
ggplot() +
    geom_line(mapping = aes(rt, intensity, colour = group, group = sample),
              data = tic.df, lwd = 1.2) +
    scale_colour_manual(values = pal_group) +
    labs(x= "Retention time (min)", y = "Intensity (a.u.)") +
    facet_grid(rows = vars(group), scales = "free_y") +
    theme_classic() +
    theme(strip.text = element_blank(),
          axis.text = element_text(size = 12,family = "Avenir"),
          axis.title = element_text(size = 14,family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          legend.text = element_text(size = 6, family = "Avenir"),
          axis.line = element_blank())
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

#check chromatograms
n = pd$source %>% unique() %>% length()
sources <- pd$source %>% unique()


chr_hex1 <- chromatogram(data_peaks, mz = c(215.03279-0.001,215.03279+0.001))
par(mfrow=c(n/2,2))
for (i in 1:n){
    plot(chr_hex1[,which(pd$source == sources[i])],
         main = sources[i])
}
chr_hex2 <- chromatogram(data_peaks, mz = c(377.08561-0.001,377.08561+0.001))
par(mfrow=c(n/2,2))
for (i in 1:n){
    plot(chr_hex2[,which(pd$source == sources[i])],
         main = sources[i])
}
chr_hex3 <- chromatogram(data_peaks, mz = c(539.1384-0.001,539.1384+0.001))
par(mfrow=c(n/2,2))
for (i in 1:n){
    plot(chr_hex3[,which(pd$source == sources[i])],
         main = sources[i])
}

hex2_deox2_s1_chr <- chromatogram(res, mz = c(389.0760-0.001, 389.0760+0.001))
par(mfrow=c(n/2,2))
for (i in 1:n){
    plot(hex2_deox2_s1_chr[,which(pd$source == sources[i])],
         main = sources[i])
}
hex2_deox2_chr <- chromatogram(res, mz = c(345.0958-0.001, 345.0958+0.001))
par(mfrow=c(n/2,2))
for (i in 1:n){
    plot(hex2_deox2_chr[,which(pd$source == sources[i])],
         main = sources[i])
}
hex3_deox3_s3_chr <- chromatogram(res, mz = c(231.0110-0.001, 231.0110+0.001))
par(mfrow=c(n/2,2))
for (i in 1:n){
    plot(hex3_deox3_s3_chr[,which(pd$source == sources[i])],
         main = sources[i])
}
hex3_deox3_s2_chr <- chromatogram(res, mz = c(307.0417-0.001, 307.0417+0.001))
par(mfrow=c(n/2,2))
for (i in 1:n){
    plot(hex3_deox3_s2_chr[,which(pd$source == sources[i])],
         main = sources[i])
}
hex1_deox1_chr <- chromatogram(res, mz = c(199.0379-0.001, 199.0379+0.001))
par(mfrow=c(n/2,2))
for (i in 1:n){
    plot(hex1_deox1_chr[,which(pd$source == sources[i])],
         main = sources[i])
}





#save as RData objects
save(data_peaks, file = "./analysis/RData/data_peaks.RData")
save(data_ms2, file = "./analysis/RData/data_ms2.RData")

#6: Group peaks to create "features"---------
#parameters
pdp <- PeakDensityParam(sampleGroups = data$source,
                        binSize = 0.005, bw = 6, minFraction = 0.2) 

#check parameters
#extract and plot chromatograms to test settings
dir.create("./analysis/processing_plots/peakgrouping", showWarnings = FALSE)

chr_pdp1 <- chromatogram(data_peaks, rt = c(100, 1000),  
                         mz = c(377.08561 - 0.001, 377.08561 + 0.001))

pal_sample.by.group <- pal_group[pd$source]
names(pal_sample.by.group) <- pd$name

cairo_pdf("./analysis/processing_plots/peakgrouping/peakgrouping_binsize0.005_bw6_dihexose.pdf",
          family = "Avenir", width = 12, height = 9)
par(mar=c(4,4,3,10))
plotChromPeakDensity(chr_pdp1, col = pal_sample.by.group, param = pdp,
                     peakBg = pal_sample.by.group[chromPeaks(chr_pdp1)[, "sample"]],
                     peakCol = pal_sample.by.group[chromPeaks(chr_pdp1)[, "sample"]],
                     peakPch = 16)
legend("topright", legend = paste0(seq(1,length(pal_sample.by.group),1), "=", 
                                   names(pal_sample.by.group)),
       inset=c(0,0), fill = pal_sample.by.group, pt.cex = 0.3,cex = 0.5, 
       bty = "n",horiz = FALSE, xpd=TRUE, ncol = 3)
dev.off()


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

#make rownames from rt and mz of features
rownames(pl_be)<-paste(round(pl_be$rt,1), round(pl_be$mz,3), sep="_")
#change NA to 0
pl_be[is.na(pl_be)] <- 0

#change name
peaks <- pl_be

#add rounded retention time as first colum
peaks <- cbind(rt_min = round(peaks$rt/60, 1), peaks)

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
peaks_iso <- unique(peaks_iso, by = "isotope_group")
#merge to get concatenated isotope lists
peaks_iso <- merge(peaks_iso, iso_concat, by = "isotope_group")
#clean up df
peaks_iso <- peaks_iso %>% 
    select(-c("isotope_group", "isotope_number","isotopes.x"))
names(peaks_iso)[names(peaks_iso) == 'isotopes.y'] <- 'isotopes'

#merge features with and without isotopes
peaks <- rbind.fill(peaks_noiso, peaks_iso)

#13: Annotate features based on predictions----
#predict sugars - function executes python script
predictSugars <- function(dp1,dp2,ESI_mode, scan_range1, scan_range2,
                          pent_option=NULL,modifications=NULL,
                          nmod_max=NULL,double_sulphate=NULL,label=NULL){
    py_install("pandas", "numpy")
    source_python("../../sugarMassesPredict-r.py")
    dp1 = as.integer(dp1)
    dp2 = as.integer(dp2)
    scan_range1 = as.integer(scan_range1)
    scan_range2 = as.integer(scan_range2)
    if(!is.null(pent_option)){
        pent_option = as.integer(pent_option)
    } else{
        pent_option = as.integer(0)
    }
    if(!is.null(nmod_max)){
        nmod_max = as.integer(nmod_max)
    } else{
        nmod_max = as.integer(1)
    }
    if(!is.null(double_sulphate)){
        double_sulphate = as.integer(double_sulphate)
    } else{
        double_sulphate = as.integer(0)
    }
    if(!is.null(modifications)){
        if(is.vector(modifications)){
            modifications = as.list(modifications)
        }
    } else{
        modifications = "none"
    }
    if(is.null(label)){
        label = "none"
    }
    df <- predict_sugars(dp1 = dp1, dp2 = dp2, ESI_mode = ESI_mode,
                         scan_range1 = scan_range1, scan_range2 = scan_range2,
                         pent_option = pent_option, modifications = modifications,
                         label = label, nmod_max = nmod_max, 
                         double_sulphate = double_sulphate)
    return(df)
}
mz_predicted <- predictSugars(dp1 = 1, dp2 = 8, ESI_mode = 'neg', scan_range1 = 175, 
                           scan_range2 = 1400, pent_option = 1, 
                           modifications = c("sulphate", "deoxy", "carboxyl"),
                           nmod_max = 2, double_sulphate = 1)

#remove "extra" columns
extraCol <- c('mass','formula')

mz_predicted <- mz_predicted %>% select(-all_of(extraCol))

#make long format
predicted <- mz_predicted %>% 
    gather(key = "ion", value = "mz", -name, -dp)

#make data.table
setDT(predicted); setDT(peaks)
#create interval to overlap with (same width as for peak grouping)
predicted$mz <- as.numeric(predicted$mz)
predicted$mzmin <- predicted$mz-0.001
predicted$mzmax <- predicted$mz+0.001

#remove NA rows
predicted <- na.omit(predicted)

#match using foverlaps from data.table (very fast)
setkey(predicted, mzmin, mzmax)
peaks_nopred <- peaks
peaks <- foverlaps(peaks,predicted)

#change NA values created during matching (features with no match) to be blank
#remove extra columns
peaks$mzmin <- NULL
peaks$mzmax <- NULL

peaks <- peaks %>% 
    replace_na(list("name"="", "ion"= "", "mz" = "","dp" = ""))
#only keep matched features
peaks_matched <- peaks[!peaks$name=="",]

#remove annotations with sulphate or carboxyl and M+Cl or M+CHOO
a <- grep("sulphate|carboxyl", peaks_matched$name)
b <- grep("Cl|CHOO", peaks_matched$ion)
c <- intersect(a, b)
peaks_matched <- peaks_matched[-c,]

#order by retention time
peaks_matched <- peaks_matched[order(rt_min),]

#remove everything before 1 min and after 22 min
peaks_matched <- peaks_matched %>% filter(between(rt_min, 1, 22))

#make id and ion column
peaks_matched$id_ion_individual <- paste0(peaks_matched$name, ": ", peaks_matched$ion)

#aggregate so that if there are multiple predictions for one feature
#they are shown in the same row. delete all of the other extra columns added 
#during matching
names <- setdiff(names(peaks_matched), names(predicted))
names <- names[!names == "id_ion_individual"]
setDF(peaks_matched)
peaks_matched_uncollapsed <- peaks_matched
peaks_matched <- peaks_matched_uncollapsed %>% 
    dplyr::group_by(across(all_of(names))) %>% 
    dplyr::mutate(id_ion = paste0(id_ion_individual, collapse = ", ")) %>% 
    ungroup() %>% 
    distinct(across(all_of(c(names, "id_ion"))))
names(peaks_matched)[names(peaks_matched) == 'i.mz'] <- 'mz'
names(peaks_matched)[names(peaks_matched) == 'i.mzmin'] <- 'mzmin'
names(peaks_matched)[names(peaks_matched) == 'i.mzmax'] <- 'mzmax'



peaks_matched_old <- peaks_matched

colOrder <- names(peaks_matched)
colOrder <- colOrder[1:length(colOrder)-1]
colOrder <- c("id_ion", colOrder)

setDF(peaks_matched)
peaks_matched <- peaks_matched[,colOrder]

#14: Extract MS2 associated with annotated features in final peak list-----
#get peak info
cp <- chromPeaks(data_ms2)

#format final feature list
peaksFiltered <- peaks_matched
peaksFiltered <- peaksFiltered[,c("mz", "mzmin", "mzmax")]

#filter peaks to match identified features
cp <- as.data.frame(cp)
setDT(cp); setDT(peaksFiltered)
setkey(cp, mzmin, mzmax)
setkey(peaksFiltered, mzmin, mzmax)
cp.matched <- foverlaps(cp,peaksFiltered)
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

#extract ms2 associated with chromPeaks
ms2_features <- chromPeakSpectra(data_ms2, msLevel = 2, expandRt = 2.5,
                                 expandMz = 0.005, skipFilled = FALSE,
                                 method = "all")

#remove zero intensity masses
ms2_features@listData <- lapply(ms2_features, clean,  all = TRUE)

#combine spectra by rtime within each file
mcols(ms2_features) <- c(mcols(ms2_features),
                         "file_rt_precursor" = paste0(fromFile(ms2_features), "_",
                                                      round(rtime(ms2_features)/60, 1), "_",
                                                      round(precursorMz(ms2_features),3)))
names(mcols(ms2_features)) <- c("peak_id", "file_rt_precursor")
ms2_features_comb <- combineSpectra(ms2_features, fcol = 'file_rt_precursor',
                                    mzd = 0.005, intensityFun = mean)

#see how many spectra per sample
table(fromFile(ms2_features_comb))

#normalise with respect to ion with highest intensity
ms2_features_comb_old <- ms2_features_comb
ms2_features_comb@listData <- lapply(ms2_features_comb_old, normalise, 
                                     method = "max")

#see which features have ms2 spectra
precursors <- precursorMz(ms2_features_comb)
rtime <- rtime(ms2_features_comb)
ms2_samples <- fromFile(ms2_features_comb)
names(ms2_samples) <- data_ms2$name[ms2_samples]

precursors.df <- data.frame(peakId = names(precursors), precursorMz = precursors,
                            rtime = rtime, sample = names(ms2_samples),
                            mzmin = precursors - 0.0025,
                            mzmax = precursors + 0.0025)

setDT(precursors.df); setDT(peaks_matched)
setkey(precursors.df, mzmin, mzmax)
sample_peaks_matched_ms2 <- foverlaps(peaks_matched, precursors.df)
sample_peaks_matched_ms2 <- sample_peaks_matched_ms2 %>% 
    replace_na(list("peakId"="","precursorMz" = "", 
                    "rtime"= "","sample" = "",
                    "mzmin" = "", "mzmax" = ""))
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

#15: Look at annotated fingerprints of glycans, EIC for all features in all samples------
par(mfrow=c(1,1))

    #F. vesiculosus fucoidan-----

ves <- filter(peaks_matched, F..vesiculosus >= 2)
ves.mz <- round(ves$mz, 3)
names(ves.mz) <- ves$id_ion

ves.mz.unique <- unique(ves.mz)
names(ves.mz.unique) <- names(ves.mz)[match(ves.mz.unique, ves.mz)]


#extract chromatograms
ves_chr <- list()
for(i in 1:length(ves.mz.unique)){
    ves_chr[[i]] <- chromatogram(res, mz = c(ves.mz.unique[i]-0.001, 
                                             ves.mz.unique[i]+0.001))
}

save(ves_chr, file = "analysis/RData/ves_chr.RData")

ves_chr.df <- data.frame(sample = as.character(),
                         substrate = as.character(),
                         source = as.character(),
                         replicate = as.character(),
                         mz = as.numeric(), rt = as.numeric(), 
                         intensity = as.numeric(),
                         ion = as.character())
for(i in 1:length(ves.mz.unique)){
    chr <- ves_chr[[i]]
    mz <- ves.mz.unique[i]
    annot <- names(ves.mz.unique)[i]
    
    lengths <- unlist(lapply(chr, function(x) length(x)))
    
    sample <- rep(chr$name, times = lengths)
    substrate <- rep(chr$substrate, times = lengths)
    source <- rep(chr$source, times = lengths)
    replicate <- rep(chr$replicate, times = lengths)
    mz <- rep(mz, sum(lengths))
    rt <- unlist(lapply(chr, function(x) x@rtime/60))
    intensity <- unlist(lapply(chr, function(x) x@intensity))
    ion <- rep(annot, sum(lengths))
    
    tmp <- data.frame(sample = sample, substrate = substrate,
                      source = source, replicate = replicate,
                      mz = mz, rt = rt, intensity = intensity,
                      ion = ion) 
    ves_chr.df <- rbind(ves_chr.df, tmp)
}

ves_chr.df$plot.group <- paste0(ves_chr.df$sample, "_", ves_chr.df$mz)

#remove any features for which the max intensity is not >= 1e5

ves_chr.df <- ves_chr.df %>% 
    group_by(mz) %>% 
    filter(any(intensity >= 1e5)) 

#manually check chromatograms (peak picking was not stringent, don't know how
#to get around doing this at the moment)
setDF(ves)
ves <- ves[order(ves$rt_min),]
ves_chr.df$mz <- factor(ves_chr.df$mz, levels = unique(round(ves$mz,3))) 

p <- ggplot(data = ves_chr.df %>% filter(source == "F. vesiculosus"), 
       aes(x = rt, y = intensity, group = replicate,
                              colour = replicate)) +
    geom_line() +
    facet_grid(rows = vars(as.factor(mz)), scales = "free_y") +
    xlim(0,15) +
    scale_y_continuous(labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 360),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/f_ves_v1.png",
    height = 20, width = 10, units = "in", res = 300)
print(p)
dev.off()

mz_nope <- c(245.022, 311.005, #assume same compound as 300.976
             565.054, 357.049, 451.022, 473.004, 391.071, 219.007, 302 #not really peaks
             )
ves_chr.df2 <- ves_chr.df %>% filter(!mz %in% mz_nope) 

ves_chr.df2$ion2[ves_chr.df2$mz == 405.07] <- "heoxse-deoxyhexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 387.06] <- "dehydrated heoxse-deoxyhexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 300.976] <- "unknown, m/z = 300.976"
ves_chr.df2$ion2[ves_chr.df2$mz == 259.013] <- "hexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 241.002] <- "dehydrated hexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 389.076] <- "di-deoxyhexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 371.065] <- "dehydrated di-deoxyhexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 214.986] <- "unknown, m/z = 214.986"
ves_chr.df2$ion2[ves_chr.df2$mz == 243.018] <- "deoxyhexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 225.007] <- "dehydrated deoxyhexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 229.002] <- "pentose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 210.991] <- "dehydrated pentose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 535.134] <- "tri-deoxyhexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 322.974] <- "unknown, m/z = 322.974 (sulphated)"
ves_chr.df2$ion2[ves_chr.df2$mz == 345.095] <- "di-deoxyhexose: [M+Cl]-"

p <- ggplot(data = ves_chr.df2 %>% filter(source == "F. vesiculosus")) +
    geom_line(aes(x = rt, y = intensity, group = plot.group,
                  colour = replicate)) +
    facet_grid(rows = vars(ion2), scales = "free") +
    scale_colour_manual(values = pal_3) +
    xlim(0,7) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/f_ves_v2.png",
    height = 6, width = 10, units = "in", res = 300)
print(p)
dev.off()

#plot for all samples
pal_3 <- c("#FEC000", "#2BB6AF", "#D1D1FE")
ves_chr.df2$source <- factor(ves_chr.df2$source, 
                             levels = c("E. bicyclis","L. digitata",
                                        "F. vesiculosus", "L. hyperborea",
                                        "S. cerevisisae", "T. weissflogii",
                                        "solvent blank"))


p <- ggplot(data = ves_chr.df2) +
    geom_line(aes(x = rt, y = intensity, group = plot.group,
                  colour = replicate)) +
    facet_grid(rows = vars(ion2), cols = vars(source), scales = "free") +
    scale_colour_manual(values = pal_3) +
    xlim(0,7) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          strip.text.x = element_text(vjust = 1, size = 10),
          legend.position = "none",
          strip.background = element_rect(colour = NA))

png(filename = "analysis/analysis_plots/fingerprint_chromatograms_all-samples/20210114_f_ves.png",
    height = 7, width = 12, units = "in", res = 300)
print(p)
dev.off()


rm(ves_chr)


    #L. hyperborea fucoidan-----

hyp <- filter(peaks_matched, L..hyperborea >= 2)
hyp.mz <- round(hyp$mz, 3)
names(hyp.mz) <- hyp$id_ion

hyp.mz.unique <- unique(hyp.mz)
names(hyp.mz.unique) <- names(hyp.mz)[match(hyp.mz.unique, hyp.mz)]


#extract chromatograms
hyp_chr <- list()
for(i in 1:length(hyp.mz.unique)){
    hyp_chr[[i]] <- chromatogram(res, mz = c(hyp.mz.unique[i]-0.001, 
                                             hyp.mz.unique[i]+0.001))
}

save(hyp_chr, file = "analysis/RData/hyp_chr.RData")


hyp_chr.df <- data.frame(sample = as.character(),
                         substrate = as.character(),
                         source = as.character(),
                         replicate = as.character(),
                         mz = as.numeric(), rt = as.numeric(), 
                         intensity = as.numeric(),
                         ion = as.character())
for(i in 1:length(hyp.mz.unique)){
    chr <- hyp_chr[[i]]
    mz <- hyp.mz.unique[i]
    annot <- names(hyp.mz.unique)[i]
    
    lengths <- unlist(lapply(chr, function(x) length(x)))
    
    sample <- rep(chr$name, times = lengths)
    substrate <- rep(chr$substrate, times = lengths)
    source <- rep(chr$source, times = lengths)
    replicate <- rep(chr$replicate, times = lengths)
    mz <- rep(mz, sum(lengths))
    rt <- unlist(lapply(chr, function(x) x@rtime/60))
    intensity <- unlist(lapply(chr, function(x) x@intensity))
    ion <- rep(annot, sum(lengths))
    
    tmp <- data.frame(sample = sample, substrate = substrate,
                      source = source, replicate = replicate,
                      mz = mz, rt = rt, intensity = intensity,
                      ion = ion) 
    hyp_chr.df <- rbind(hyp_chr.df, tmp)
}

hyp_chr.df$plot.group <- paste0(hyp_chr.df$sample, "_", hyp_chr.df$mz)

#remove any features for which the max intensity is not >= 1e5
hyp_chr.df <- hyp_chr.df %>% 
    group_by(mz) %>% 
    filter(any(intensity >= 1e5)) 

#manually check chromatograms (peak picking was not stringent, don't know how
#to get around doing this at the moment)
setDF(hyp)
hyp <- hyp[order(hyp$rt_min),]
hyp_chr.df$mz <- factor(hyp_chr.df$mz, levels = unique(round(hyp$mz,3))) 

p <- ggplot(data = hyp_chr.df %>% filter(source == "L. hyperborea"), 
            aes(x = rt, y = intensity, group = replicate,
                colour = replicate)) +
    geom_line() +
    facet_grid(rows = vars(as.factor(mz)), scales = "free_y") +
    xlim(0,15) +
    scale_y_continuous(labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 360),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/l_hyp_v1.png",
    height = 15, width = 8, units = "in", res = 300)
print(p)
dev.off()

mz_nope <- c(245.022, 311.005, #assume same compound as 300.976
             451.022, #assume same compound as 473.004
             601.639 #not really peaks
)

hyp_chr.df2 <- hyp_chr.df %>% filter(!mz %in% mz_nope) 

hyp_chr.df2$ion2[hyp_chr.df2$mz == 243.018] <- "deoxyhexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 300.976] <- "unknown, m/z = 300.976"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 357.049] <- "unknown, m/z = 357.049"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 214.986] <- "unknown, m/z = 214.986"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 387.06] <- "dehydrated heoxse-deoxyhexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 473.004] <- "unknown, m/z = 473.004"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 229.002] <- "pentose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 225.007] <- "dehydrated deoxyhexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 322.974] <- "unknown, m/z = 322.974 (sulphated)"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 389.076] <- "di-deoxyhexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 259.013] <- "hexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 234.012] <- "di-deoxyhexose disulphate: [M-2H]-2"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 412.944] <- "unknown, m/z = 412.944"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 210.991] <- "dehydrated pentose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 302] <- "unknown, m/z = 302.000"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 198.991] <- "unknown, m/z = 198.991 (sulphated)"


p <- ggplot(data = hyp_chr.df2 %>% filter(source == "L. hyperborea")) +
    geom_line(aes(x = rt, y = intensity, group = plot.group,
                  colour = replicate)) +
    facet_grid(rows = vars(ion2), scales = "free") +
    scale_colour_manual(values = pal_3) +
    xlim(0,7) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/l_hyp_v2.png",
    height = 6, width = 10, units = "in", res = 300)
print(p)
dev.off()

#plot for all samples
hyp_chr.df2$source <- factor(hyp_chr.df2$source, 
                             levels = c("E. bicyclis","L. digitata",
                                        "F. vesiculosus", "L. hyperborea",
                                        "S. cerevisisae", "T. weissflogii",
                                        "solvent blank"))


p <-ggplot(data = hyp_chr.df2) +
    geom_line(aes(x = rt, y = intensity, group = plot.group,
                  colour = replicate)) +
    facet_grid(rows = vars(ion2), cols = vars(source), scales = "free") +
    scale_colour_manual(values = pal_3) +
    xlim(0,10) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          strip.text.x = element_text(vjust = 1, size = 10),
          legend.position = "none",
          strip.background = element_rect(colour = NA))

png(filename = "analysis/analysis_plots/fingerprint_chromatograms_all-samples/20210114_l_hyp.png",
    height = 7, width = 12, units = "in", res = 300)
print(p)
dev.off()

rm(hyp_chr)

    #L. digitata laminarin-----

dig <- filter(peaks_matched, L..digitata >= 2)
dig.mz <- round(dig$mz, 3)
names(dig.mz) <- dig$id_ion

dig.mz.unique <- unique(dig.mz)
names(dig.mz.unique) <- names(dig.mz)[match(dig.mz.unique, dig.mz)]


#extract chromatograms
dig_chr <- list()
for(i in 1:length(dig.mz.unique)){
    dig_chr[[i]] <- chromatogram(res, mz = c(dig.mz.unique[i]-0.001, 
                                             dig.mz.unique[i]+0.001))
}

dig_chr.df <- data.frame(sample = as.character(),
                         substrate = as.character(),
                         source = as.character(),
                         replicate = as.character(),
                         mz = as.numeric(), rt = as.numeric(), 
                         intensity = as.numeric(),
                         ion = as.character())
for(i in 1:length(dig.mz.unique)){
    chr <- dig_chr[[i]]
    mz <- dig.mz.unique[i]
    annot <- names(dig.mz.unique)[i]
    
    lengths <- unlist(lapply(chr, function(x) length(x)))
    
    sample <- rep(chr$name, times = lengths)
    substrate <- rep(chr$substrate, times = lengths)
    source <- rep(chr$source, times = lengths)
    replicate <- rep(chr$replicate, times = lengths)
    mz <- rep(mz, sum(lengths))
    rt <- unlist(lapply(chr, function(x) x@rtime/60))
    intensity <- unlist(lapply(chr, function(x) x@intensity))
    ion <- rep(annot, sum(lengths))
    
    tmp <- data.frame(sample = sample, substrate = substrate,
                      source = source, replicate = replicate,
                      mz = mz, rt = rt, intensity = intensity,
                      ion = ion) 
    dig_chr.df <- rbind(dig_chr.df, tmp)
}

dig_chr.df$plot.group <- paste0(dig_chr.df$sample, "_", dig_chr.df$mz)

#remove any features for which the max intensity is not >= 1e5
dig_chr.df <- dig_chr.df %>% 
    group_by(mz) %>% 
    filter(any(intensity >= 1e5)) 

#manually check chromatograms (peak picking was not stringent, don't know how
#to get around doing this at the moment)
setDF(dig)
dig <- dig[order(dig$rt_min),]
dig_chr.df$mz <- factor(dig_chr.df$mz, levels = unique(round(dig$mz,3))) 

p <- ggplot(data = dig_chr.df %>% filter(source == "L. digitata"), 
            aes(x = rt, y = intensity, group = replicate,
                colour = replicate)) +
    geom_line() +
    facet_grid(rows = vars(as.factor(mz)), scales = "free_y") +
    xlim(0,15) +
    scale_y_continuous(labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 360),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/l_dig_v1.png",
    height = 15, width = 8, units = "in", res = 300)
print(p)
dev.off()

mz_nope <- c(245.022, 311.005, #assume same compound as 300.976
             451.022, #assume same compound as 473.004
             601.639 #not really peaks
)

hyp_chr.df2 <- hyp_chr.df %>% filter(!mz %in% mz_nope) 

hyp_chr.df2$ion2[hyp_chr.df2$mz == 243.018] <- "deoxyhexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 300.976] <- "unknown, m/z = 300.976"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 357.049] <- "unknown, m/z = 357.049"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 214.986] <- "unknown, m/z = 214.986"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 387.06] <- "dehydrated heoxse-deoxyhexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 473.004] <- "unknown, m/z = 473.004"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 229.002] <- "pentose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 225.007] <- "dehydrated deoxyhexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 322.974] <- "unknown, m/z = 322.974 (sulphated)"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 389.076] <- "di-deoxyhexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 259.013] <- "hexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 234.012] <- "di-deoxyhexose disulphate: [M-2H]-2"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 412.944] <- "unknown, m/z = 412.944"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 210.991] <- "dehydrated pentose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 302] <- "unknown, m/z = 302.000"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 198.991] <- "unknown, m/z = 198.991 (sulphated)"


p <- ggplot(data = hyp_chr.df2 %>% filter(source == "L. hyperborea")) +
    geom_line(aes(x = rt, y = intensity, group = plot.group,
                  colour = replicate)) +
    facet_grid(rows = vars(ion2), scales = "free") +
    scale_colour_manual(values = pal_3) +
    xlim(0,7) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/l_hyp_v2.png",
    height = 6, width = 10, units = "in", res = 300)
print(p)
dev.off()

#plot for all samples
hyp_chr.df2$source <- factor(hyp_chr.df2$source, 
                             levels = c("E. bicyclis","L. digitata",
                                        "F. vesiculosus", "L. hyperborea",
                                        "S. cerevisisae", "T. weissflogii",
                                        "solvent blank"))


p <-ggplot(data = hyp_chr.df2) +
    geom_line(aes(x = rt, y = intensity, group = plot.group,
                  colour = replicate)) +
    facet_grid(rows = vars(ion2), cols = vars(source), scales = "free") +
    scale_colour_manual(values = pal_3) +
    xlim(0,10) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          strip.text.x = element_text(vjust = 1, size = 10),
          legend.position = "none",
          strip.background = element_rect(colour = NA))

png(filename = "analysis/analysis_plots/fingerprint_chromatograms_all-samples/20210114_l_hyp.png",
    height = 7, width = 12, units = "in", res = 300)
print(p)
dev.off()


