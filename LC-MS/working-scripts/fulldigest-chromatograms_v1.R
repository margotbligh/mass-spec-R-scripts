#setwd("C:/Users/mbligh/ownCloud/marglyco/Data/LC-MS_data/alpha-mannan/MS31_20220203")
setwd("~/ownCloud/marglyco/Data/LC-MS_data/alpha-mannan/MS31_20220203")
load("RData.RData")

#1: Install packages --------------------------------------------------------
library(xcms)
library(ggplot2)
library(tidyverse)
library(scales)
library(data.table)
library(CAMERA)
library(plyr)


#2: Import and inspect MS data --------------------------------------------------------
#get file paths to mzML files
fp <- dir(path = "mzML", all.files = FALSE,full.names = TRUE)
fp <- fp[grep("full|blank", fp)]

#create phenodata data.frame
#each sample must have a unique name!
pd <- data.frame(name = basename(fp) %>%
                     gsub("MS31_20220203_|_100.*", "", .) %>% 
                     sub("_\\d\\d.mzML", "", .),
                 sample_type = basename(fp) %>% 
                     sub(".*digest.*", "GH99 full digest", .) %>% 
                     sub(".*blank.*", "solvent blank", .),
                 month = basename(fp) %>% 
                     sub(".*april.*", "April", .) %>% 
                     sub(".*june.*", "June", .) %>% 
                     sub(".*solvent.*", "", .),
                 stringsAsFactors = FALSE)

#read in data
all_data <- readMSData(files = fp, pdata = new("NAnnotatedDataFrame", pd), 
                       mode = "onDisk")





#3: Extract chromatograms (targeted)-----
mz.possible <- c(331.0341,
                 412.0605,
                 301.0235,
                 277.0150)

names(mz.possible) <- c("disulphated trisaccharide: [M-2H]-2",
                        "disulphated tetrasaccharide: [M-2H]-2",
                        "trisulphated tetrasaccharide: [M-2H]-3",
                        "pentasulphated hexasaccharide: [M-5H]-5")

#extract chromatograms
chr.list <- list()
for(i in 1:length(mz.possible)){
    chr.list[[i]] <- chromatogram(all_data,
                                  mz = c(mz.possible[i] - 0.001,
                                         mz.possible[i] + 0.001))
}

#build dataframe for plotting with ggplot
chr.df <- data.frame(sample = as.character(),
                     sampletype = as.character(),
                     month = as.character(),
                     mz = as.numeric(), ion = as.character(),
                     rt = as.numeric(), intensity = as.numeric())

for(i in 1:length(mz.possible)){
    chr <- chr.list[[i]]
    mz <- mz.possible[i]
    ion <- names(mz.possible)[i]
    
    lengths <- unlist(lapply(chr, function(x) length(x)))

    sample <- rep(chr$name, times = lengths)
    sampletype <- rep(chr$sample_type, times = lengths)
    month <- rep(chr$month, times = lengths)
    mz <- rep(mz, sum(lengths))
    ion <- rep(ion, sum(lengths))
    rt <- unlist(lapply(chr, function(x) x@rtime/60))
    intensity <- unlist(lapply(chr, function(x) x@intensity))
    
    
    temp <- data.frame(sample = sample, sampletype = sampletype,
                       month = month, mz = mz, ion = ion,
                       rt = rt, intensity = intensity)
    
    chr.df <- rbind(chr.df, temp)
}

chr.df$colourgroup <- paste(chr.df$sampletype, chr.df$month, sep = ", ") %>%
    sub(",\\s$", "", .)

#make palette
pal_3 <- c("#FEC000", "#2BB6AF", "#909193")
names(pal_3) <- chr.df$colourgroup %>% unique()
#show_col(pal_3)

#plot
scientific_function <- function(x){
    text <- gsub("E0", "", gsub("e\\+0", "E", scales::scientific_format()(x)))
    text
}


chr.df$ion <- factor(chr.df$ion,
                     levels = c("disulphated trisaccharide: [M-2H]-2",
                                "disulphated tetrasaccharide: [M-2H]-2",
                                "trisulphated tetrasaccharide: [M-2H]-3",
                                "pentasulphated hexasaccharide: [M-5H]-5"))
chr.df <- chr.df %>% filter(between(rt, 1, 40))

p <- ggplot(data = chr.df, aes(x = rt, y = intensity, group = sample,
                          colour = colourgroup)) +
    geom_line(lwd = 1.2) +
    scale_colour_manual(values = pal_3, name = "") +
    facet_grid(rows = vars(ion), scales = "free_y") +
    xlim(10, 30) +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function) +
    labs(x = "Retention time (min)", y= "Intensity (a.u.)") +
    theme_classic() +
    theme(text = element_text(family = "Arial"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "top",
          strip.background = element_rect(colour = NA))
png(filename = "20220203_fulldigest-chromatograms.png",
    width = 12, height = 6, res = 300, units = "in")    
print(p)
dev.off()





#4: Peak picking----
cwp<-CentWaveParam()
cwp@ppm<-3
cwp@peakwidth<-c(7,65)
cwp@snthresh<-5
cwp@prefilter <- c(3, 1000)

peaks <- findChromPeaks(all_data, param = cwp)

#check chromatograms
chr.list.cwp <- list()
for(i in 1:length(mz.possible)){
    chr.list.cwp[[i]] <- chromatogram(peaks,
                                      mz = c(mz.possible[i] - 0.001,
                                             mz.possible[i] + 0.001))
}

par(mfrow=c(2,2))
for(i in 1:length(mz.possible)){
    plot(chr.list.cwp[[i]], main = names(mz.possible)[i],
         xlim = c(10*60, 25*60))
}

    
    
    
    
    
    
#5: Peak grouping----
pdp <- PeakDensityParam(sampleGroups = peaks$sample_type,
                        binSize = 0.005,
                        bw = 6,
                        minFraction = 0.5) 

par(mar=c(4,4,3,10))
plotChromPeakDensity(chr.list.cwp[[1]], param = pdp, peakPch = 16)
plotChromPeakDensity(chr.list.cwp[[2]], param = pdp, peakPch = 16)
plotChromPeakDensity(chr.list.cwp[[3]], param = pdp, peakPch = 16)
plotChromPeakDensity(chr.list.cwp[[4]], param = pdp, peakPch = 16)

peaks_grouped <- groupChromPeaks(peaks, param = pdp)


#6: Fill in missing peaks----------
fpp <- FillChromPeaksParam()
peaks_filled <- fillChromPeaks(peaks_grouped)


#7: Save diffreport of xdata -----
xset <- as(peaks_filled, "xcmsSet")
sampnames(xset) <- pData(peaks_filled)$name
sampclass(xset) <- pData(peaks_filled)$sample_type

#8. Isotope picking----
an <- xsAnnotate(xset)
an <- groupFWHM(an, perfwhm = 0.6)
an <- findIsotopes(an, mzabs=0.01)
an <- groupCorr(an, cor_eic_th=0.75)
an <- findAdducts(an, polarity="negative")

#9. Peak list filtering and formatting----
#get peak list
pl <-getPeaklist(an)

#filter by blank exclusion (detected peaks)
pl_be <-pl[pl$solvent.blank==0,]

#change NA to 0
pl_be[is.na(pl_be)] <- 0

#add rounded retention time as first colum
pl_be <- cbind(rt_min = round(pl_be$rt/60, 1), pl_be)

#filter by retention time
pl_be_rt <- pl_be %>% filter(between(rt_min, 1, 40))


#10: Collapse features with multiple isotopes -----
peaks.cwp <- peaks
peaks <- pl_be_rt

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
peaks_iso <- peaks_iso[order(isotope_group,isotope_number),]
#get concatenated list of isotopes per group
iso_concat <- peaks_iso[, 
                        list(isotopes = paste(isotopes, 
                                              collapse = ', ')), 
                        by = isotope_group]
#remove duplicates within each isotope group (will keep [M] isotope)
#because of ordering
peaks_iso <- unique(peaks_iso,by = "isotope_group")
#merge to get concatenated isotope lists
peaks_iso <- merge(peaks_iso, iso_concat, by = "isotope_group")
#clean up df
peaks_iso <- peaks_iso %>% 
    select(-c("isotope_group", "isotope_number","isotopes.x"))
names(peaks_iso)[names(peaks_iso) == 'isotopes.y'] <- 'isotopes'

#merge features with and without isotopes
peaks <- rbind.fill(peaks_noiso, peaks_iso)

#11: Annotate features based on predictions----
predicted <- fread("../../predicted_sugars_20220204.txt")

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

#remove annotations with sulphate and M+Cl or M+CHOO
a <- grep("sulphate", peaks_matched$name)
b <- grep("Cl|CHOO", peaks_matched$ion)
c <- intersect(a, b)
peaks_matched <- peaks_matched[-c,]

#order by retention time
peaks_matched <- peaks_matched[order(rt_min),]

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




#12: Extract chromatograms to check annotated peaks-----
peaks_matched.uniq <- peaks_matched %>% distinct(by = id_ion, .keep_all = T)

chr.list.pks <- list()
for(i in 1:nrow(peaks_matched.uniq)){
    chr.list.pks[[i]] <- chromatogram(peaks_filled,
                                      mz = c(peaks_matched.uniq$mz[i] - 0.001,
                                             peaks_matched.uniq$mz[i] + 0.001))
}

#build dataframe for plotting with ggplot
chr.df.pks <- data.frame(sample = as.character(),
                         sampletype = as.character(),
                         month = as.character(),
                         mz = as.numeric(), ion = as.character(),
                         rt = as.numeric(), intensity = as.numeric())

for(i in 1:nrow(peaks_matched.uniq)){
    chr <- chr.list.pks[[i]]
    mz <- peaks_matched.uniq$mz[i]
    ion <- peaks_matched.uniq$id_ion[i]
    
    lengths <- unlist(lapply(chr, function(x) length(x)))
    
    sample <- rep(chr$name, times = lengths)
    sampletype <- rep(chr$sample_type, times = lengths)
    month <- rep(chr$month, times = lengths)
    mz <- rep(mz, sum(lengths))
    ion <- rep(ion, sum(lengths))
    rt <- unlist(lapply(chr, function(x) x@rtime/60))
    intensity <- unlist(lapply(chr, function(x) x@intensity))
    
    
    temp <- data.frame(sample = sample, sampletype = sampletype,
                       month = month, mz = mz, ion = ion,
                       rt = rt, intensity = intensity)
    
    chr.df.pks <- rbind(chr.df.pks, temp)
}

chr.df.pks$colourgroup <- paste(chr.df.pks$sampletype, chr.df.pks$month, sep = ", ") %>%
    sub(",\\s$", "", .)

#remove features that obviously do not have real peaks
chr.df.pks2 <- chr.df.pks %>% 
    filter(!ion %in% c("hex-1-sulphate-2: [M-2H]-2",
                       "hex-5: [M-2H]-2",
                       "hex-7-sulphate-10: [M-9H]-9",
                       "hex-7-sulphate-8: [M-7H]-7"))


#plot
p <- ggplot(data = chr.df.pks2, aes(x = rt, y = intensity, group = sample,
                               colour = colourgroup)) +
    geom_line(lwd = 1.2) +
    scale_colour_manual(values = pal_3, name = "") +
    facet_grid(rows = vars(ion), scales = "free_y") +
    xlim(1, 30) +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function) +
    labs(x = "Retention time (min)", y= "Intensity (a.u.)") +
    theme_classic() +
    theme(text = element_text(family = "Arial"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "top",
          strip.background = element_rect(colour = NA))
png(filename = "20220203_fulldigest-chromatograms_peakpicking.png",
    width = 12, height = 9, res = 300, units = "in")    
print(p)
dev.off()

 x <- cbind(mz = round(unique(chr.df.pks2$mz), 4), name = unique(chr.df.pks2$ion))
fwrite(x, file = "peaks_for_DIA.txt", sep = "\t")

p <- ggplot(data = chr.df.pks2, aes(x = rt, y = intensity, group = sample,
                                    colour = colourgroup)) +
    geom_line(lwd = 1.2) +
    scale_colour_manual(values = pal_3, name = "") +
    facet_wrap(~ion, ncol = 1, scales = "free", strip.position = "right") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function) +
    scale_x_continuous(breaks = seq(1,30, by = 2), limits = c(1,30)) +
    labs(x = "Retention time (min)", y= "Intensity (a.u.)") +
    theme_classic() +
    theme(text = element_text(family = "Arial"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "top",
          strip.background = element_rect(colour = NA))
png(filename = "20220207_fulldigest-chromatograms_peakpicking.png",
    width = 12, height = 9, res = 300, units = "in")    
print(p)
dev.off()
