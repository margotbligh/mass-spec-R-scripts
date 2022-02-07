#analyse full MS run and MS/MS run data together
#running on linux server of MPI

#setwd("/home/mbligh/ownCloud2/marglyco/Data/LC-MS_data/FITDOG/20211129_20211216_analysis")
setwd("~/ownCloud/marglyco/Data/LC-MS_data/FITDOG/20211129_20211216_analysis")
load("analysis/RData/RData_20220126.RData")

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
library(vegan)

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

#12: Collapse features with multiple isotopes -----
setDT(peaks)
#split out features without an isotope detected
peaks_noiso <- peaks[peaks$isotopes=="",]
peaks_iso <- peaks[!peaks$isotopes=="",]
#make column for the isotope group
peaks_iso$isotope_group <- peaks_iso$isotopes %>% 
    sub("\\[M.+", "", .)
#order isotopes within each group correctly
peaks_iso$isotope_number <- peaks_iso$isotopes %>% 
    sub(".*\\[M\\].*", "0", .) %>% 
    sub(".*\\[M\\+", "", .) %>% 
    sub("\\].*", "", .) %>% 
    as.numeric()
setDT(peaks_iso)
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











#13: Annotate features based on predictions ----
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
peaks <- foverlaps(peaks_nopred,predicted)

#change NA values created during matching (features with no match) to be blank
#remove extra columns
peaks$mz <- NULL
peaks$mzmin <- NULL
peaks$mzmax <- NULL

peaks <- peaks %>% 
    replace_na(list("name"="", "ion"= "", "dp" = ""))

#remove annotations with sulphate or carboxyl and M+Cl or M+CHOO
a <- grep("sulphate|carboxyl", peaks$name)
b <- grep("Cl|CHOO", peaks$ion)
c <- intersect(a, b)
peaks <- peaks[-c,]

#order by retention time
peaks <- peaks[order(rt_min),]

#make id and ion column
peaks$id_ion_individual <- paste0(peaks$name, ": ", peaks$ion)
peaks$id_ion_individual[peaks$name == ""] <- paste0("unknown, m/z=",
                                                    round(peaks$i.mz[peaks$name == ""], 3))

#aggregate so that if there are multiple predictions for one feature
#they are shown in the same row. delete all of the other extra columns added 
#during matching
names <- setdiff(names(peaks), names(predicted))
names <- names[!names == "id_ion_individual"]
setDF(peaks)
peaks_uncollapsed <- peaks
peaks <- peaks_uncollapsed %>% 
    dplyr::group_by(across(all_of(names))) %>% 
    dplyr::mutate(id_ion = paste0(id_ion_individual, collapse = ", ")) %>% 
    ungroup() %>% 
    distinct(across(all_of(c(names, "id_ion"))))
names(peaks)[names(peaks) == 'i.mz'] <- 'mz'
names(peaks)[names(peaks) == 'i.mzmin'] <- 'mzmin'
names(peaks)[names(peaks) == 'i.mzmax'] <- 'mzmax'


#reorder columns
peaks_old <- peaks

colOrder <- names(peaks)
colOrder <- colOrder[1:length(colOrder)-1]
colOrder <- c("id_ion", colOrder)

setDF(peaks)
peaks <- peaks_old[,colOrder]



#15: Look at annotated fingerprints of glycans, EIC for all features in all samples------
par(mfrow=c(1,1))

pal_6 <- rep(c("#FEC000", "#2BB6AF", "#D1D1FE"), 2)
pal_3 <- c("#FEC000", "#2BB6AF", "#D1D1FE")

    #F. vesiculosus fucoidan-----
ves <- peaks %>% 
    filter(F..vesiculosus >= 3 & if_all(contains("_ves"), ~ . >= 1e6) |
               F..vesiculosus >= 4) 

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
ves_chr.dfb <- ves_chr.df %>% 
    group_by(mz) %>% 
    filter(any(intensity >= 1e5)) 

#manually check chromatograms (peak picking was not stringent, don't know how
#to get around doing this at the moment)
setDF(ves)
ves <- ves[order(ves$rt_min),]
ves_chr.dfb$mz <- factor(ves_chr.dfb$mz, levels = unique(round(ves$mz,3))) 

ves_chr.dfb$runtype <- ves_chr.dfb$sample %>% 
    sub("20211129_.*", "full MS", .) %>% 
    sub("20211216_.*", "MS/MS", .)
ves_chr.dfb$runtype <- factor(ves_chr.dfb$runtype, levels = c("full MS",
                                                              "MS/MS"))
ves_chr.dfb$replicate <- ves_chr.dfb$sample %>% 
    sub(".*_\\D+", "", .)

p <- ggplot(data = ves_chr.dfb %>% filter(source == "F. vesiculosus"), 
            aes(x = rt, y = intensity, group = sample,
                colour = replicate, linetype = runtype)) +
    geom_line() +
    scale_colour_manual(values = pal_3) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    facet_grid(rows = vars(as.factor(mz)), scales = "free_y") +
    xlim(0,10) +
    scale_y_continuous(labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 360),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/f_ves_v1.png",
    height = 25, width = 10, units = "in", res = 300)
print(p)
dev.off()

mz_nope <- c(
    286.895, 525.826, 584.793, 236.953, 314.895, 616.783, 452.805,#not peaks
    406.074, #assume same as 405.071
    344.956, #assume same compound as 322.974
    391.071, 313.908, #assume same compound as 447.035
    473.004, #assume same compound as 451.022,
    214.986, #assume same compound as 229.002
    300.976, 245.014, 302.973,#same as 243.018
    282.023, #assume same compound as 565.055
    222.09 #assume same as 360.079
)


ves_chr.df2 <- ves_chr.dfb %>% filter(!mz %in% mz_nope)

ves_chr.df2$ion2[ves_chr.df2$mz == 198.991] <- "unknown, m/z = 198.991 (sulphated) *"
ves_chr.df2$ion2[ves_chr.df2$mz == 221.066] <- "unknown, m/z = 221.066"
ves_chr.df2$ion2[ves_chr.df2$mz == 301.023] <- "tetra-hexose trisulphate: [M-3H]-3"
ves_chr.df2$ion2[ves_chr.df2$mz == 331.08] <- "deoxyhexose-pentose: [M+Cl]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 210.991] <- "dehydrated pentose monosulphate: [M-H]-*"
ves_chr.df2$ion2[ves_chr.df2$mz == 234.013] <- "di-deoxyhexose disulphate: [M-2H]-2 *"
ves_chr.df2$ion2[ves_chr.df2$mz == 199.022] <- "unknown, m/z = 199.022"
ves_chr.df2$ion2[ves_chr.df2$mz == 215.033] <- "hexose: [M+Cl]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 345.095] <- "di-deoxyhexose: [M+Cl]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 259.013] <- "hexose monosulphate: [M-H]- *"
ves_chr.df2$ion2[ves_chr.df2$mz == 225.007] <- "dehydrated deoxyhexose monosulphate: [M-H]- *"
ves_chr.df2$ion2[ves_chr.df2$mz == 229.002] <- "pentose monosulphate: [M-H]- *"
ves_chr.df2$ion2[ves_chr.df2$mz == 405.071] <- "heoxse-deoxyhexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 208.065] <- "unknown, m/z = 208.065"
ves_chr.df2$ion2[ves_chr.df2$mz == 322.974] <- "unknown, m/z = 322.974 (sulphated) *"
ves_chr.df2$ion2[ves_chr.df2$mz == 535.135] <- "tri-deoxyhexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 457.064] <- "unknown, m/z = 457.064"
ves_chr.df2$ion2[ves_chr.df2$mz == 219.007] <- "unknown, m/z = 219.007"
ves_chr.df2$ion2[ves_chr.df2$mz == 447.035] <- "unknown, m/z = 447.035"
ves_chr.df2$ion2[ves_chr.df2$mz == 451.022] <- "unknown, m/z = 451.022"
ves_chr.df2$ion2[ves_chr.df2$mz == 243.018] <- "deoxyhexose monosulphate: [M-H]- *"
ves_chr.df2$ion2[ves_chr.df2$mz == 241.002] <- "dehydrated hexose monosulphate: [M-H]- *"
ves_chr.df2$ion2[ves_chr.df2$mz == 358.935] <- "unknown, m/z = 358.935"
ves_chr.df2$ion2[ves_chr.df2$mz == 371.065] <- "dehydrated di-deoxyhexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 565.055] <- "unknown, m/z = 565.055"
ves_chr.df2$ion2[ves_chr.df2$mz == 389.076] <- "di-deoxyhexose monosulphate: [M-H]- *"
ves_chr.df2$ion2[ves_chr.df2$mz == 385.045] <- "unknown, m/z = 385.045"
ves_chr.df2$ion2[ves_chr.df2$mz == 311.005] <- "unknown, m/z = 311.005"
ves_chr.df2$ion2[ves_chr.df2$mz == 387.06] <- "dehydrated heoxse-deoxyhexose monosulphate: [M-H]-"
ves_chr.df2$ion2[ves_chr.df2$mz == 360.079] <- "unknown, m/z = 360.079"


p <- ggplot(data = ves_chr.df2 %>% filter(source == "F. vesiculosus")) +
    geom_line(aes(x = rt, y = intensity, group = sample,
                  colour = replicate, linetype = runtype)) +
    facet_grid(rows = vars(ion2), scales = "free") +
    scale_colour_manual(values = pal_3, name = "Replicate") +
    scale_linetype_manual(values = c("solid", "dashed"), name = "Run type") +
    xlim(0,10) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "top",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/analysis_plots/fingerprint_chromatograms_all-samples/20210201_f_ves.png",
    height = 12, width = 12, units = "in", res = 300)
print(p)
dev.off()


rm(ves_chr)

    #L. hyperborea fucoidan-----
hyp <- peaks %>% 
    filter(L..hyperborea >= 3 & if_all(contains("_hyp"), ~ . >= 1e6) |
               L..hyperborea >= 4)
hyp.mz <- round(hyp$mz, 3)
names(hyp.mz) <- hyp$id_ion

hyp.mz.unique <- unique(hyp.mz)
names(hyp.mz.unique) <- names(hyp.mz)[match(hyp.mz.unique, hyp.mz)]

#extract chromatograms
#takes wayyy to long actually to extract for all files, from
#now on only extract for specific samples
hyp_chr <- list()
res.hyp <- filterFile(res, which(pd$source == "L. hyperborea"))
for(i in 1:length(hyp.mz.unique)){
    hyp_chr[[i]] <- chromatogram(res.hyp, mz = c(hyp.mz.unique[i]-0.001, 
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
hyp_chr.dfb <- hyp_chr.df %>% 
    group_by(mz) %>% 
    filter(any(intensity >= 1e5)) 

#manually check chromatograms (peak picking was not stringent, don't know how
#to get around doing this at the moment)
setDF(hyp)
hyp <- hyp[order(hyp$rt_min),]
hyp_chr.dfb$mz <- factor(hyp_chr.dfb$mz, levels = unique(round(hyp$mz,3))) 

hyp_chr.dfb$runtype <- hyp_chr.dfb$sample %>% 
    sub("20211129_.*", "full MS", .) %>% 
    sub("20211216_.*", "MS/MS", .)
hyp_chr.dfb$runtype <- factor(hyp_chr.dfb$runtype, levels = c("full MS",
                                                              "MS/MS"))
hyp_chr.dfb$replicate <- hyp_chr.dfb$sample %>% 
    sub(".*_\\D+", "", .)

ggplot(data = hyp_chr.dfb, 
            aes(x = rt, y = intensity, group = sample,
                colour = replicate, linetype = runtype)) +
    geom_line() +
    scale_colour_manual(values = pal_3) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    facet_grid(rows = vars(as.factor(mz)), scales = "free_y") +
    xlim(0,20) +
    scale_y_continuous(labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 360),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/l_hyp_v1.png",
    height = 25, width = 10, units = "in", res = 300)
print(p)
dev.off()


mz_nope <- c(
    433.831, 334.938, 298.937, 236.953, 629.747, 525.826,522.842, 493.81, 461.838, #not real peaks
    313.908, 257.925, #not real peaks
    324.97, 412.944, 346.952, 404.912, 462.871, #assume same compound as 402.915
    234.514, 235.01, 469.033, #assume same compound as 234.012 (isotopes + M-H)
    344.956, #assume same compound as 322.974
    473.004, #assume same compound as 451.022,
    292.964, #assume same compound as 314.946
    222.09, 246.069, #assume same as 360.079
    300.976, 245.014, 302.973,#same as 243.018
    358.935, #assume same as 230.998
    286.961, #assume same as 214.986
    253.987, 530.964 #assume same as 226.005
)

hyp_chr.df2 <- hyp_chr.dfb %>% filter(!mz %in% mz_nope) 


hyp_chr.df2$ion2[hyp_chr.df2$mz == 284.982] <- "unknown, m/z = 284.982"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 198.991] <- "unknown, m/z = 198.991 (sulphated) *"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 535.135] <- "tri-deoxyhexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 210.991] <- "dehydrated pentose monosulphate: [M-H]-*"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 215.033] <- "hexose: [M+Cl]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 402.915] <- "unknown, m/z = 402.915"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 234.012] <- "di-deoxyhexose disulphate: [M-2H]-2 *"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 259.013] <- "hexose monosulphate: [M-H]- *"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 322.974] <- "unknown, m/z = 322.974 (sulphated) *"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 225.007] <- "dehydrated deoxyhexose monosulphate: [M-H]- *"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 243.018] <- "deoxyhexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 389.076] <- "di-deoxyhexose monosulphate: [M-H]- *"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 451.022] <- "unknown, m/z = 451.022"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 314.946] <- "unknown, m/z = 314.946"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 360.079] <- "unknown, m/z = 360.079"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 387.06] <- "dehydrated heoxse-deoxyhexose monosulphate: [M-H]-"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 194.96] <- "unknown, m/z = 194.96"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 311.005] <- "unknown, m/z = 311.005"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 357.049] <- "unknown, m/z = 357.049"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 342.941] <- "unknown, m/z = 342.941"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 180.981] <- "unknown, m/z = 180.981"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 233.005] <- "unknown, m/z = 233.005"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 229.002] <- "pentose monosulphate: [M-H]- *"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 230.998] <- "unknown, m/z = 230.998"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 214.986] <- "unknown, m/z = 214.986"
hyp_chr.df2$ion2[hyp_chr.df2$mz == 226.005] <- "unknown, m/z = 226.005"


#plot again
p <- ggplot(data = hyp_chr.df2) +
    geom_line(aes(x = rt, y = intensity, group = plot.group,
                  colour = replicate, linetype = runtype)) +
    facet_grid(rows = vars(ion2), scales = "free") +
    scale_colour_manual(values = pal_3, name = "Replicate") +
    scale_linetype_manual(values = c("solid", "dotted"), name = "Run type") +
    xlim(0,15) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "top",
          strip.background = element_rect(colour = NA))

png(filename = "analysis/analysis_plots/fingerprint_chromatograms_all-samples/20220202_l_hyp.png",
    height = 7, width = 12, units = "in", res = 300)
print(p)
dev.off()

rm(hyp_chr)


    #L. digitata laminarin-----
dig <- peaks %>% 
    filter(L..digitata >= 3 & if_all(contains("_dig"), ~ . >= 1e6) |
               L..digitata >= 4) 

dig.mz <- round(dig$mz, 3)
names(dig.mz) <- dig$id_ion

dig.mz.unique <- unique(dig.mz)
names(dig.mz.unique) <- names(dig.mz)[match(dig.mz.unique, dig.mz)]

#extract chromatograms
res.dig <- filterFile(res, which(pd$source == "L. digitata"))

dig_chr <- list()
for(i in 1:length(dig.mz.unique)){
    dig_chr[[i]] <- chromatogram(res.dig, mz = c(dig.mz.unique[i]-0.001, 
                                                 dig.mz.unique[i]+0.001))
}

save(dig_chr, file = "analysis/RData/dig_chr.RData")


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
dig_chr.dfb <- dig_chr.df %>% 
    group_by(mz) %>% 
    filter(any(intensity >= 1e5)) 

#manually check chromatograms (peak picking was not stringent, don't know how
#to get around doing this at the moment)
setDF(dig)
dig <- dig[order(dig$rt_min),]
dig_chr.dfb$mz <- factor(dig_chr.dfb$mz, levels = unique(round(dig$mz,3))) 

dig_chr.dfb$runtype <- dig_chr.dfb$sample %>% 
    sub("20211129_.*", "full MS", .) %>% 
    sub("20211216_.*", "MS/MS", .)
dig_chr.dfb$runtype <- factor(dig_chr.dfb$runtype, levels = c("full MS",
                                                              "MS/MS"))
dig_chr.dfb$replicate <- dig_chr.dfb$sample %>% 
    sub(".*_\\D+", "", .)

p <- ggplot(data = dig_chr.dfb, 
            aes(x = rt, y = intensity, group = sample,
                colour = replicate, linetype = runtype)) +
    geom_line() +
    scale_colour_manual(values = pal_3) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    facet_grid(rows = vars(as.factor(mz)), scales = "free_y") +
    xlim(0,15) +
    scale_y_continuous(labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 360),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/l_dig_v1.png",
    height = 25, width = 10, units = "in", res = 300)
print(p)
dev.off()


mz_nope <- c(
    382.861, 553.778, 346.917, 421.845, 422.842, 360.879, 493.81, 460.839, #not real peaks
    196.903, 522.842, 525.826, 551.77, 286.95, #not real peaks
    584.793, 236.953, 286.895, #not real peaks
    207.004, 217.033, #assume same compound as 323.01
    335.902, #assume same compound as 257.925
    222.09, #assume same as 360.079
    194.994 #assume same as 205.023
) 

dig_chr.df2 <- dig_chr.dfb %>% filter(!mz %in% mz_nope)

dig_chr.df2$ion2[dig_chr.df2$mz == 268.169] <- "unknown, m/z = 268.169"
dig_chr.df2$ion2[dig_chr.df2$mz == 323.01] <- "unknown, m/z = 323.01"
dig_chr.df2$ion2[dig_chr.df2$mz == 377.085] <- "di-hexose: [M+Cl]-"
dig_chr.df2$ion2[dig_chr.df2$mz == 539.139] <- "tri-hexose: [M+Cl]-"
dig_chr.df2$ion2[dig_chr.df2$mz == 360.079] <- "unknown, m/z = 360.079"
dig_chr.df2$ion2[dig_chr.df2$mz == 257.925] <- "unknown, m/z = 257.925"
dig_chr.df2$ion2[dig_chr.df2$mz == 205.023] <- "unknown, m/z = 205.023"
dig_chr.df2$ion2[dig_chr.df2$mz == 208.065] <- "unknown, m/z = 208.065"
dig_chr.df2$ion2[dig_chr.df2$mz == 375.07] <- "unknown, m/z = 375.070"
dig_chr.df2$ion2[dig_chr.df2$mz == 215.032] <- "hexose: [M+Cl]-"
dig_chr.df2$ion2[dig_chr.df2$mz == 187.99] <- "unknown, m/z = 187.99"
dig_chr.df2$ion2[dig_chr.df2$mz == 347.075] <- "hexose-pentose: [M+Cl]-"
dig_chr.df2$ion2[dig_chr.df2$mz == 237.061] <- "unknown, m/z = 237.061"


#plot again
p <- ggplot(data = dig_chr.df2) +
    geom_line(aes(x = rt, y = intensity, group = plot.group,
                  colour = replicate, linetype = runtype)) +
    facet_grid(rows = vars(ion2), scales = "free") +
    scale_colour_manual(values = pal_3, name = "Replicate") +
    scale_linetype_manual(values = c("solid", "dotted"), name = "Run type") +
    xlim(0,15) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "top",
          strip.background = element_rect(colour = NA))

png(filename = "analysis/analysis_plots/fingerprint_chromatograms_all-samples/20220202_l_dig.png",
    height = 7, width = 12, units = "in", res = 300)
print(p)
dev.off()

rm(dig_chr)


    #E. bicyclis laminarin-----
bic <- peaks %>% 
    filter(E..bicyclis >= 3 & if_all(contains("_bic"), ~ . >= 1e6) |
               E..bicyclis >= 4) 

bic.mz <- round(bic$mz, 3)
names(bic.mz) <- bic$id_ion

bic.mz.unique <- unique(bic.mz)
names(bic.mz.unique) <- names(bic.mz)[match(bic.mz.unique, bic.mz)]


#extract chromatograms
bic_chr <- list()
res.bic <- filterFile(res, which(pd$source == "E. bicyclis"))

for(i in 1:length(bic.mz.unique)){
    bic_chr[[i]] <- chromatogram(res.bic, mz = c(bic.mz.unique[i]-0.001, 
                                                 bic.mz.unique[i]+0.001))
}

save(bic_chr, file = "analysis/RData/bic_chr.RData")


bic_chr.df <- data.frame(sample = as.character(),
                         substrate = as.character(),
                         source = as.character(),
                         replicate = as.character(),
                         mz = as.numeric(), rt = as.numeric(), 
                         intensity = as.numeric(),
                         ion = as.character())
for(i in 1:length(bic.mz.unique)){
    chr <- bic_chr[[i]]
    mz <- bic.mz.unique[i]
    annot <- names(bic.mz.unique)[i]
    
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
    bic_chr.df <- rbind(bic_chr.df, tmp)
}

bic_chr.df$plot.group <- paste0(bic_chr.df$sample, "_", bic_chr.df$mz)

#remove any features for which the max intensity is not >= 1e5
bic_chr.dfb <- bic_chr.df %>% 
    group_by(mz) %>% 
    filter(any(intensity >= 1e5)) 

#chromatograms, check manually
setDF(bic)
bic <- bic[order(bic$rt_min),]
bic_chr.dfb$mz <- factor(bic_chr.dfb$mz, levels = unique(round(bic$mz,3))) 

bic_chr.dfb$runtype <- bic_chr.dfb$sample %>% 
    sub("20211129_.*", "full MS", .) %>% 
    sub("20211216_.*", "MS/MS", .)
bic_chr.dfb$runtype <- factor(bic_chr.dfb$runtype, levels = c("full MS",
                                                              "MS/MS"))
bic_chr.dfb$replicate <- bic_chr.dfb$sample %>% 
    sub(".*_\\D+", "", .)

p <- ggplot(data = bic_chr.dfb, 
            aes(x = rt, y = intensity, group = sample,
                colour = replicate, linetype = runtype)) +
    geom_line() +
    scale_colour_manual(values = pal_3) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    facet_grid(rows = vars(as.factor(mz)), scales = "free_y") +
    xlim(0,12) +
    scale_y_continuous(labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 360),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/e_bic_v1.png",
    height = 25, width = 10, units = "in", res = 300)
print(p)
dev.off()



mz_nope <- c(
    493.81, 551.77, 322.893, 310.886, 334.937, 369.899, 584.793, 620.766, #not real peaks
    395.867, #assume same compound as 257.925
    300.976, 245.014, #same as 243.018
    222.09 #assume same as 360.079
    )

bic_chr.df2 <- bic_chr.dfb %>% filter(!mz %in% mz_nope)

#bic_chr.dfb$ion[bic_chr.dfb$mz == 222.09] %>% unique()

bic_chr.df2$ion2[bic_chr.df2$mz == 377.085] <- "di-hexose: [M+Cl]-"
bic_chr.df2$ion2[bic_chr.df2$mz == 237.061] <- "unknown, m/z = 237.061"
bic_chr.df2$ion2[bic_chr.df2$mz == 215.032] <- "hexose: [M+Cl]-"
bic_chr.df2$ion2[bic_chr.df2$mz == 257.925] <- "unknown, m/z = 257.925"
bic_chr.df2$ion2[bic_chr.df2$mz == 229.002] <- "pentose monosulphate: [M-H]- *"
bic_chr.df2$ion2[bic_chr.df2$mz == 243.018] <- "deoxyhexose monosulphate: [M-H]- *"
bic_chr.df2$ion2[bic_chr.df2$mz == 360.079] <- "unknown, m/z = 360.079"

#plot again
p <- ggplot(data = bic_chr.df2) +
    geom_line(aes(x = rt, y = intensity, group = plot.group,
                  colour = replicate, linetype = runtype)) +
    facet_grid(rows = vars(ion2), scales = "free") +
    scale_colour_manual(values = pal_3, name = "Replicate") +
    scale_linetype_manual(values = c("solid", "dotted"), name = "Run type") +
    xlim(0,15) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "top",
          strip.background = element_rect(colour = NA))

png(filename = "analysis/analysis_plots/fingerprint_chromatograms_all-samples/20220202_e_bic.png",
    height = 7, width = 12, units = "in", res = 300)
print(p)
dev.off()

rm(bic_chr)


    #S. cerevisiae mannan-----
yea <- peaks %>% 
    filter(S..cerevisisae >= 3 & if_all(contains("_yea"), ~ . >= 1e6) |
               S..cerevisisae >= 4) 

yea.mz <- round(yea$mz, 3)
names(yea.mz) <- yea$id_ion

yea.mz.unique <- unique(yea.mz)
names(yea.mz.unique) <- names(yea.mz)[match(yea.mz.unique, yea.mz)]


#extract chromatograms
yea_chr <- list()
res.yea <- filterFile(res, which(pd$source == "S. cerevisisae"))

for(i in 1:length(yea.mz.unique)){
    yea_chr[[i]] <- chromatogram(res.yea, mz = c(yea.mz.unique[i]-0.001, 
                                                 yea.mz.unique[i]+0.001))
}

save(yea_chr, file = "analysis/RData/yea_chr.RData")


yea_chr.df <- data.frame(sample = as.character(),
                         substrate = as.character(),
                         source = as.character(),
                         replicate = as.character(),
                         mz = as.numeric(), rt = as.numeric(), 
                         intensity = as.numeric(),
                         ion = as.character())
for(i in 1:length(yea.mz.unique)){
    chr <- yea_chr[[i]]
    mz <- yea.mz.unique[i]
    annot <- names(yea.mz.unique)[i]
    
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
    yea_chr.df <- rbind(yea_chr.df, tmp)
}

yea_chr.df$plot.group <- paste0(yea_chr.df$sample, "_", yea_chr.df$mz)

#remove any features for which the max intensity is not >= 1e5
yea_chr.dfb <- yea_chr.df %>% 
    group_by(mz) %>% 
    filter(any(intensity >= 1e5)) 

#chromatograms, check manually
setDF(yea)
yea <- yea[order(yea$rt_min),]
yea_chr.dfb$mz <- factor(yea_chr.dfb$mz, levels = unique(round(yea$mz,3))) 

yea_chr.dfb$runtype <- yea_chr.dfb$sample %>% 
    sub("20211129_.*", "full MS", .) %>% 
    sub("20211216_.*", "MS/MS", .)
yea_chr.dfb$runtype <- factor(yea_chr.dfb$runtype, levels = c("full MS",
                                                              "MS/MS"))
yea_chr.dfb$replicate <- yea_chr.dfb$sample %>% 
    sub(".*_\\D+", "", .)

p <- ggplot(data = yea_chr.dfb, 
            aes(x = rt, y = intensity, group = sample,
                colour = replicate, linetype = runtype)) +
    geom_line() +
    scale_colour_manual(values = pal_3) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    facet_grid(rows = vars(as.factor(mz)), scales = "free_y") +
    xlim(0,15) +
    scale_y_continuous(labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 360),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/s_cer_v1.png",
    height = 25, width = 10, units = "in", res = 300)
print(p)
dev.off()

mz_nope <- c(
    314.895, 334.938, 312.905,584.793,353.857,525.826,522.842,412.824,
    196.903,420.846,364.889, 394.847, 236.953, 443.866,
    411.027, #same as 237.061
    459.045,387.114, #same as 377.085
    349.072, #assume isotope of 347.074
    222.09, 349.106, 246.069 #assume same as 360.079 
)

yea_chr.df2 <- yea_chr.dfb %>% filter(!mz %in% mz_nope)

yea_chr.df2$ion2[yea_chr.df2$mz == 399.114] <- "unknown, m/z = 399.114"
yea_chr.df2$ion2[yea_chr.df2$mz == 237.061] <- "unknown, m/z = 237.061"
yea_chr.df2$ion2[yea_chr.df2$mz == 377.085] <- "di-hexose: [M+Cl]-"
yea_chr.df2$ion2[yea_chr.df2$mz == 347.075] <- "hexose-pentose: [M+Cl]-"
yea_chr.df2$ion2[yea_chr.df2$mz == 215.032] <- "hexose: [M+Cl]-"
yea_chr.df2$ion2[yea_chr.df2$mz == 256.059] <- "unknown, m/z=256.059"
yea_chr.df2$ion2[yea_chr.df2$mz == 284.054] <- "unknown, m/z=284.054"
yea_chr.df2$ion2[yea_chr.df2$mz == 360.079] <- "unknown, m/z=360.079"
yea_chr.df2$ion2[yea_chr.df2$mz == 365.164] <- "unknown, m/z=365.164"

#plot again
p <- ggplot(data = yea_chr.df2) +
    geom_line(aes(x = rt, y = intensity, group = plot.group,
                  colour = replicate, linetype = runtype)) +
    facet_grid(rows = vars(ion2), scales = "free") +
    scale_colour_manual(values = pal_3, name = "Replicate") +
    scale_linetype_manual(values = c("solid", "dotted"), name = "Run type") +
    xlim(0,15) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "top",
          strip.background = element_rect(colour = NA))

png(filename = "analysis/analysis_plots/fingerprint_chromatograms_all-samples/20220207_s_cer.png",
    height = 4.5, width = 12, units = "in", res = 300)
print(p)
dev.off()

rm(yea_chr)


    #T. weisflogii mannan-----
twa<- peaks %>% 
    filter(T..weissflogii >= 3 & if_all(contains("_twa"), ~ . >= 1e6) |
               T..weissflogii >= 4) 

twa.mz <- round(twa$mz, 3)
names(twa.mz) <- twa$id_ion

twa.mz.unique <- unique(twa.mz)
names(twa.mz.unique) <- names(twa.mz)[match(twa.mz.unique, twa.mz)]


#extract chromatograms
twa_chr <- list()
res.twa <- filterFile(res, which(pd$source == "T. weissflogii"))
for(i in 1:length(twa.mz.unique)){
    twa_chr[[i]] <- chromatogram(res.twa, mz = c(twa.mz.unique[i]-0.001, 
                                                 twa.mz.unique[i]+0.001))
}

save(twa_chr, file = "analysis/RData/twa_chr.RData")


twa_chr.df <- data.frame(sample = as.character(),
                         substrate = as.character(),
                         source = as.character(),
                         replicate = as.character(),
                         mz = as.numeric(), rt = as.numeric(), 
                         intensity = as.numeric(),
                         ion = as.character())
for(i in 1:length(twa.mz.unique)){
    chr <- twa_chr[[i]]
    mz <- twa.mz.unique[i]
    annot <- names(twa.mz.unique)[i]
    
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
    twa_chr.df <- rbind(twa_chr.df, tmp)
}

twa_chr.df$plot.group <- paste0(twa_chr.df$sample, "_", twa_chr.df$mz)

#remove any features for which the max intensity is not >= 1e5
twa_chr.dfb <- twa_chr.df %>% 
    group_by(mz) %>% 
    filter(any(intensity >= 1e5)) 

#chromatograms, check manually
setDF(twa)
twa <- twa[order(twa$rt_min),]
twa_chr.dfb$mz <- factor(twa_chr.dfb$mz, levels = unique(round(twa$mz,3))) 

twa_chr.dfb$runtype <- twa_chr.dfb$sample %>% 
    sub("20211129_.*", "full MS", .) %>% 
    sub("20211216_.*", "MS/MS", .)
twa_chr.dfb$runtype <- factor(twa_chr.dfb$runtype, levels = c("full MS",
                                                              "MS/MS"))
twa_chr.dfb$replicate <- twa_chr.dfb$sample %>% 
    sub(".*_\\D+", "", .)

p <- ggplot(data = twa_chr.dfb, 
            aes(x = rt, y = intensity, group = sample,
                colour = replicate, linetype = runtype)) +
    geom_line() +
    scale_colour_manual(values = pal_3) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    facet_grid(rows = vars(as.factor(mz)), scales = "free_y") +
    xlim(0,15) +
    scale_y_continuous(labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 360),
          legend.position = "none",
          strip.background = element_rect(colour = NA))
png(filename = "analysis/processing_plots/eic_fingerprints/t_wei_v1.png",
    height = 25, width = 10, units = "in", res = 300)
print(p)
dev.off()

mz_nope <- c(
    334.937, 426.88, 334.863, 196.903, 236.953,372.831, 382.861, 199.954,
    313.908, 393.87, 339.926,
    379.082, 387.114, #same as 377.085
    368.027, #assume same as 405.071
    474.984, #assume same as 225.997
    261.009, #assume as same as 259.013
    191.019, #assume same as 327.987
    286.961, #assume same as 230.998
    521.119, 338.016, 241.998, #assume same 375.06
    472.968, 450.986, 342.049, #assume same as 224.989
    300.976, 245.014, #assume same as 311.005
    210.007, #assume same as 225.491
    222.09, 286.039 #assume same as 360.079
    
)

twa_chr.df2 <- twa_chr.dfb %>% filter(!mz %in% mz_nope)

twa_chr.df2$ion2[twa_chr.df2$mz == 250.008] <- "di-hexose disulphate: [M-2H]-2"
twa_chr.df2$ion2[twa_chr.df2$mz == 285.971] <- "unknown, m/z = 285.971"
twa_chr.df2$ion2[twa_chr.df2$mz == 184.976] <- "unknown, m/z = 184.976"
twa_chr.df2$ion2[twa_chr.df2$mz == 377.085] <- "di-hexose: [M+Cl]-"
twa_chr.df2$ion2[twa_chr.df2$mz == 421.066] <- "di-hexose monosulphate: [M-H]-"
twa_chr.df2$ion2[twa_chr.df2$mz == 287.007] <- "unknown, m/z = 287.007"
twa_chr.df2$ion2[twa_chr.df2$mz == 391.055] <- "hexose-pentose monosulphate: [M-H]-"
twa_chr.df2$ion2[twa_chr.df2$mz == 215.032] <- "hexose: [M+Cl]-"
twa_chr.df2$ion2[twa_chr.df2$mz == 405.071] <- "hexose-deoxyhexose monosulphate: [M-H]-"
twa_chr.df2$ion2[twa_chr.df2$mz == 551.13] <- "hexose-di-deoxyhexose monosulphate: [M-H]-"
twa_chr.df2$ion2[twa_chr.df2$mz == 256.059] <- "unknown, m/z = 256.059"
twa_chr.df2$ion2[twa_chr.df2$mz == 225.997] <- "unknown, m/z = 225.997"
twa_chr.df2$ion2[twa_chr.df2$mz == 208.065] <- "unknown, m/z = 208.065"
twa_chr.df2$ion2[twa_chr.df2$mz == 259.013] <- "hexose monosulphate: [M-H]- *"
twa_chr.df2$ion2[twa_chr.df2$mz == 357.998] <- "unknown, m/z = 357.998"
twa_chr.df2$ion2[twa_chr.df2$mz == 316.971] <- "unknown, m/z = 316.971"
twa_chr.df2$ion2[twa_chr.df2$mz == 371.029] <- "unknown, m/z = 371.029"
twa_chr.df2$ion2[twa_chr.df2$mz == 408.102] <- "unknown, m/z = 408.102"
twa_chr.df2$ion2[twa_chr.df2$mz == 229.002] <- "pentose monosulphate: [M-H]- *"
twa_chr.df2$ion2[twa_chr.df2$mz == 241.002] <- "dehydrated hexose monosulphate: [M-H]- *"
twa_chr.df2$ion2[twa_chr.df2$mz == 327.987] <- "unknown, m/z = 327.987"
twa_chr.df2$ion2[twa_chr.df2$mz == 230.998] <- "unknown, m/z = 230.998"
twa_chr.df2$ion2[twa_chr.df2$mz == 375.06] <- "unknown, m/z = 375.07"
twa_chr.df2$ion2[twa_chr.df2$mz == 224.989] <- "unknown, m/z = 224.989"
twa_chr.df2$ion2[twa_chr.df2$mz == 311.005] <- "unknown, m/z = 311.005"
twa_chr.df2$ion2[twa_chr.df2$mz == 282.029] <- "unknown, m/z = 282.029"
twa_chr.df2$ion2[twa_chr.df2$mz == 225.491] <- "unknown, m/z = 225.491"
twa_chr.df2$ion2[twa_chr.df2$mz == 360.079] <- "unknown, m/z = 360.079"

#plot again
p <- ggplot(data = twa_chr.df2) +
    geom_line(aes(x = rt, y = intensity, group = plot.group,
                  colour = replicate, linetype = runtype)) +
    facet_grid(rows = vars(ion2), scales = "free") +
    scale_colour_manual(values = pal_3, name = "Replicate") +
    scale_linetype_manual(values = c("solid", "dotted"), name = "Run type") +
    xlim(0,12) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function)+
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          legend.position = "top",
          strip.background = element_rect(colour = NA))

png(filename = "analysis/analysis_plots/fingerprint_chromatograms_all-samples/20220207_t_wei.png",
    height = 7, width = 12, units = "in", res = 300)
print(p)
dev.off()

rm(twa_chr)


#17: Bubble plot of all fingerprints-----
all_feat <- rbind(ves, hyp, dig, bic, yea, twa)
all_feat <- all_feat %>% 
    group_by(mz, rt) %>% 
    distinct(.keep_all = T)

#only keep features that are also in final dataframes
all_feat$mz_round <- round(all_feat$mz, 3)
all_feat2 <- all_feat[all_feat$mz_round %in% c(ves_chr.df2$mz, hyp_chr.df2$mz, 
                                               dig_chr.df2$mz, bic_chr.df2$mz,
                                               yea_chr.df2$mz, twa_chr.df2$mz),] 

##add annotations
#format df
all_feat3 <- all_feat2
all_feat3$mz <-NULL
names(all_feat3)[names(all_feat3) == "mz_round"] <- "mz"

#vesiculosus
ves_chr.df3 <- ves_chr.df2 %>% distinct(ion2, .keep_all = TRUE)
ves_chr.df3$rt <- NULL
ves_chr.df3$mz <- as.numeric(paste(ves_chr.df3$mz))
all_feat3 <- left_join(all_feat3, ves_chr.df3, by = "mz" ) %>% 
    select(all_of(c(names(all_feat3), "ion2")))
names(all_feat3)[names(all_feat3) == "ion2"] <- "annot"

#hyperborea
hyp_chr.df3 <- hyp_chr.df2 %>% distinct(ion2, .keep_all = TRUE)
hyp_chr.df3 <- hyp_chr.df3[!hyp_chr.df3$mz %in% ves_chr.df3$mz, ]
hyp_chr.df3$rt <- NULL
hyp_chr.df3$mz <- as.numeric(paste(hyp_chr.df3$mz))
all_feat3 <- left_join(all_feat3, hyp_chr.df3, by = "mz" ) %>% 
    select(all_of(c(names(all_feat3), "ion2")))
all_feat3 <- all_feat3 %>% 
    mutate(annot = coalesce(annot,ion2))
all_feat3$ion2 <- NULL

#digitata
dig_chr.df3 <- dig_chr.df2 %>% distinct(ion2, .keep_all = TRUE)
dig_chr.df3 <- dig_chr.df3[!dig_chr.df3$mz %in% c(ves_chr.df3$mz,
                                                  hyp_chr.df3$mz), ]
dig_chr.df3$rt <- NULL
dig_chr.df3$mz <- as.numeric(paste(dig_chr.df3$mz))
all_feat3 <- left_join(all_feat3, dig_chr.df3, by = "mz" ) %>% 
    select(all_of(c(names(all_feat3), "ion2")))
all_feat3 <- all_feat3 %>% 
    mutate(annot = coalesce(annot,ion2))
all_feat3$ion2 <- NULL

#bicyclis
bic_chr.df3 <- bic_chr.df2 %>% distinct(ion2, .keep_all = TRUE)
bic_chr.df3 <- bic_chr.df3[!bic_chr.df3$mz %in% c(ves_chr.df3$mz,
                                                  hyp_chr.df3$mz,
                                                  dig_chr.df3$mz), ]
bic_chr.df3$rt <- NULL
bic_chr.df3$mz <- as.numeric(paste(bic_chr.df3$mz))
all_feat3 <- left_join(all_feat3, bic_chr.df3, by = "mz" ) %>% 
    select(all_of(c(names(all_feat3), "ion2")))
all_feat3 <- all_feat3 %>% 
    mutate(annot = coalesce(annot,ion2))
all_feat3$ion2 <- NULL

#yeast
yea_chr.df3 <- yea_chr.df2 %>% distinct(ion2, .keep_all = TRUE)
yea_chr.df3 <- yea_chr.df3[!yea_chr.df3$mz %in% c(ves_chr.df3$mz,
                                                  hyp_chr.df3$mz,
                                                  dig_chr.df3$mz,
                                                  bic_chr.df3$mz), ]
yea_chr.df3$rt <- NULL
yea_chr.df3$mz <- as.numeric(paste(yea_chr.df3$mz))
all_feat3 <- left_join(all_feat3, yea_chr.df3, by = "mz" ) %>% 
    select(all_of(c(names(all_feat3), "ion2")))
all_feat3 <- all_feat3 %>% 
    mutate(annot = coalesce(annot,ion2))
all_feat3$ion2 <- NULL

#t. weissflogii
twa_chr.df3 <- twa_chr.df2 %>% distinct(ion2, .keep_all = TRUE)
twa_chr.df3 <- twa_chr.df3[!twa_chr.df3$mz %in% c(ves_chr.df3$mz,
                                                  hyp_chr.df3$mz,
                                                  dig_chr.df3$mz,
                                                  bic_chr.df3$mz,
                                                  yea_chr.df3$mz), ]
twa_chr.df3$rt <- NULL
twa_chr.df3$mz <- as.numeric(paste(twa_chr.df3$mz))
all_feat3 <- left_join(all_feat3, twa_chr.df3, by = "mz" ) %>% 
    select(all_of(c(names(all_feat3), "ion2")))
all_feat3 <- all_feat3 %>% 
    mutate(annot = coalesce(annot,ion2))
all_feat3$ion2 <- NULL

#add time to annotation
all_feat3$annot2 <- paste0(all_feat3$annot, " (", 
                           round(all_feat3$rt/60, 2),
                           " min)")
all_feat3$annot2 <- all_feat3$annot2 %>% 
    sub("=", " = ", .) %>% gsub("\\s{2}", " ", .)
all_feat4 <- all_feat3[order(all_feat3$annot2),]

#get intensity data and names only
all_feat4.int <- all_feat4 %>% 
    ungroup() %>% 
    select(ends_with("1") | ends_with("2") | ends_with("3") | "annot2")

#make long
all_feat4.int <- all_feat4.int %>% 
    pivot_longer(cols = !contains("annot"), names_to = "sample",
                 values_to = "intensity")

#add sample groups
all_feat4.int$group <- all_feat4.int$sample %>% 
    sub(".*yea\\d", "S. cerevisiae mannan", .) %>% 
    sub(".*twa\\d", "T. wiessflogii sulphated mannan", .) %>% 
    sub(".*bic\\d", "E. bicyclis laminarin", .) %>% 
    sub(".*dig\\d", "L. digitata laminarin", .) %>%
    sub(".*ves\\d", "F. vesiculosus fucoidan", .) %>%
    sub(".*hyp\\d", "L. hyperborea fucoidan", .) %>%
    sub(".*solvent.*", "solvent blank", .)

all_feat4.int.samples <- all_feat4.int %>% filter(group != "solvent blank")

all_feat4.int.samples$group <- factor(all_feat4.int.samples$group, 
                                      levels = c("E. bicyclis laminarin",
                                                 "L. digitata laminarin",
                                                 "F. vesiculosus fucoidan",
                                                 "L. hyperborea fucoidan",
                                                 "S. cerevisiae mannan",
                                                 "T. wiessflogii sulphated mannan"))

all_feat4.int.samples$runtype <- all_feat4.int.samples$sample %>% 
    sub("X20211129_.*", "full MS", .) %>% 
    sub("X20211216_.*", "MS/MS", .)

all_feat4.int.samples$intensity[all_feat4.int.samples$intensity == 0] <- NA
all_feat4.int.samples$unique <- ""

a <- all_feat4.int.samples$annot2 %>% unique()
for(i in 1:length(a)){
    df1 <- all_feat4.int.samples %>% filter(annot2 == !!a[i])
    df1 <- df1[complete.cases(df1),]
    df2 <- df1 %>% 
        group_by(group) %>% 
        dplyr::summarise(n = n())
    if(nrow(df2) > 1){
        all_feat4.int.samples$unique[all_feat4.int.samples$annot2 == a[i]] <- "no"
    } else if(nrow(df2) == 1){
        all_feat4.int.samples$unique[all_feat4.int.samples$annot2 == a[i]] <- "yes"
    }
}


all_feat4.int.samples$sample2 <- all_feat4.int.samples$sample %>% 
    sub("X20211129_|X20211216_", "", .)
all_feat4.int.samples$sample2  <- paste(all_feat4.int.samples$sample2, 
                                        all_feat4.int.samples$runtype)

p <- ggplot(all_feat4.int.samples, aes(x = sample2, y = annot2, size = intensity, 
                                       colour = unique, fill = group)) +
    geom_point(shape = 21) + labs(y = "Feature") +
    scale_size_continuous(name = "Integrated intensity") +
    scale_fill_manual(name = "",
                      values = hcl.colors(palette = "Dynamic", 6)) +
    scale_colour_manual(values = c("darkgrey", "black"),
                        guide = "none") +
    theme_bw() +
    theme(text = element_text(family = "Avenir"),legend.position = "top",
          axis.title.x = element_blank(), axis.text.y = element_blank(),
          axis.text.x = element_text(angle =  45, hjust =1))

png(filename = "analysis/analysis_plots/20220207_bubble-plot_v1.png",
    height = 7, width = 12, units = "in", res = 300)
print(p)
dev.off()

#18: Clustering-----
feat1 <- all_feat4 %>% 
    ungroup() %>% 
    select(ends_with("1") | ends_with("2") | ends_with("3") | "annot2") %>% 
    select(!contains("solvent"))
feat1 <- as.data.frame(feat1)
rownames(feat1) <- feat1$annot2
feat1$annot2 <- NULL

bray <- vegdist(t(feat1), method = "bray")
clust <- hclust(bray)
plot(clust)

library(ggdendro)
clust.gg <- dendro_data(clust)

clust.gg.label <- label(clust.gg)
clust.gg.label$group <- clust.gg.label$label %>% 
    sub(".*yea\\d", "S. cerevisiae mannan", .) %>% 
    sub(".*twa\\d", "T. wiessflogii sulphated mannan", .) %>% 
    sub(".*bic\\d", "E. bicyclis laminarin", .) %>% 
    sub(".*dig\\d", "L. digitata laminarin", .) %>%
    sub(".*ves\\d", "F. vesiculosus fucoidan", .) %>%
    sub(".*hyp\\d", "L. hyperborea fucoidan", .)

clust.gg.label$group  <- factor(clust.gg.label$group, 
                                levels = c("E. bicyclis laminarin",
                                           "L. digitata laminarin",
                                           "F. vesiculosus fucoidan",
                                           "L. hyperborea fucoidan",
                                           "S. cerevisiae mannan",
                                           "T. wiessflogii sulphated mannan"))

clust.gg.label$runtype <- clust.gg.label$label %>% 
    sub("X20211129_.*", "full MS", .) %>% 
    sub("X20211216_.*", "MS/MS", .)


#plot
p <- ggplot()+
    geom_segment(data = segment(clust.gg), 
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_point(data = clust.gg.label, aes(x = x, y = y,
                                          fill = group, shape = runtype),
               colour = "black", size = 3) +
    scale_fill_manual(name = "",
                      values = hcl.colors(palette = "Dynamic", 6)) + 
    scale_shape_manual(values = c(21, 24), name = "") +
    coord_flip() +
    theme_classic() +
    scale_y_reverse(expand = c(0.2, 0)) +
    guides(fill = guide_legend(override.aes=list(shape=21))) +
    theme(text = element_text(family = "Avenir"),
          legend.position = "top",
          axis.title = element_blank(),
          axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank())

png(filename = "analysis/analysis_plots/20220207_clustering_v1.png",
    height = 5, width = 8, units = "in", res = 300)
print(p)
dev.off()





