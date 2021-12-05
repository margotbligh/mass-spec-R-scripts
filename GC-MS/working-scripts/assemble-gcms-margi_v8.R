setwd("/Users/margotbligh/ownCloud/ASSEMBLE-margi")
load("Rdata.Rdata")
#1: Libraries
library(tidyverse)
library(dplyr)
library(data.table)
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(vegan)
library(grid)
library(ropls)
library(viridis)
library(pairwiseComparisons)
library(rstatix)
library(ggprism)
library(xcms)
library(CAMERA)
library(metaMS)
library(MSnbase)
library(scales)
library(scico)
library(FinCal)
library(gplots)
library(robustbase)
library(scatterplot3d)


#clean environment from RData loaded (v7 of script)
# rm(list=ls(pattern = "feat"))
# rm(list=ls(pattern = "pspectra.q"))
# rm(ps_nr143, no.q, q.simmat.1, )

font_family<-"sans"
font_size<-10
stats_text_size<-3.2

#shut up dplyr
options(dplyr.summarise.inform = FALSE)

#2: Import data ----
fp = list.files("mzML-files", recursive = TRUE, 
                full.names = TRUE, pattern = ".mzML")


  #get metadata----
samplenames <- basename(fp)

#make dataframe with most variables filled in
pd <-data.frame(SampleName = samplenames,
                Plant = samplenames %>% 
                    sub("^A_.*", "ALGAE", .) %>% 
                    sub("^S_.*", "SEAGRASS", .) %>% 
                    sub("^BLK_.*|^ASW_.*", "NA", .), 
                WaterBody = samplenames %>% 
                    sub("^BLK_.*", "MQ", .) %>% 
                    sub("^ASW_.*", "ASW", .) %>% 
                    sub("^[AS]_WC.*", "WC", .) %>% 
                    sub("^S_PW.*", "PW", .),
                SamplingTime = samplenames %>% 
                    sub("^[AS]_[PW][WC]\\d{1,2}[Il]_.*", "I", .) %>% 
                    sub("^[AS]_[PW][WC]\\d{1,2}E_.*", "E", .) %>% 
                    sub("^BLK_.*|^ASW_.*", "NA", .),
                Condition = samplenames %>% 
                    sub("^A_\\D{2}[123].*", "light", .) %>% 
                    sub("^A_\\D{2}[456].*", "dark", .) %>% 
                    sub("^BLK_.*", "NA", .) %>% 
                    sub("^ASW_blank.*", "blank", .) %>% 
                    sub("^ASW_met.*", "metabolitestd", .) %>% 
                    sub("^S_\\D{2}[123789]\\D.*", "light", .) %>% 
                    sub("^S_\\D{2}[456]\\D.*", "dark", .) %>% 
                    sub("^S_\\D{2}\\d{2}.*", "dark", .),
                Enrichment13C = samplenames %>% 
                    sub("^A_.*|^BLK_.*|^ASW_.*", "NO", .) %>% 
                    sub("^S_\\D{2}[123456]\\D_.*", "NO", .) %>% 
                    sub("^S_.*", "YES", .))

#create SampleGroup column
pd$SampleGroup <- NA
pd$SampleGroup[pd$Condition == "blank"] <- "blank"
pd$SampleGroup[pd$Condition == "metabolitestd"] <- "standard"
pd$SampleGroup[pd$WaterBody == "MQ"] <- "MQ"
pd$SampleGroup[pd$Plant == "ALGAE" & pd$SamplingTime == "I"] <- "algae seawater"
pd$SampleGroup[pd$Plant == "ALGAE" & pd$SamplingTime == "E" & 
                 pd$Condition == "light"] <- "algae light"
pd$SampleGroup[pd$Plant == "ALGAE" & pd$SamplingTime == "E" & 
                 pd$Condition == "dark"] <- "algae dark"
pd$SampleGroup[pd$Plant == "SEAGRASS" & pd$SamplingTime == "I" &
                 pd$WaterBody == "WC"] <- "seagrass seawater"
pd$SampleGroup[pd$Plant == "SEAGRASS" & pd$SamplingTime == "E" & 
                 pd$Condition == "light" & pd$WaterBody == "WC"] <- "seagrass light"
pd$SampleGroup[pd$Plant == "SEAGRASS" & pd$SamplingTime == "E" & 
                 pd$Condition == "dark" & pd$WaterBody == "WC"] <- "seagrass dark"
pd$SampleGroup[pd$Plant == "SEAGRASS" & pd$SamplingTime == "I" &
                 pd$WaterBody == "PW"] <- "seagrass porewater"
pd$SampleGroup[pd$Plant == "SEAGRASS" & pd$SamplingTime == "E" & 
                 pd$Condition == "light" & pd$WaterBody == "PW"] <- "seagrass porewater light"
pd$SampleGroup[pd$Plant == "SEAGRASS" & pd$SamplingTime == "E" & 
                 pd$Condition == "dark" & pd$WaterBody == "PW"] <- "seagrass porewater dark"

#with new lines for plots
pd$SampleGroup.n <- pd$SampleGroup
pd$SampleGroup.n <- pd$SampleGroup.n %>% gsub("\\s", "\n", .)


  #create data object----
raw_data <- readMSData(files = fp, 
                       pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk", 
                       msLevel. = 1, 
                       centroided. = TRUE) 
#3: Peak picking ----

# parameters that worked in terms of number of features for algae
cwp<-CentWaveParam()
cwp@peakwidth<-c(2,10)
cwp@ppm<-300
cwp@integrate<-1L
cwp@noise<-0
cwp@prefilter<-c(3,100)
cwp@fitgauss<-F
cwp@snthresh=10
cwp@mzdiff=-0.001
cwp@mzCenterFun=c("wMean")

data_peaks <- findChromPeaks(raw_data, param = cwp)


#4: Initial XCMS peak grouping -----
pdp <- PeakDensityParam(sampleGroups = data_peaks$SampleGroup2,
                        minFraction = 0.5, bw=2, binSize = 0.25)
grouped_peaks <- groupChromPeaks(data_peaks, param = pdp)

#5: Rt alignment ----
rtalign_params<-PeakGroupsParam(minFraction = 0.5, span = 0.4)
peaks_RTaligned <- adjustRtime(grouped_peaks, param = rtalign_params)


#6: Final XCMS peak grouping ----
pdp <- PeakDensityParam(sampleGroups = peaks_RTaligned$SampleGroup2,
                        minFraction = 0.5, binSize = 0.25, bw = 2)
grouped_peaks_RTaligned<-groupChromPeaks(peaks_RTaligned, 
                                         param = pdp)
#7: Fill peaks----
fpp <- FillChromPeaksParam(expandMz = 0.5, expandRt = 0.5, ppm = 20)
peaks_final <- fillChromPeaks(grouped_peaks_RTaligned, param = fpp)

#8: Create xsAnnotate object----
xset <- as(peaks_final, "xcmsSet")
sampnames(xset) <- pData(peaks_final)$SampleName
sampclass(xset) <- pData(peaks_final)$SampleGroup2
an <- xsAnnotate(xset)
#9: CAMERA peak grouping ----
#keep in mind that each PC group = one pseudospectra
#therefore - check that there are not multiple PC groups with the same rt...
#this happened in previous analysis attempt!

#peak grouping by rt
an.group <- groupFWHM(an, perfwhm = 3)

#find isotopic peaks
an.iso <- findIsotopes(an.group, mzabs = 0.01)

#re-group by correlation
#an.grp.corr <- groupCorr(an.iso, cor_eic_th=0.9, calcIso = T)

#10: Format conversion of sample psuedospectra ----
  #nested list (matrices) ----
#get intensities of samples as dataframe
allpks <- an.iso@groupInfo
allpks <- as.data.frame(allpks)
allpksInt <- select(allpks, contains("Copy"))

#filter pseudospectra to have minimum 5 associated features
pspectra <- an.iso@pspectra
npeaks <- sapply(pspectra, length)
pspectra <- pspectra[npeaks >= 5]

#this makes a nested list as is the output from the metaMS wrappers
#each sample is a single element - within that each pseudospetra is an element
spclist.by.sample <- vector("list", length(pd$SampleName))
for (i in 1:length(pd$SampleName)) {
    spclist.by.sample[[i]] <- lapply(pspectra, 
                                     function(x)
                                         cbind(mz = round(allpks[x,"mz"], 
                                                          digits = 1),
                                               intensity = allpksInt[x, i],
                                               rt = allpks$rt[x]))
}

#get retention time list
rt.spclist <- sapply(spclist.by.sample[[1]], "[", 1, 3) # subset from ps list, row 1, column 3 (contains rt of ps)
#rt.spclist.sort <- sort(rt.spclist)
#rt.spclist.sort <- floor(rt.spclist.sort)
#rt.spclist.sort <- unique(rt.spclist.sort)

#n >= 5 with an.grp.corr = 176
#only 168 when rounded to nearest second... these are too close!

#n >= 5 with an.group = 239
#also 239 when rounded to nearest second... they are all different times! SUCCESS

#use grouping only by retention time

rm(rt.spclist.sort)
rm(an.grp.corr, grouped_peaks, grouped_peaks_RTaligned, data_peaks, peaks_RTaligned)

  #nested list (spectra)-----
list1 <- vector("list", length(pspectra))
spclist.by.sample.spec <- vector("list", nrow(pd))
spclist.by.sample.spec[[1]] <- list1 
rm(list1)

for (i in 1:length(spclist.by.sample)){
    for (j in 1:length(rt.spclist)){
        ps <- spclist.by.sample[[i]][[j]]
        ps <- as.data.frame(ps)
        ps <- ps[complete.cases(ps),]
        ps <- as.matrix(ps)
        if (any(class(ps) == "numeric")){
            spclist.by.sample.spec[[i]][[j]] <- new("Spectrum1", centroided = T)
            next}
        if (nrow(ps) == 0){
            spclist.by.sample.spec[[i]][[j]] <- new("Spectrum1", centroided = T)
            next}
        else {
            spclist.by.sample.spec[[i]][[j]] <- new("Spectrum1",
                                                    mz = ps[,"mz"],
                                                    intensity = ps[,"intensity"],
                                                    centroided = T)}
    }
}
for(i in 1:length(spclist.by.sample.spec)){
    spclist.by.sample.spec[[i]] <- lapply(spclist.by.sample.spec[[i]], 
                                          clean, all = T)
}

spclist.by.sample.spec.norm <- list()
for(i in 1:length(spclist.by.sample.spec)){
    spclist.by.sample.spec.norm[[i]] <- lapply(spclist.by.sample.spec[[i]], 
                                               normalise, method = "max")
}

#11: Load reference spectra and standards table----
  #mass bank of north america (mona) ----
load("mona.spec.RData")
load("mona.names.RData")
load("mona.msp.RData")
#previous code to generate:

#mona.msp <- metaMS::read.msp(file = "MoNA-export-GC-MS_Spectra.msp")
#mona.names <- unlist(sapply(mona.msp, "[", 1))
#create spectrum objects for each db pseudospectrum
### CAREFUL: takes super long ###
#mona.spec <- lapply(mona.msp,
#                    function(x) new("Spectrum1", 
#                                    mz = x$pspectrum[,1], 
#                                    intensity = x$pspectrum[,2],
#                                    centroided = TRUE))

  #standards table formatted - from Marvin's table on confluence ----
mlist <- fread("MS57-salty-compounds_v3.csv")

#previous code used to generate (UGH):
# mlist <- fread("MS57-salty-compounds_v2.csv", na.strings="")
# names(mlist) <- c("compound", "top3ions","rt")
# mlist <- mlist[!is.na(mlist$top3ions)]
# mlist$index <- rownames(mlist)
# 
# #split out rows that need ions fixes
# mlist1 <- mlist
# mlist1 <- mlist1[grep(":", mlist1$top3ions)]
# 
# 
# #split out rows where ions are fine
# mlist2 <- mlist[!mlist$index %in% mlist1$index,]
# 
# #fix retention times for these
# mlist2$rt <- mlist2$rt %>% 
#     gsub("[[:alpha:]]", "", .) %>% 
#     gsub("\\s{2}.*", "", .) %>% 
#     sub("\\s5\\s,$", "", .) %>% 
#     sub(",$", "", .)
# mlist2 <- mlist2 %>% filter(rt != "")
# 
# #fix Adenosine -3´,5´ cyclophosphoric acid ions
# mlist2$top3ions[mlist2$top3ions== "Adenine"] <- mlist2$top3ions[mlist2$compound== "Adenine Hydrochloride"]
# 
# #split out top 3 ions
# mlist2$ion1 <- str_split_fixed(mlist2$top3ions, ",", 3)[,1] %>% as.numeric()
# mlist2$ion2 <- str_split_fixed(mlist2$top3ions, ",", 3)[,2] %>% as.numeric()
# mlist2$ion3 <- str_split_fixed(mlist2$top3ions, ",", 3)[,3] %>% as.numeric()
# 
# #split out retention times
# mlist2 <- cbind(mlist2, str_split(mlist2$rt, ",", simplify = T))
# names(mlist2) <- names(mlist2) %>% sub("V", "rt", .)
# mlist2 <- mlist2 %>% pivot_longer(., cols = grep("rt", names(mlist2)),
#                                   names_to = "k", values_to = "rt")
# mlist2$k <- NULL
# mlist2$rt <- as.numeric(mlist2$rt)
# mlist2 <- mlist2[complete.cases(mlist2),]
# mlist2 <- distinct(mlist2)
# mlist2$top3ions <- NULL
# 
# #deal with mlist1 - where no + is presenet
# mlist1.1 <- mlist1[!grepl("\\+", mlist1$top3ions),]
# mlist1.1.1 <- mlist1.1[grepl("\n", mlist1.1$top3ions),]
# mlist1.1.2 <- mlist1.1[!mlist1.1$index %in% mlist1.1.1$index,]
# 
# #with multiple lines
# mlist1.1.1 <- cbind(mlist1.1.1, 
#                     mlist1.1.1$top3ions %>% str_split(., "\n\n", simplify = T))
# names(mlist1.1.1) <- sub("V", "top3ions", names(mlist1.1.1))
# mlist1.1.1$top3ions1[mlist1.1.1$compound == "Fructose"] <- "18.34: 307.3, 217.2, 103.1"
# mlist1.1.1$top3ions2[mlist1.1.1$compound == "Fructose"] <- "18.44: 307.3, 217.2, 103.1"
# mlist1.1.1$top3ions <- NULL
# mlist1.1.1$rt <- NULL
# mlist1.1.1 <- mlist1.1.1 %>% pivot_longer(., 
#                                           cols = c("top3ions1", "top3ions2", "top3ions3"),
#                                           names_to = "k", values_to = "top3ions")
# mlist1.1.1$k <- NULL
# mlist1.1.1$rt <-mlist1.1.1$top3ions %>% sub(":.*", "", .) %>% as.numeric()
# mlist1.1.1 <- mlist1.1.1[complete.cases(mlist1.1.1),]
# mlist1.1.1$top3ions <- mlist1.1.1$top3ions %>% sub(".*:\\s*", "", .)
# mlist1.1.1 <- cbind(mlist1.1.1,
#                     str_split_fixed(mlist1.1.1$top3ions, ", ", n = 3))
# names(mlist1.1.1)[grep("^\\d", names(mlist1.1.1))] <- paste0("ion", 
#                                                              names(mlist1.1.1)[grep("^\\d", 
#                                                                                     names(mlist1.1.1))])
# mlist1.1.1$top3ions <- NULL
# mlist1.1.1 <- mlist1.1.1[names(mlist2)]
# 
# 
# #without multiple lines
# mlist1.1.2$rt <- mlist1.1.2$top3ions %>% sub(":.*", "", .) %>% as.numeric()
# mlist1.1.2$top3ions <- mlist1.1.2$top3ions %>% sub(".*:\\s*", "", .)
# mlist1.1.2 <- cbind(mlist1.1.2,
#                     str_split_fixed(mlist1.1.2$top3ions, ", ", n = 3))
# names(mlist1.1.2) <- names(mlist1.1.2) %>% sub("V", "ion", .)
# mlist1.1.2$top3ions <- NULL
# setDF(mlist1.1.2)
# mlist1.1.2 <- mlist1.1.2[names(mlist2)]
# 
# mlist1.1 <- rbind(mlist1.1.1, mlist1.1.2)
# 
# #deal with mlist1 - where + is presenet
# mlist1.2 <- mlist1[!mlist1$index %in% mlist1.1$index,]
# mlist1.2 <- cbind(mlist1.2,
#                   str_split(mlist1.2$top3ions, "\n\n", simplify = T))
# mlist1.2$V4 <- mlist1.2$V1 %>% sub(".*\\+", "", .)
# mlist1.2$V1 <- mlist1.2$V1 %>% sub("\\+.*:", ":", .)
# mlist1.2$top3ions <- NULL
# mlist1.2$rt <- NULL
# names(mlist1.2) <- sub("V", "top3ions", names(mlist1.2))
# mlist1.2 <- mlist1.2 %>% pivot_longer(., 
#                                       cols = grep("top3ions", names(mlist1.2)),
#                                       names_to = "k", values_to = "top3ions")
# mlist1.2$k <- NULL
# mlist1.2$rt <-mlist1.2$top3ions %>% sub(":.*", "", .) %>% as.numeric()
# mlist1.2 <- mlist1.2[complete.cases(mlist1.2),]
# mlist1.2$top3ions <- mlist1.2$top3ions %>% sub(".*:\\s*", "", .)
# mlist1.2 <- cbind(mlist1.2,
#                   str_split_fixed(mlist1.2$top3ions, ", ", n = 3))
# names(mlist1.2)[grep("^\\d", names(mlist1.2))] <- paste0("ion", 
#                                                          names(mlist1.2)[grep("^\\d", 
#                                                                               names(mlist1.2))])
# mlist1.2$top3ions <- NULL
# mlist1.2 <- mlist1.2[names(mlist2)]
# mlist1.2$ion2[mlist1.2$ion2 == "341,3"] <- "341.3"
# mlist1.2$ion3[mlist1.2$ion3 == "341,3"] <- "341.3"
# 
# 
# #bind all together
# mlist.3 <- rbind(mlist1.1, mlist1.2, mlist2)
# 
# #make ion values numeric
# x <- grep("ion", names(mlist.3))
# mlist.3[x] <- sapply(mlist.3[x],as.numeric)
# 
# #write to file and reload
# fwrite(mlist.3, "MS57-salty-compounds_v3.csv")
# rm(list=ls(pattern="^mlist"))
# mlist <- fread("MS57-salty-compounds_v3.csv")


#12: Find ribitol pseudospectra (RIBITOL = 23)----
#get ribitol ref spec
ref_spec_ribitol<-mona.spec[which(mona.names=="Ribitol")]

#get ribitol rt
ribitol_rt.m <- mlist$rt[mlist$compound == "Adonitol (=Ribitol/Xylitol)"] #min
ribitol_rt.s <- ribitol_rt.m * 60 #sec

#get ribitol mz (top 3 ions)
ribitol_mz <- c(mlist$ion1[mlist$compound == "Adonitol (=Ribitol/Xylitol)"],
                mlist$ion2[mlist$compound == "Adonitol (=Ribitol/Xylitol)"],
                mlist$ion3[mlist$compound == "Adonitol (=Ribitol/Xylitol)"])

#check ps by rt
names(rt.spclist) <- 1:length(pspectra)
ribitol.ps_nr <- names(rt.spclist)[rt.spclist > ribitol_rt.s-15 &
                                       rt.spclist < ribitol_rt.s+15] #5 options
ribitol.ps_nr <- as.numeric(ribitol.ps_nr)

#check for top 3 ions, exact matches of all 3 (sample 16 = standard)
for (i in 1:length(ribitol.ps_nr)){
    ps_nr <- ribitol.ps_nr[i]
    ps <- spclist.by.sample[[16]][[ps_nr]] %>% as.data.frame()
    if(all(ribitol_mz %in% ps$mz)){
        print(paste(ps_nr, "is a match"))
    }
    else{print(paste(ps_nr, "is NOT a match"))}
}

#only ps_nr 23 is a match! 
#check for top 3 ions, exact matches of any 3 (sample 16 = standard)
for (i in 1:length(ribitol.ps_nr)){
    ps_nr <- ribitol.ps_nr[i]
    ps <- spclist.by.sample[[16]][[ps_nr]] %>% as.data.frame()
    if(any(ribitol_mz %in% ps$mz)){
        print(paste(ps_nr, "is a match"))
    }
    else{print(paste(ps_nr, "is NOT a match"))}
}

#only ps_nr 23 is a match! looks promising 
#now final check against ref spec

for(i in 1:length(ref_spec_ribitol)){
    dot_prod = compareSpectra(ref_spec_ribitol[[i]],
                              spclist.by.sample.spec.norm[[16]][[23]],
                              fun = "dot")
    plot(ref_spec_ribitol[[i]],
         spclist.by.sample.spec.norm[[16]][[23]],     
         main = paste0("dotpr=", round(dot_prod, 3),
                       " for ref spec ", i))
}

#looks pretty good - I will assume that this ribitol

#13: Choose ribitol quant ion (mz = 160.2) ----
#get data for ps 23
names(spclist.by.sample) <- 1:nrow(pd)
for(i in 1:nrow(pd)){
  spclist.by.sample[[i]] <- lapply(spclist.by.sample[[i]], function(x) as.data.frame(x))
}
ribitol_all <- bind_rows(sapply(spclist.by.sample, "[", 23), .id = "sample_nr")
#add metadata
pd$sample_nr <- 1:nrow(pd)
ribitol_all$sample_nr <- as.numeric(ribitol_all$sample_nr)
ribitol_all <- dplyr::left_join(ribitol_all, pd)
#calculate mean, sd, and cv per group
ribitol_groups <- ribitol_all %>% 
  dplyr::group_by(mz, SampleGroup) %>% 
  dplyr::summarise(n = n(), mean = mean(intensity), sd = sd(intensity))
ribitol_groups$cv <- coefficient.variation(ribitol_groups$sd, ribitol_groups$mean)
#remove MQ
ribitol_groups<- ribitol_groups %>% filter(SampleGroup != "MQ")
#order first by mean (high to low) then CV (low to high)
ribitol_groups <- ribitol_groups[order(ribitol_groups$mean, decreasing = T),]
ribitol_groups <- ribitol_groups[order(ribitol_groups$cv),]
#exclude derivative ions
ribitol_groups <- ribitol_groups[!ribitol_groups$mz %in% c(147.1, 73.1, 74.1, 75.1),] 
#pick first ion per group
ribitol_groups.y <- ribitol_groups %>% dplyr::group_by(SampleGroup) %>% dplyr::slice(1)
ribitol_groups.y <- na.omit(ribitol_groups.y)
#get cv for all picked
ribitol_groups.yy <- ribitol_groups[ribitol_groups$mz %in% ribitol_groups.y$mz,]

#plot intensities for each
r <- unique(ribitol_groups.y$mz)
for(i in 1:length(r)){
  x <- ribitol_all[ribitol_all$mz == r[i],]
  
  p <- ggplot(data = x, aes(x=SampleGroup, y=intensity)) +
    geom_boxplot(fill = "#AAABAF",position = position_dodge(width = 0.8), 
                 alpha = 0.7) + 
    geom_point(fill = "#AAABAF", pch=21, size=2,position = position_dodge(width = 0.8), 
               alpha = 0.7) +
    scale_y_continuous(name=expression(Ribitol~intensity~(a.u.))) +
    theme_bw() + 
    theme(text=element_text(size=font_size, family = font_family, colour = "#262626"),
          axis.text = element_text(size=8, family = font_family, colour = "#262626",
                                   angle = 45, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "top",
          strip.background.x =element_blank())
  
  tiff(paste0("plots/ribitol/quant-ion-",r[i],".tiff"),
       height=12, width=30, units = "cm", res=300)
  print(p)
  dev.off()
}

#choose based on plots
ribitol_q <- 160.2
#check chromatograms on MPI servers

#14: Normalise each sample by intensity of ribitol quant ion ----
spclist.by.sample.raw <- spclist.by.sample
for(i in 1:length(spclist.by.sample)){
    x <- spclist.by.sample.raw[[i]] #subset
    r <- x[[23]]$intensity[x[[23]]$mz == ribitol_q] #ribitol quant ion intensity
    #get those with only 1 row (weirdly form 3x1 df)
    n <- which(sapply(spclist.by.sample.raw[[i]], function(y) is.null(nrow(y))))
    if(length(n) > 0){
        for(j in 1:length(n)){
            index = n[j]
            x[[index]] <- t(x[[index]]) #transpose
        }
    }
    x <- lapply(x, transform, norm_intensity = intensity/r) #normalise
    spclist.by.sample[[i]] <- x #add back in
}

#15: Relative quantification after ribitol normalisation by weighted regression, algae + stds + blanks only----
#adapted from metaMS source code, "relInt" function: https://github.com/yguitton/metaMS/blob/master/R/relInt.R
relInt <- function(spec, refspec) {
  common.masses <- spec[spec[, 1] %in% refspec[, 1], 1]
  sampI <- spec[spec[, 1] %in% refspec[, 1], 4]
  refI <- refspec[refspec[, 1] %in% spec[, 1], 4]
  
  if (length(sampI) == 0) {
    ## can happen when comparing unknowns...
    return(0)
  }
  
  if (length(sampI) == 1) {
    ## only one mass in common, should not happen very often
    sampI/refI
  } else {
    if (length(sampI) < 5) {
      ## here we should use weighted regression! 
      ## Agreement in high masses is much more important than agreement in low masses.
      relI <- coef(lm(sampI ~ refI, weights = sqrt(common.masses)))[2]
    } else {
      ## if the robust estimator breaks down, we use the regular one... 
      ## this sometimes happens when the two vectors are really very similar
      relI <- try(coef(lmrob(sampI ~ refI))[2], silent = TRUE)
      if (is(relI)[1] == "try-error")
        relI <- coef(lm(sampI ~ refI, weights = sqrt(common.masses)))[2]
    }
    relI
  }
}
removeDoubleMasses <- function(spec){
  if(length(spec$mz) == length(unique(spec$mz))){
    return(spec)
  }
  else if (length(spec$mz) != length(unique(spec$mz))){
    rt <- spec %>% distinct(mz, .keep_all = T)
    rt <- rt$rt
    spec <- spec %>% dplyr::group_by(mz) %>% 
      dplyr::summarise(intensity = mean(intensity), 
                       norm_intensity = mean(norm_intensity))
    spec <- as.data.frame(spec)
    spec <- cbind(spec, rt)
    spec <- spec[,c("mz", "intensity", "rt", "norm_intensity")]
  }
  spec
}
pd.alg <- pd[pd$Plant == "ALGAE" | pd$Condition == "blank" | pd$Condition == "metabolitestd" ,]
spclist.by.sample.alg <- spclist.by.sample[which(pd$SampleName %in% pd.alg$SampleName)]
names(spclist.by.sample.alg) <- which(pd$SampleName %in% pd.alg$SampleName)

pspectra.r.alg <- data.frame(sample_name = rep(pd.alg$SampleName, length(pspectra)),
                             sample_nr = rep(names(spclist.by.sample.alg), length(pspectra)),
                             ps_nr = rep(1:length(pspectra), each = nrow(pd.alg)),
                             rel_intensity = NA, ref_sample = NA)

std.n <- which(pd.alg$Condition == "metabolitestd")
blank.n <- which(pd.alg$Condition == "blank")
seawater.n  <- which(pd.alg$SamplingTime == "I")
inc.n <- which(pd.alg$SamplingTime == "E")

#use as a reference spectra the spectra that has all ions not NA
#priority: standards, blanks, samples
for (i in 1:length(pspectra)){
  x <- lapply(spclist.by.sample.alg, "[[", i)
  index <- vector()
  for (n in 1:length(std.n)){
    a <- x[[std.n[n]]]
    if(any(is.na(a))){
      next
    } else if(!any(is.na(a))){index <- c(index, std.n[n])
    }
  }
  if (length(index) == 1){
    refspec <- spclist.by.sample.alg[[index]][[i]]
    refspec <- removeDoubleMasses(refspec)
  } else if (length(index)==2){
    spec1 <- spclist.by.sample.alg[[index[1]]][[i]]
    spec1 <- removeDoubleMasses(spec1)
    spec2 <- spclist.by.sample.alg[[index[2]]][[i]]
    spec2 <- removeDoubleMasses(spec2)
    rel <- relInt(spec1, spec2)
    if (rel < 1) {
      index <- index[2]
      refspec <- spec2
    } else if (rel > 1) {
      index <- index[1]
      refspec <- spec1
    }
  } else if (length(index)==3){
    spec1 <- spclist.by.sample.alg[[index[1]]][[i]]
    spec1 <- removeDoubleMasses(spec1)
    spec2 <- spclist.by.sample.alg[[index[2]]][[i]]
    spec2 <- removeDoubleMasses(spec2)
    spec3 <- spclist.by.sample.alg[[index[3]]][[i]]
    spec3 <- removeDoubleMasses(spec3)
    rel1v2 <- relInt(spec1, spec2)
    rel1v3 <- relInt(spec1, spec3)
    rel2v3 <- relInt(spec2, spec3)
    if (rel1v2 > 1 & rel1v3 > 1) {
      index <- index[1]
      refspec <- spec1
    } else if (rel1v2 < 1 & rel2v3 > 1) {
      index <- index[2]
      refspec <- spec2
    } else if (rel1v3 < 1 & rel2v3 < 1) {
      index <- index[3]
      refspec <- spec3
    } else {
      index <- index[3]
      refspec <- spec3
    }
  }
  if (length(index) == 0){
    specs <- x[sapply(x, function(x) !any(is.na(x)))]
    if(length(specs) == 0){
      next
    }   else if(any(names(specs) %in% blank.n)){
      index <- names(specs)[names(specs) %in% blank.n] %>% as.numeric()
      index <- index[1]
      refspec <- spclist.by.sample.alg[[index]][[i]]
      refspec <- removeDoubleMasses(refspec)
    } else if (!any(names(specs) %in% blank.n)){
      if(any(names(specs) %in% seawater.n)){
        index <- names(specs)[names(specs) %in% seawater.n] %>% as.numeric()
        index <- index[1]
        refspec <- spclist.by.sample.alg[[index]][[i]]
        refspec <- removeDoubleMasses(refspec)
      } else if(!any(names(specs) %in% seawater.n)){
        index <- names(specs) %>% as.numeric()
        index <- index[1]
        refspec <- spclist.by.sample.alg[[index]][[i]]
        refspec <- removeDoubleMasses(refspec)
      }
    }
  }
  pspectra.r.alg$ref_sample[pspectra.r.alg$ps_nr == i] <- index
  for (j in 1:nrow(pd.alg)){
    s = pd.alg$SampleName[j]
    spec <- na.omit(spclist.by.sample.alg[[j]][[i]])
    spec <- removeDoubleMasses(spec)
    pspectra.r.alg$rel_intensity[pspectra.r.alg$ps_nr == i & 
                                   pspectra.r.alg$sample_name == s] <-  relInt(spec, refspec)
  }
}
names(pspectra.r.alg)[names(pspectra.r.alg) == "sample_name"] <- "SampleName"
pspectra.r.alg$sample_nr <- as.numeric(pspectra.r.alg$sample_nr)
pspectra.r.alg.met <- left_join(pspectra.r.alg, pd, by = c("SampleName", "sample_nr"))


#16: Plot intensities after relative quantification----
pspectra.r.alg.met$SampleGroup.n <- factor(pspectra.r.alg.met$SampleGroup.n,
                                          levels= c("blank", "standard",
                                                    "algae\nseawater",
                                                    "algae\ndark", "algae\nlight"))

pal.alg <-c("white", "grey", "#DEAF3A",
            scico(2, palette = "cork",begin=0.75,end=0.6, alpha=0.7))
names(pal.alg) <- levels(pspectra.r.alg.met$SampleGroup.n)

for(i in 1:length(pspectra)){
  x <- pspectra.r.alg.met %>% filter(ps_nr == i)
  if(any(is.na(x$ref_sample))){next}
  p <- ggplot(x, aes(x = SampleGroup.n, y = rel_intensity)) +
    geom_boxplot(mapping = aes(fill=SampleGroup.n), 
                 position = position_dodge(width = 0.8), alpha = 0.7) + 
    geom_point(mapping = aes(fill=SampleGroup.n), pch=21, size=2, 
               position = position_dodge(width = 0.8), alpha = 0.7) +
    scale_fill_manual(values = pal.alg, name = "")+
    labs(x = "", title = paste("Pseudospectra", i, "- intensities relative to",
                               pd.alg$SampleGroup[unique(x$ref_sample)],
                               "( sample number", unique(x$ref_sample), ")")) +
    theme_bw() + 
    theme(text=element_text(size=font_size, family = font_family, colour = "#262626"),
          axis.text = element_text(size=font_size, family = font_family, colour = "#262626"),
          panel.grid = element_blank(),
          plot.title = element_text(size=font_size, family = font_family, colour = "#262626",
                                    hjust = 0.5, face = "bold"),
          legend.position = "none")
  
  png(filename = paste0("plots/relative-quantification-algae/", i, ".png"),
      height = 12, width = 9, res = 300, units = "in")
  print(p)
  dev.off()
}

#17: Blank subtraction----
pspectra.r.alg.met.all <- pspectra.r.alg.met

#remove any spectra with no ref sample
pspectra.r.alg.met <- pspectra.r.alg.met.all[complete.cases(pspectra.r.alg.met.all),]
pspectra.r.alg.met$ps_nr %>% unique %>% length() #204
a <- pspectra.r.alg.met$ps_nr %>% unique
#do wilcox tests for all groups vs blanks (algae only)
#get column names for table
i = 1
ps_nr = a[i]
df1 <- pspectra.r.alg.met %>% filter(ps_nr == !!ps_nr)
df1 <- df1 %>% 
  rstatix::wilcox_test(rel_intensity ~ SampleGroup, ref.group = "blank",
                       alternative = "less") %>% 
  rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position()
cols <- cbind(ps_nr = NA, cols)
blank.wilcox.df <- cols[!any(is.na(cols)),]
rm(cols)

#do tests and p-value adjustment
for(i in 1:length(a)){
  ps_nr = a[i]
  if(ps_nr == 23){next} #ribitol
  else{
    df1 <- pspectra.r.alg.met %>% filter(ps_nr == !!ps_nr)
    df2 <- df1 %>% 
      rstatix::wilcox_test(rel_intensity ~ SampleGroup, ref.group = "blank",
                           alternative = "less") %>% 
      rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
      rstatix::add_significance(p.col = "p.adj") %>% 
      rstatix::add_xy_position()
    df2 <- cbind(ps_nr = ps_nr, df2)
    blank.wilcox.df <- rbind(blank.wilcox.df, df2)
  }
}

blank.wilcox.df.sig <- blank.wilcox.df %>% filter(p <= 0.05)
blank.wilcox.df.sig$ps_nr %>% unique() %>% length() #n = 127

#18: Significance testing (wilcox)----
#do wilcox tests for dark and light incubations vs seawater (algae only)
pspectra.r.alg.met.samples <- pspectra.r.alg.met %>% 
  filter(SampleGroup != "blank" & SampleGroup != "standard") %>% 
  filter(ps_nr %in% blank.wilcox.df.sig$ps_nr)
ps_nr.tests <- pspectra.r.alg.met.samples$ps_nr %>% unique()
#get column names for table
i = 1
ps_nr = ps_nr.tests[i]
df1 <- pspectra.r.alg.met.samples %>% filter(ps_nr == !!ps_nr)
df1 <- df1 %>% 
  rstatix::wilcox_test(rel_intensity ~ SampleGroup, ref.group = "algae seawater") %>% 
  rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position()
cols <- cbind(ps_nr = NA, cols)
wilcox.df <- cols[!any(is.na(cols)),]
rm(cols)


#do tests and p-value adjustment
for(i in 1:length(ps_nr.tests)){
  ps_nr = ps_nr.tests[i]
  if(ps_nr == 23){next} #ribitol
  else{
    df1 <- pspectra.r.alg.met.samples %>% filter(ps_nr == !!ps_nr)
    df2 <- df1 %>% 
      rstatix::wilcox_test(rel_intensity ~ SampleGroup, ref.group = "algae seawater") %>% 
      rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
      rstatix::add_significance(p.col = "p.adj") %>% 
      rstatix::add_xy_position()
    df2 <- cbind(ps_nr = ps_nr, df2)
    wilcox.df <- rbind(wilcox.df, df2)
  }
}

wilcox.df.sig <- wilcox.df %>% filter(p.adj <= 0.05)
wilcox.df.sig$ps_nr %>% unique() %>% length() #n = 15

#19: Plot intensities of ps with significant differences-----
#dataframe
pspectra.r.alg.sig <- pspectra.r.alg.met[pspectra.r.alg.met$ps_nr %in% wilcox.df.sig$ps_nr,]
pspectra.r.alg.sig <- pspectra.r.alg.sig %>% filter(SampleGroup != "standard")
pspectra.r.alg.sig$SampleGroup <- factor(pspectra.r.alg.sig$SampleGroup,
                                         levels = c("blank",  "algae seawater", 
                                                    "algae dark", "algae light"))
#to remove if lower than seawater or blank
ps_nr.no <- c(14, 25, 132)
pspectra.r.alg.sig <- pspectra.r.alg.sig[!pspectra.r.alg.sig$ps_nr %in% ps_nr.no,]
wilcox.df.sig <- wilcox.df.sig[!wilcox.df.sig$ps_nr %in% ps_nr.no,]

#palette
time_cat_5colors<-c("white", scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                    scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))

#significance dataframe
wilcox.df.sig$xmin <- 2
wilcox.df.sig$xmax[wilcox.df.sig$group2 == "algae dark"] <- 3
wilcox.df.sig$xmax[wilcox.df.sig$group2 == "algae light"] <- 4


#plot
ggplot(data = pspectra.r.alg.sig, aes(x=SampleGroup, y=rel_intensity)) +
  geom_boxplot(mapping = aes(fill=SampleGroup),
               position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(mapping = aes(fill=SampleGroup), pch=21, size=2, 
             position = position_dodge(width = 0.8), alpha = 0.7) +
  scale_fill_manual(values = time_cat_5colors, name = "") +
  scale_y_continuous(name=expression(Relative~ribitol~normalised~intensity~(a.u.))) +
  facet_wrap(~ps_nr, strip.position = "bottom", scales = "free_y") +
  labs(x = "Pseudospectra") +
  add_pvalue(wilcox.df.sig, xmin = "xmin", xmax = "xmax", 
             y.position = wilcox.df.sig$y.position,
             label = "p.adj.signif",
             label.size = 5, fontfamily = "Arial") +
  theme_bw() + 
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"),
        strip.text.x = element_text(size=font_size, family = font_family, colour = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text.y = element_text(size=font_size, family = font_family, colour = "#262626"),
        axis.title = element_text(size=font_size, family = font_family, colour = "#262626"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.x =element_blank())

png(filename = "relative-intensities-sig_v1", height = 9, width = 12, res = 600,
    units = "in")
print(p)
dev.off()


#20: Get putative annotations-----
#annotations based on standard compound table
algae.mlist1 <- list()
a <- pspectra.r.alg.sig$ps_nr %>% unique()
for (i in 1:length(a)){
  ps_nr = a[i]
  ps <- bind_rows(lapply(spclist.by.sample, "[[", ps_nr), .id = "sample_nr")
  ps$mzmin <- ps$mz - 0.1
  ps$mzmax <- ps$mz + 0.1
  rt = rt.spclist[ps_nr]/60
  x <- mlist %>% filter(between(rt, !!rt-0.25, !!rt+0.25))
  x <- x %>% pivot_longer(cols = starts_with("ion"),
                          names_to = "ion", values_to = "mz")
  x$mzmin <- x$mz-0.1
  x$mzmax <- x$mz+0.1
  x$index <- NULL
  setDT(ps);setDT(x)
  setkey(x, mzmin, mzmax)
  df1 <- foverlaps(ps,x)
  df1 <- na.omit(df1)
  algae.mlist1[[i]] <- df1 %>% dplyr::group_by(compound, ion) %>% 
    dplyr::summarise(n=n())
}
names(algae.mlist1) <- as.character(a)
algae.mlist1 <- bind_rows(algae.mlist1, .id = "ps_nr")
fwrite(algae.mlist1,file = "algae.mlist1.txt")
#annotations by comparison to mona db
algae.mona1 <- data.frame(ps_nr = NA, rt = NA,
                          annot1 = NA, dotprod1 = NA, 
                          annot2 = NA,dotprod2 = NA,
                          annot3 = NA,dotprod3 = NA, 
                          annot4 = NA, dotprod4 = NA,
                          annot5 = NA, dotprod5 = NA)
for (i in 1:length(a)){
  ps_nr = a[i]
  l <- lapply(mona.spec,
              function(x) compareSpectra(x,
                                         spclist.by.sample.spec.norm[[18]][[ps_nr]],
                                         fun = "dot"))
  l <- unlist(l)
  names(l) <- 1:length(l)
  l <- sort(l, decreasing = TRUE)
  tmp.df <- data.frame(ps_nr = ps_nr, rt = rt.spclist[ps_nr],
                       annot1 = mona.names[as.numeric(names(l)[1])], dotprod1 = l[1],
                       annot2 = mona.names[as.numeric(names(l)[2])], dotprod2 = l[2],
                       annot3 = mona.names[as.numeric(names(l)[3])], dotprod3 = l[3],
                       annot4 = mona.names[as.numeric(names(l)[4])], dotprod4 = l[4],
                       annot5 = mona.names[as.numeric(names(l)[5])], dotprod5 = l[5])
  algae.mona1 <- rbind(algae.mona1, tmp.df)
  
}
algae.mona1 <- na.omit(algae.mona1)
fwrite(algae.mona1,file = "algae.mona1.txt")

algae.mona2 <- algae.mona1 %>% 
  pivot_longer(cols = contains("annot"), names_to = "annot_n", values_to = "annot") %>% 
  pivot_longer(cols = contains("dotprod"), names_to = "dotprod_n", values_to = "dotprod")
algae.mona2$annot_n <- sub("annot", "", algae.mona2$annot_n) %>% as.numeric()
algae.mona2$dotprod_n <- sub("dotprod", "", algae.mona2$dotprod_n) %>% as.numeric()
algae.mona2 <- algae.mona2[algae.mona2$annot_n == algae.mona2$dotprod_n,]

fwrite(algae.mona2,file = "algae.mona2.txt")




#21: Re-plot after checking chromatograms on GC-MS computer-----
ps_nr.no <- c(ps_nr.no, 10, 177, 163, 177, 220, 224)
pspectra.r.alg.sig <- pspectra.r.alg.sig[!pspectra.r.alg.sig$ps_nr %in% ps_nr.no,]
wilcox.df.sig <- wilcox.df.sig[!wilcox.df.sig$ps_nr %in% ps_nr.no,]

pspectra.r.alg.sig$compound <- NA
pspectra.r.alg.sig$compound[pspectra.r.alg.sig$ps_nr == 3] <- "small organic acid-like"
pspectra.r.alg.sig$compound[pspectra.r.alg.sig$ps_nr == 13] <- "mannitol"
pspectra.r.alg.sig$compound[pspectra.r.alg.sig$ps_nr == 19] <- "carbamate-like"
pspectra.r.alg.sig$compound[pspectra.r.alg.sig$ps_nr == 20] <- "small aromatic-like"
pspectra.r.alg.sig$compound[pspectra.r.alg.sig$ps_nr == 44] <- "small aromatic-like"
pspectra.r.alg.sig$compound[pspectra.r.alg.sig$ps_nr == 117] <- "small organic acid-like"
pspectra.r.alg.sig$compound[pspectra.r.alg.sig$ps_nr == 136] <- "alcohol-like"


wilcox.df.sig$compound <- NA
wilcox.df.sig$compound[wilcox.df.sig$ps_nr == 3] <- "small organic acid-like"
wilcox.df.sig$compound[wilcox.df.sig$ps_nr == 13] <- "mannitol"
wilcox.df.sig$compound[wilcox.df.sig$ps_nr == 19] <- "carbamate-like"
wilcox.df.sig$compound[wilcox.df.sig$ps_nr == 20] <- "small aromatic-like"
wilcox.df.sig$compound[wilcox.df.sig$ps_nr == 44] <- "small aromatic-like"
wilcox.df.sig$compound[wilcox.df.sig$ps_nr == 117] <- "small organic acid-like"
wilcox.df.sig$compound[wilcox.df.sig$ps_nr == 136] <- "alcohol-like"

compound_names <- list(`3`="small organic acid-like","13"="mannitol","19"="carbamate-like",
  "20"="small aromatic-like", "44"="small aromatic-like", "117"="small organic acid-like",
  "136"="alcohol-like")
compound_names <- unlist(compound_names)

png(filename = "plots/gmeeting/gc-ms.plot.png", height = 6, width = 8, units = "in", res = 600)
ggplot(data = pspectra.r.alg.sig, aes(x=SampleGroup, y=rel_intensity)) +
  geom_boxplot(mapping = aes(fill=SampleGroup),
               position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(mapping = aes(fill=SampleGroup), pch=21, size=2, 
             position = position_dodge(width = 0.8), alpha = 0.7) +
  scale_fill_manual(values = time_cat_5colors, name = "") +
  scale_y_continuous(name=expression(Relative~ribitol-normalised~intensity~(a.u.))) +
  facet_wrap(~ps_nr, strip.position = "bottom", scales = "free_y", nrow = 2,
             labeller = as_labeller(compound_names)) +
  labs(x = "") +
  add_pvalue(wilcox.df.sig, xmin = "xmin", xmax = "xmax", 
             y.position = wilcox.df.sig$y.position,
             label = "p.adj.signif",
             label.size = 5, fontfamily = "Arial") +
  theme_bw() + 
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"),
        strip.text.x = element_text(size=font_size, family = font_family, colour = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text.y = element_text(size=font_size, family = font_family, colour = "#262626"),
        axis.title = element_text(size=font_size, family = font_family, colour = "#262626"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.x =element_blank())
dev.off()




#22: Plots for presentation-----
peaks_final.stds <- filterFile(peaks_final, file = c(16, 17, 18))
chr1 <- chromatogram(peaks_final.stds, mz = c(160.1, 160.3))
par(oma= c(0,0,0,0))

png(filename = "plots/gmeeting/peaks.png", height = 4, width = 4, units = "in", res = 300)
par(oma= c(0,0,0,0))
plot(chr1, xlim = c(1100,1500), main = "", xlab = "Retention time (s)", 
     ylab = "Intensity (a.u.)")
dev.off()

chr1.sht <- chromatogram(peaks_final.stds, mz = c(160.1, 160.3), rt = c(1100,1500))
pdp <- PeakDensityParam(sampleGroups = peaks_final.stds$SampleGroup,
                        minFraction = 0.5, bw=2, binSize = 0.25)

png(filename = "plots/gmeeting/groups.png", height = 4, width = 4, units = "in", res = 300)
plotChromPeakDensity(chr1.sht, param = pdp, peakPch = 16)
dev.off()

chr1.sht.all <- chromatogram(peaks_final, mz = c(160.1, 160.3), rt = c(1100,1500))

par(mfrow=c(1,1), oma=c(2,0,0,0))
png(filename = "plots/gmeeting/rtime.png", height = 4, width = 4, units = "in", res = 300)
plotAdjustedRtime(peaks_final)
dev.off()

png(filename = "plots/gmeeting/peaks-filled.png", height = 4, width = 4, units = "in", res = 300)
plot(chr1.sht.all[,c(16:18)], xlim = c(1100,1500), main = "", xlab = "Retention time (s)", 
     ylab = "Intensity (a.u.)")
dev.off()

chr.2 <- chromatogram(peaks_final.stds, mz = c(217.1, 217.3), rt = c(900,1100))
chr.3 <- chromatogram(peaks_final.stds, mz = c(205.1, 205.3), rt = c(900,1100))
chr.4 <- chromatogram(peaks_final.stds, mz = c(319.2, 319.4), rt = c(900,1100))


rib <- list(chr1, chr.2, chr.3, chr.4)
rib.df <- data.frame(sample_nr = as.numeric(), rt = as.numeric(), 
                     mz = as.numeric(), intensity = as.numeric())
mz.rib <- c(160.2, 217.2, 205.2, 319.3)

for(i in 1:4){
  x <- rib[[i]]
  rt <- unlist(lapply(x, function(x) x@rtime))
  mz <- rep(mz.rib[i], length(rt))
  intensity <- unlist(lapply(x, function(x) x@intensity))
  sample_nr <- rep(c(16,17,18), each = length(rt)/3)
  tmp.df <- data.frame(sample_nr = sample_nr, rt = rt, 
                       mz = mz, intensity = intensity)
  rib.df <- rbind(rib.df, tmp.df)
}

png(filename = "plots/gmeeting/psuedospec.png", height = 4, width = 4, units = "in", res = 300)
plot(rib.df$rt[rib.df$mz == 160.2], log2(rib.df$intensity[rib.df$mz == 160.2]),
     type = "l", xlim = c(1000,1020), ylim = c(0,25), col = "#868451", lwd = 2,
     xlab = "Retention time (s)", ylab = "Log(Intensity (a.u.))")
lines(rib.df$rt[rib.df$mz == 217.2], log2(rib.df$intensity[rib.df$mz == 217.2]),
     type = "l", col = "#DAAC3A", lwd = 2)
lines(rib.df$rt[rib.df$mz == 319.3], log2(rib.df$intensity[rib.df$mz == 319.3]),
      type = "l", col = "#2BB6AF", lwd = 2)
dev.off()

png(filename = "plots/gmeeting/psuedospec.spec.png", height = 4, width = 4, units = "in", res = 300)
ggplot(spclist.by.sample[[16]][[13]]) +
  geom_segment(aes(y = 0, yend = intensity, x = mz, xend = mz)) +
  labs(y = "Intensity (a.u.)", x = expression(italic(m/z))) +
  theme_classic()
dev.off()


refspec = spclist.by.sample[[7]][[44]]
spec = spclist.by.sample[[9]][[44]]
common.masses <- spec[spec[, 1] %in% refspec[, 1], 1]
sampI <- spec[spec[, 1] %in% refspec[, 1], 4]
refI <- refspec[refspec[, 1] %in% spec[, 1], 4]

png(filename = "plots/gmeeting/lm.png", height = 4, width = 4, units = "in", res = 300)
plot(sampI, refI, xlab = "Intensity of common ions in sample spectrum",
     ylab = "Intensity of common ions in reference spectrum")
abline(lmrob(sampI ~ refI))
text(x = 1, y = 4, labels = "coefficient = relative intensity = 0.8333525")
dev.off()


blank.wilcox.df.sig$xmin <- 1
blank.wilcox.df.sig$xmax[blank.wilcox.df.sig$group2 == "algae seawater"] <- 2
blank.wilcox.df.sig$xmax[blank.wilcox.df.sig$group2 == "algae dark"] <- 3
blank.wilcox.df.sig$xmax[blank.wilcox.df.sig$group2 == "algae light"] <- 4

png(filename = "plots/gmeeting/blank.png", height = 4, width = 4, units = "in", res = 300)
ggplot(data = pspectra.r.alg.met[pspectra.r.alg.met$ps_nr == 135 &
                                   pspectra.r.alg.met$SampleGroup != "standard",], 
       aes(x=SampleGroup, y=rel_intensity)) +
  geom_boxplot(mapping = aes(fill=SampleGroup),
               position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(mapping = aes(fill=SampleGroup), pch=21, size=2, 
             position = position_dodge(width = 0.8), alpha = 0.7) +
  scale_fill_manual(values = time_cat_5colors, name = "",
                    labels = c("blank", "seawater", "algae dark", "algae light")) +
  scale_y_continuous(name=expression(Relative~ribitol-normalised~intensity~(a.u.))) +
  add_pvalue(blank.wilcox.df.sig[blank.wilcox.df.sig$group2 != "stanbdard" &
                                   blank.wilcox.df.sig$ps_nr == 35,], 
             xmin = "xmin", xmax = "xmax", 
             y.position = 2,
             label = "p.adj.signif",
             label.size = 5, fontfamily = "Arial") +
  labs(x = "Pseudospectra") +
  theme_bw() + 
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"),
        strip.text.x = element_text(size=font_size, family = font_family, colour = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text.y = element_text(size=font_size, family = font_family, colour = "#262626"),
        axis.title = element_text(size=font_size, family = font_family, colour = "#262626"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.x =element_blank())
dev.off()

png(filename = "plots/gmeeting/hypothesis.png", height = 4, width = 4, units = "in", res = 300)
ggplot(data = pspectra.r.alg.sig[pspectra.r.alg.sig$ps_nr == 117,], 
       aes(x=SampleGroup, y=rel_intensity)) +
  geom_boxplot(mapping = aes(fill=SampleGroup),
               position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(mapping = aes(fill=SampleGroup), pch=21, size=2, 
             position = position_dodge(width = 0.8), alpha = 0.7) +
  scale_fill_manual(values = time_cat_5colors, name = "") +
  scale_y_continuous(name=expression(Relative~ribitol-normalised~intensity~(a.u.))) +
  facet_wrap(~ps_nr, strip.position = "bottom", scales = "free_y", nrow = 2) +
  labs(x = "Pseudospectra") +
  add_pvalue(wilcox.df.sig[wilcox.df.sig$ps_nr == 117,], xmin = "xmin", xmax = "xmax", 
             y.position = wilcox.df.sig$y.position[wilcox.df.sig$ps_nr == 117],
             label = "p.adj.signif",
             label.size = 5, fontfamily = "Arial") +
  theme_bw() + 
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"),
        strip.text.x = element_text(size=font_size, family = font_family, colour = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text.y = element_text(size=font_size, family = font_family, colour = "#262626"),
        axis.title = element_text(size=font_size, family = font_family, colour = "#262626"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.x =element_blank())
dev.off()

png(filename = "plots/gmeeting/annot.png", height = 4, width = 4, units = "in", res = 300)
plot(spclist.by.sample.spec.norm[[18]][[23]],
     spclist.by.sample.spec.norm[[17]][[23]])
dev.off()


mannitol_3d_plot.df <- dplyr::bind_rows(lapply(spclist.by.sample.alg, "[[", 13), 
                                        .id = "sample_nr") 
xpositions <- c(5.4, 2.8, 5.6, 3, 5.8, 3.2, 4.4, 3.4, 4.6, 3.6, 4.8, 3.8,
                0.8,1,1.2, 1.8, 2, 2.2)
mannitol_3d_plot.df$xposition <- NA
#x = intensity
#y = sample
#z = mz

#blank = 13, 14, 15 <- 0.8,1,1.2
#standard = 16, 17, 18 <- 1.8, 2, 2.2
#seawater = 2,4,6,8,10,12 <- 2.8,3,3.2,3.4,3.6,3.8
#dark = 7,9,11 <- 4.4,4.6,4.8
#light=1,3,5 <- 5.4, 5.6, 5.8

for (i in 1:18){
  mannitol_3d_plot.df$xposition[mannitol_3d_plot.df$sample_nr == i] <- xpositions[i]
}

scatterplot3d(x = mannitol_3d_plot.df$xposition, y = mannitol_3d_plot.df$mz, 
              z =  mannitol_3d_plot.df$intensity, type = "h",
              pch = "", grid = F, angle=50)

mannitol_3d_plot.df$sample_group <- pd$SampleGroup[match(mannitol_3d_plot.df$sample_nr, pd$sample_nr)]
mannitol_3d_plot.df$sample_group <- sub("algae\\s", "", mannitol_3d_plot.df$sample_group)
mannitol_3d_plot.df$sample_group <- factor(mannitol_3d_plot.df$sample_group,
                                           levels = c("blank", "standard", "seawater",
                                                      "dark", "light"))
colors_3dplot <- c("grey", "black", time_cat_5colors[2:4])
colors_3dplot <- colors_3dplot[as.numeric(mannitol_3d_plot.df$sample_group)]


source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

png("plots/gmeeting/mannitol_3d.png", height = 6, width = 7, units = "in",
    res = 600)
png("plots/gmeeting/mannitol_3d.png", height = 6, width = 7, units = "in",
    res = 600)
scatterplot3d(x = mannitol_3d_plot.df$xposition, y = mannitol_3d_plot.df$mz, 
              z =  sqrt(mannitol_3d_plot.df$intensity), type = "h",
              pch = "", grid = F, color = colors_3dplot, box = F, xlab = "",
              zlab = expression(sqrt(Intensity~(a.u.))), ylab = expression(italic(m/z)),
              x.ticklabs = c("", "blank", "standard", "seawater", "", "dark",
                             "light"), lwd = 1.3)
# addgrids3d(x = mannitol_3d_plot.df$xposition, y = mannitol_3d_plot.df$mz, 
#            z =  mannitol_3d_plot.df$intensity, grid = c("xy", "xz", "yz"))
dev.off()



#23: FIGURE FOR PAPER #1-------
# Fuc_O2 <- ggplot(data = Fucus)+
#   geom_boxplot(aes(x = Treatment, y = O2, fill=Treatment), position = position_dodge(width = 1),
#                color="#262626")+
#   geom_point(aes(x=Treatment, y=O2, group=Treatment, fill=Treatment),
#              pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
#   # annotate("text", x=0.5, y=max(Fucus$O2), label = "italic(Fucus~vesiculosus)", parse=T, vjust = 1, hjust = 0,size=3)+
#   scale_fill_manual(name="", labels=c("initial", "dark", "light"),
#                     values = time_cat_colors)+
#   scale_y_continuous(name=expression('O'[2]~(mu*'mol'~L^-1)))+
#   scale_x_discrete(name="")+
#   theme_bw()+  
#   theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
#         strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
#         legend.text = element_text(size=font_size, family = font_family),
#         axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
#         panel.grid = element_blank(),
#         legend.position = "none",
#         strip.background.y = element_rect(fill = "#ECECEC"),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())+
#   labs(title = expression(O[2]))

compound_names <- list(`3`="Small organic acid-like","13"="Mannitol","19"="Carbamate-like",
                       "20"="Small aromatic-like", "44"="Small aromatic-like", 
                       "117"="Small organic acid-like","136"="Alcohol-like")
compound_names <- unlist(compound_names)

p <- ggplot(data = pspectra.r.alg.sig)+
  geom_boxplot(aes(x = SampleGroup, y = rel_intensity, fill=SampleGroup), 
               position = position_dodge(width = 1), color="#262626",
               lwd = 1.5)+
  geom_point(aes(x=SampleGroup, y=rel_intensity, group=SampleGroup, fill=SampleGroup),
             pch=21, size=1, stroke = 1.5)+
  scale_fill_manual(values = time_cat_5colors, name = "",
                    labels = c("buffer blank", "before incubation",
                               "dark incubation", "light incubation")) +
  scale_y_continuous(name="Relative ribitol-normalized intensity (a.u.)") +
  scale_x_discrete(name="")+
  facet_wrap(~ps_nr, scales = "free_y", nrow = 1,
             labeller = as_labeller(compound_names)) +
  add_pvalue(wilcox.df.sig, xmin = "xmin", xmax = "xmax", 
             y.position = wilcox.df.sig$y.position,
             label = "p.adj.signif",
             bracket.colour = "#262626", bracket.size = 1.5,
             label.size = 5, fontfamily = font_family) +
  theme_bw()+
  theme(text=element_text(size=6.1, family = font_family, colour = "#262626"),
        strip.text.x = element_text(size=6.1, family = font_family, color = "#262626",
                                    hjust = 0),
        legend.text = element_text(size=6.1, family = font_family),
        axis.text = element_text(size=6.1, family=font_family, colour = "#262626"),
        panel.grid = element_blank(), panel.border = element_blank(),
        axis.line.x = element_blank(), 
        axis.line.y = element_line(colour = "#262626", size =1.5), 
        axis.ticks.y = element_line(colour = "#262626", size =1.5),
        legend.position = "none",
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


pdf(file = "figs/gc-ms_v1.pdf", width = 8.1141732, height = 2.106299)
p
dev.off()

  

p <- ggplot(data = pspectra.r.alg.sig)+
  geom_boxplot(aes(x = SampleGroup, y = rel_intensity, fill=SampleGroup), 
               position = position_dodge(width = 1), color="#262626",
               lwd = 0.7)+
  geom_point(aes(x=SampleGroup, y=rel_intensity, group=SampleGroup, fill=SampleGroup),
             pch=21, size=1, stroke = 0.7)+
  scale_fill_manual(values = time_cat_5colors, name = "",
                    labels = c("buffer blank", "before incubation",
                               "dark incubation", "light incubation")) +
  scale_y_continuous(name="Relative ribitol-normalized intensity (a.u.)") +
  scale_x_discrete(name="")+
  facet_wrap(~ps_nr, scales = "free_y", nrow = 1,
             labeller = as_labeller(compound_names)) +
  add_pvalue(wilcox.df.sig, xmin = "xmin", xmax = "xmax", 
             y.position = wilcox.df.sig$y.position,
             label = "p.adj.signif",
             bracket.colour = "#262626", bracket.size = 0.7,
             label.size = 8, fontfamily = font_family) +
  theme_bw()+
  theme(text=element_text(size=8, family = font_family, colour = "#262626"),
        strip.text.x = element_text(size=8, family = font_family, color = "#262626",
                                    hjust = 0),
        legend.text = element_text(size=8, family = font_family),
        axis.text = element_text(size=8, family=font_family, colour = "#262626"),
        panel.grid = element_blank(), panel.border = element_blank(),
        axis.line.x = element_blank(), 
        axis.line.y = element_line(colour = "#262626", size =0.7), 
        axis.ticks.y = element_line(colour = "#262626", size =0.7),
        legend.position = "top",
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


png(file = "figs/gc-ms_v2.png", width = 10, height = 3, units = "in",
    res = 300)
p
dev.off()

#24: FIGURE FOR PAPER #2----
pspectra.sig <- pspectra.r.alg.met %>% filter(ps_nr %in% pspectra.r.alg.sig$ps_nr)
pspectra.sig$SampleGroup <- factor(pspectra.sig$SampleGroup,
                                   levels = c("blank",  "standard",
                                              "algae seawater", 
                                              "algae dark", "algae light"))


pspectra.sig.sum <- pspectra.sig %>% 
  dplyr::group_by(ps_nr, ref_sample, SampleGroup) %>% 
  dplyr::summarise(mean_rel_intensity = mean(rel_intensity))

time_cat_5colors<-c("white", "#FEC00099",
                    scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                    scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))
show_col(time_cat_5colors)

ggplot(data = pspectra.sig[pspectra.sig$ps_nr == 13,])+
  geom_boxplot(aes(x = SampleGroup, y = rel_intensity, fill=SampleGroup), 
               position = position_dodge(width = 1), color="#262626",
               lwd = 0.7)+
  geom_point(aes(x=SampleGroup, y=rel_intensity, group=SampleGroup, fill=SampleGroup),
             pch=21, size=1, stroke = 0.7)+
  scale_fill_manual(values = time_cat_5colors, name = "",
                    labels = c("buffer blank", "standard", "before incubation",
                               "dark incubation", "light incubation")) +
  scale_y_continuous(name="Relative ribitol-normalized intensity (a.u.)") +
  scale_x_discrete(name="")+
  facet_wrap(~ps_nr, scales = "free_y", nrow = 1,
             labeller = as_labeller(compound_names)) +
  add_pvalue(wilcox.df.sig[wilcox.df.sig$ps_nr == 13,], xmin = "xmin", xmax = "xmax", 
             y.position = wilcox.df.sig$y.position[wilcox.df.sig$ps_nr == 13],
             label = "p.adj.signif",
             bracket.colour = "#262626", bracket.size = 0.7,
             label.size = 8, fontfamily = font_family) +
  theme_bw()+
  theme(text=element_text(size=8, family = font_family, colour = "#262626"),
        strip.text.x = element_text(size=8, family = font_family, color = "#262626",
                                    hjust = 0),
        legend.text = element_text(size=8, family = font_family),
        axis.text = element_text(size=8, family=font_family, colour = "#262626"),
        panel.grid = element_blank(), panel.border = element_blank(),
        axis.line.x = element_blank(), 
        axis.line.y = element_line(colour = "#262626", size =0.7), 
        axis.ticks.y = element_line(colour = "#262626", size =0.7),
        legend.position = "top",
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

mannitol_3d_plot.df2 <- dplyr::bind_rows(lapply(spclist.by.sample.alg, "[[", 13), 
                                        .id = "sample_nr") 
xpositions <- c(5.4, 2.8, 5.6, 3, 5.8, 3.2, 4.4, 3.4, 4.6, 3.6, 4.8, 3.8,
                0.8,1,1.2, 1.8, 2, 2.2)
mannitol_3d_plot.df$xposition <- NA
#x = intensity
#y = sample
#z = mz

#blank = 13, 14, 15 <- 0.8,1,1.2
#standard = 16, 17, 18 <- 1.8, 2, 2.2
#seawater = 2,4,6,8,10,12 <- 2.8,3,3.2,3.4,3.6,3.8
#dark = 7,9,11 <- 4.4,4.6,4.8
#light=1,3,5 <- 5.4, 5.6, 5.8

for (i in 1:18){
  mannitol_3d_plot.df$xposition[mannitol_3d_plot.df$sample_nr == i] <- xpositions[i]
}

scatterplot3d(x = mannitol_3d_plot.df$xposition, y = mannitol_3d_plot.df$mz, 
              z =  mannitol_3d_plot.df$intensity, type = "h",
              pch = "", grid = F, angle=50)

mannitol_3d_plot.df$sample_group <- pd$SampleGroup[match(mannitol_3d_plot.df$sample_nr, pd$sample_nr)]
mannitol_3d_plot.df$sample_group <- sub("algae\\s", "", mannitol_3d_plot.df$sample_group)
mannitol_3d_plot.df$sample_group <- factor(mannitol_3d_plot.df$sample_group,
                                           levels = c("blank", "standard", "seawater",
                                                      "dark", "light"))
colors_3dplot <- c("grey", "black", time_cat_5colors[2:4])
colors_3dplot <- colors_3dplot[as.numeric(mannitol_3d_plot.df$sample_group)]


source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

par(family = "Helvetica")

png("plots/gmeeting/mannitol_3d.png", height = 6, width = 7, units = "in",
    res = 600)
svg("figs/mannitol_3d_v1.svg", height = 6, width = 7)
scatterplot3d(x = mannitol_3d_plot.df$xposition, y = mannitol_3d_plot.df$mz, 
              z =  sqrt(mannitol_3d_plot.df$norm_intensity), type = "h",
              pch = "", grid = F, color = colors_3dplot, box = F, xlab = "",
              zlab = expression(sqrt(Ribitol-normalized~untensity~(a.u.))), 
              ylab = expression(italic(m/z)), xlim = c(0, 6.6), ylim = c(0,600),
              x.ticklabs = c("", "blank", "standard", "before incubation", "", 
                             "dark incubation",
                             "light incubation"), lwd = 1)
dev.off()


mannitol_plot.df <- dplyr::bind_rows(lapply(spclist.by.sample.alg, "[[", 13), 
                                     .id = "sample_nr") 
mannitol_plot.df$sample_group <- pd$SampleGroup[match(mannitol_plot.df$sample_nr, 
                                                      pd$sample_nr)]
mannitol_plot.df$sample_group <- sub("algae\\s", "", mannitol_plot.df$sample_group)
mannitol_plot.df$sample_group <- factor(mannitol_plot.df$sample_group,
                                           levels = c("blank", "standard", "seawater",
                                                      "dark", "light"))
mannitol_plot.df$norm_intensity[is.na(mannitol_plot.df$norm_intensity)] <- 0
mannitol_plot.df <- mannitol_plot.df[order(mannitol_plot.df$mz),]




