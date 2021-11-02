setwd("/Users/margotbligh/Google_Drive/MPI_Masters/Assemble")
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
pd$SampleGroup<-paste(pd$WaterBody,
                      pd$Plant,
                      pd$Condition, 
                      pd$SamplingTime,
                      sep="_")
pd$SampleGroup2 <- paste(pd$WaterBody,
                         pd$Plant,
                         pd$Condition, 
                         sep="_")

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
##mass bank of north america (mona) ----
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



#13: Normalise all ions by ribitol quant ion intensity (per sample)----
#choose which ion to use as ribitol quant ion (highest intensity across most samples)
x <- lapply(spclist.by.sample, "[[", 23)
x <- lapply(x, function(x) as.data.frame(x))
x <- lapply(x, function(x) replace_na(x, list(mz = 0, intensity = 0, rt = 0)))
y <- sapply(x, function(x) x$mz[x$intensity == max(x$intensity)])
ribitol_q <- names(table(y))[table(y) == max(table(y))] %>% as.numeric()

#normalise each sample by intensity of ribitol quant ion
spclist.by.sample.raw <- spclist.by.sample
for(i in 1:length(spclist.by.sample)){
    x <- spclist.by.sample.raw[[i]] #subset
    x <- lapply(x, function(x) as.data.frame(x)) #make df
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

#14: Get quant ions for all pseudospec----
#set up df
pspectra.q.ions <- data.frame(ps_nr = seq(1, length(pspectra),1),
                              quant_mz = NA, max_sample = NA)
pspectra.q <- data.frame(sample_name = rep(names(allpksInt), length(pspectra)),
                         sample_nr = rep(1:nrow(pd), length(pspectra)),
                         ps_nr = rep(1:length(pspectra), each = nrow(pd)),
                         quant_mz = NA, intensity = NA, norm_intensity = NA)
pname <- names(pspectra.q)

#get quant ions
for(i in 1:length(pspectra)){
    #subset list and make df with no NA
    x <- lapply(spclist.by.sample, "[[", i)
    x <- lapply(x, function(x) as.data.frame(x))
    x <- lapply(x, function(x) drop_na(x))
    #remove df with only NA
    n <- sapply(x, function(x) nrow(x))
    x <- x[which(n>0)]

    #get max ion for each sample and find most common (= quant ion)
    y <- sapply(x, function(x) x$mz[x$norm_intensity == max(x$norm_intensity)])
    q <- names(table(y))[table(y) == max(table(y))] %>% as.numeric()
    pspectra.q.ions$quant_mz[i] <- q
    pspectra.q$quant_mz[pspectra.q$ps_nr == i] <- q
    #find sample with max intensity of quant ion
    x <- bind_rows(lapply(spclist.by.sample, "[[", i), .id = "sample_nr")
    x <- x[x$mz == q,]
    x <- x[complete.cases(x),]
    pspectra.q.ions$max_sample[i] <- x$sample_nr[
        x$norm_intensity == max(x$norm_intensity)]
    #add into table
    x$quant_mz <- q
    x$ps_nr <- i
    x$sample_nr <- as.numeric(x$sample_nr)
    x <- x[,c("sample_nr","ps_nr","quant_mz", "intensity", "norm_intensity")]
    x$sample_nr <- as.numeric(x$sample_nr)
    pspectra.q <- left_join(pspectra.q, x, by = c("sample_nr", "quant_mz", "ps_nr")) %>% 
        mutate(intensity = coalesce(intensity.x, intensity.y),
               norm_intensity = coalesce(norm_intensity.x, norm_intensity.y)) %>% 
        dplyr::select(all_of(pname))
}

#the max sample number for most is 23 or 29... which are blanks???
#but ribitol was not added to blanks, so the normalised intensity is
#of course super high!
#actually I won't worry about this :)

#format sample name to "..." format created by CAMERA
pd$SampleName_original <- pd$SampleName
pd$sample_name <- pd$SampleName_original %>% sub("\\s-\\s", "...", .)

#add metadata to quantification table
pspectra.q.met <- left_join(pspectra.q, pd, by = "sample_name")

#remove MQ (normalisation gives crazy values because no ribitol added)
pspectra.q.met.wMQ <- pspectra.q.met
pspectra.q.met <- pspectra.q.met %>% filter(SampleGroup != "MQ_NA_NA_NA")

#15: Significance testing----
#T-tests within groups----
#set up grouping variables
pspectra.q.met$SampleGroup3 <- pspectra.q.met$SampleGroup2 %>% 
    gsub("_dark|_light|_NA", "", .)
pspectra.q.met$SampleGroup <- pspectra.q.met$SampleGroup3 %>% 
    sub("_I|_E|_blank|_metab.*", "", .) 
pspectra.q.met$SamplingTime[pspectra.q.met$SampleGroup3=="ASW_blank"] <- "I"
pspectra.q.met$SamplingTime[pspectra.q.met$SampleGroup3=="ASW_metabolitestd"] <- "E"
pspectra.q.met$SamplingTime <- factor(pspectra.q.met$SamplingTime,
                                         levels = c("I", "E"))

#get column names for table
i = 1
ps_nr = i
df1 <- pspectra.q.met %>% filter(ps_nr == !!ps_nr)
cols <- df1 %>%
    rstatix::group_by(SampleGroup) %>%
    rstatix::t_test(norm_intensity ~ SamplingTime) %>%
    rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
    rstatix::add_significance(p.col = "p.adj") %>% 
    rstatix::add_xy_position(x = "SamplingTime", dodge = 0.8)
cols <- cbind(ps_nr = NA, cols)
ttests.df <- cols[!any(is.na(cols)),]
rm(cols)

#do t-tests and p-value adjustment
for(i in 1:length(pspectra)){
    if(i == 23){next}
    else{
        ps_nr = i
        df1 <- pspectra.q.met %>% filter(ps_nr == !!ps_nr)
        df1$intensity[is.na(df1$intensity)] <- 0
        df1$norm_intensity[is.na(df1$norm_intensity)] <- 0
        df2 <- df1 %>%
            rstatix::group_by(SampleGroup) %>%
            rstatix::t_test(norm_intensity ~ SamplingTime) %>%
            rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
            rstatix::add_significance(p.col = "p.adj") %>% 
            rstatix::add_xy_position(x = "SamplingTime", dodge = 0.8)
        df2 <- cbind(ps_nr = i, df2)
        ttests.df <- rbind(ttests.df, df2)
    }
}

ttests.df.sig <- ttests.df %>% filter(p.adj <= 0.05)
ttests.df.sig$ps_nr %>% unique() %>% length() #n = 58

#Wilcox tests between groups and ASW blank ----
#get column names for table
i = 1
ps_nr = i
df1 <- pspectra.q.met %>% filter(ps_nr == !!ps_nr)
df1 <- df1 %>% filter(SampleGroup3 != "ASW_metabolitestd")
df1$norm_intensity[is.na(df1$norm_intensity)] <- 0
df1 <- df1 %>% 
    rstatix::wilcox_test(norm_intensity ~ SampleGroup, ref.group = "ASW",
                         alternative = "greater") %>% 
    rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
    rstatix::add_significance(p.col = "p.adj") %>% 
    rstatix::add_xy_position()
cols <- cbind(ps_nr = NA, df1)
wilcox.df <- cols[!any(is.na(cols)),]
rm(cols)

#do t-tests and p-value adjustment
for(i in 1:length(pspectra)){
    if(i == 23){next} #ribitol
    ps_nr = i
    df1 <- pspectra.q.met %>% filter(ps_nr == !!ps_nr)
    df1 <- df1 %>% filter(SampleGroup3 != "ASW_metabolitestd")
    df1$norm_intensity[is.na(df1$norm_intensity)] <- 0
    k <- df1 %>% 
        dplyr::group_by(SampleGroup) %>% 
        dplyr::summarise(sum = sum(norm_intensity))
    if(nrow(k[k$sum==0,])>=2){next}
    df2 <- df1 %>% 
        rstatix::wilcox_test(norm_intensity ~ SampleGroup, ref.group = "ASW",
                             alternative = "greater") %>% 
        rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
        rstatix::add_significance(p.col = "p.adj") %>% 
        rstatix::add_xy_position()
    df2 <- cbind(ps_nr = i, df2)
    wilcox.df <- rbind(wilcox.df, df2)
}

wilcox.df.sig <- wilcox.df %>% filter(p.adj  <= 0.05)
wilcox.df.sig$ps_nr %>% unique() %>% length() #n = 77

w <- wilcox.df.sig$ps_nr %>% unique() %>% as.character()
t <- ttests.df.sig$ps_nr %>% unique() %>% as.character()

setdiff(t, w)
setdiff(w, t)
#both contain ps not in the other

#16: Get annotation for interesting pseudospectra----
#set coordinates for plots
ttests.df.sig$xmin[ttests.df.sig$SampleGroup == "ASW"] <- 0.8
ttests.df.sig$xmax[ttests.df.sig$SampleGroup == "ASW"] <- 1.2
ttests.df.sig$xmin[ttests.df.sig$SampleGroup == "PW_SEAGRASS"] <- 1.8
ttests.df.sig$xmax[ttests.df.sig$SampleGroup == "PW_SEAGRASS"] <- 2.2
ttests.df.sig$xmin[ttests.df.sig$SampleGroup == "WC_ALGAE"] <- 2.8
ttests.df.sig$xmax[ttests.df.sig$SampleGroup == "WC_ALGAE"] <- 3.2
ttests.df.sig$xmin[ttests.df.sig$SampleGroup == "WC_SEAGRASS"] <- 3.8
ttests.df.sig$xmax[ttests.df.sig$SampleGroup == "WC_SEAGRASS"] <- 4.2

wilcox.df.sig$xmin <- 0.8


#significant differences within algae incubations----
alg.t <- ttests.df.sig %>% 
    filter(SampleGroup == "WC_ALGAE")
a <- alg.t$ps_nr
#plot intensities
for (i in 1:length(a)){
    ps_nr = a[i]
    df1 <- pspectra.q.met %>% filter(ps_nr == !!ps_nr)
    df2 <- ttests.df.sig%>% filter(ps_nr == !!ps_nr)
    df3 <- wilcox.df.sig%>% filter(ps_nr == !!ps_nr)
    if (ps_nr %in% w & ps_nr %in% t){
        p <- ggplot(data = df1, aes(x=SampleGroup, y=sqrt(norm_intensity))) +
            geom_boxplot(aes(fill=SamplingTime), position = position_dodge(width = 0.8)) + 
            geom_point(pch=21, size=2,mapping =aes(fill=SamplingTime), 
                       position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("#FEC000","#2BB6AF"),
                              name = "Sampling time or sample type") +
            scale_y_continuous(name=expression(sqrt(Ribitol~normalised~intensity~(a.u.)))) +
            labs(x = "Sample group",
                 title = paste("Pseudospectra number ", ps_nr)) +
            add_pvalue(df2, xmin = "xmin", xmax = "xmax", 
                       y.position = sqrt(df2$y.position),
                       label = "p.adj.signif",
                       label.size = 5, fontfamily = "Arial") +
            add_pvalue(df3, xmin = "xmin", xmax = "xmax", 
                       y.position = sqrt(df3$y.position),
                       label = "p.adj.signif",
                       label.size = 5, fontfamily = "Arial",
                       step.increase = 0.1) +
            theme_bw() + 
            theme(text=element_text(size=12, family = "Arial", colour = "#262626"),
                  plot.title = element_text(size=14, family = "Arial", face = "bold",
                                            colour = "#262626", hjust = 0.5),
                  panel.grid = element_blank(),
                  legend.position = "top",
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "#262626"),
                  panel.spacing = unit(0.1, "lines"))
        png(filename = paste0("plots/algae/", ps_nr, ".png"),
            width = 18, height = 9, res = 300, units = "in")
        print(p)
        dev.off()
    }
    else if (ps_nr %in% w & !ps_nr %in% t){
        p <- ggplot(data = df1, aes(x=SampleGroup, y=sqrt(norm_intensity))) +
            geom_boxplot(aes(fill=SamplingTime), position = position_dodge(width = 0.8)) + 
            geom_point(pch=21, size=2,mapping =aes(fill=SamplingTime), 
                       position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("#FEC000","#2BB6AF"),
                              name = "Sampling time or sample type") +
            scale_y_continuous(name=expression(sqrt(Ribitol~normalised~intensity~(a.u.)))) +
            labs(x = "Sample group",
                 title = paste("Pseudospectra number ", ps_nr)) +
            add_pvalue(df3, xmin = "xmin", xmax = "xmax", 
                       y.position = sqrt(df3$y.position),
                       label = "p.adj.signif",
                       label.size = 5, fontfamily = "Arial",
                       step.increase = 0.1) +
            theme_bw() + 
            theme(text=element_text(size=12, family = "Arial", colour = "#262626"),
                  plot.title = element_text(size=14, family = "Arial", face = "bold",
                                            colour = "#262626", hjust = 0.5),
                  panel.grid = element_blank(),
                  legend.position = "top",
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "#262626"),
                  panel.spacing = unit(0.1, "lines"))
        png(filename = paste0("plots/algae/", ps_nr, ".png"),
            width = 18, height = 9, res = 300, units = "in")
        print(p)
        dev.off()
    }
    else if (!ps_nr %in% w & ps_nr %in% t){
        p <- ggplot(data = df1, aes(x=SampleGroup, y=sqrt(norm_intensity))) +
            geom_boxplot(aes(fill=SamplingTime), position = position_dodge(width = 0.8)) + 
            geom_point(pch=21, size=2,mapping =aes(fill=SamplingTime), 
                       position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("#FEC000","#2BB6AF"),
                              name = "Sampling time or sample type") +
            scale_y_continuous(name=expression(sqrt(Ribitol~normalised~intensity~(a.u.)))) +
            labs(x = "Sample group",
                 title = paste("Pseudospectra number ", ps_nr)) +
            add_pvalue(df2, xmin = "xmin", xmax = "xmax", 
                       y.position = sqrt(df2$y.position),
                       label = "p.adj.signif",
                       label.size = 5, fontfamily = "Arial") +
            theme_bw() + 
            theme(text=element_text(size=12, family = "Arial", colour = "#262626"),
                  plot.title = element_text(size=14, family = "Arial", face = "bold",
                                            colour = "#262626", hjust = 0.5),
                  panel.grid = element_blank(),
                  legend.position = "top",
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "#262626"),
                  panel.spacing = unit(0.1, "lines"))
        png(filename = paste0("plots/algae/", ps_nr, ".png"),
            width = 18, height = 9, res = 300, units = "in")
        print(p)
        dev.off()
    }
}
#annotations based on standard compound table
algae.mlist1 <- list()
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



#significant differences between algae incubations and ASW blanks (not used)----
alg.w <- wilcox.df.sig %>% 
    filter(group2 == "WC_ALGAE") #n=61
a <- alg.w$ps_nr
a <- a[!a %in% alg.t$ps_nr] #n=60
#plot
for (i in 1:length(a)){
    ps_nr = a[i]
    df1 <- pspectra.q.met %>% filter(ps_nr == !!ps_nr)
    df2 <- ttests.df.sig%>% filter(ps_nr == !!ps_nr)
    df3 <- wilcox.df.sig%>% filter(ps_nr == !!ps_nr)
    if (ps_nr %in% w & ps_nr %in% t){
        p <- ggplot(data = df1, aes(x=SampleGroup, y=sqrt(norm_intensity))) +
            geom_boxplot(aes(fill=SamplingTime), position = position_dodge(width = 0.8)) + 
            geom_point(pch=21, size=2,mapping =aes(fill=SamplingTime), 
                       position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("#FEC000","#2BB6AF"),
                              name = "Sampling time or sample type") +
            scale_y_continuous(name=expression(sqrt(Ribitol~normalised~intensity~(a.u.)))) +
            labs(x = "Sample group",
                 title = paste("Pseudospectra number ", ps_nr)) +
            add_pvalue(df2, xmin = "xmin", xmax = "xmax", 
                       y.position = sqrt(df2$y.position),
                       label = "p.adj.signif",
                       label.size = 5, fontfamily = "Arial") +
            add_pvalue(df3, xmin = "xmin", xmax = "xmax", 
                       y.position = sqrt(df3$y.position),
                       label = "p.adj.signif",
                       label.size = 5, fontfamily = "Arial",
                       step.increase = 0.1) +
            theme_bw() + 
            theme(text=element_text(size=12, family = "Arial", colour = "#262626"),
                  plot.title = element_text(size=14, family = "Arial", face = "bold",
                                            colour = "#262626", hjust = 0.5),
                  panel.grid = element_blank(),
                  legend.position = "top",
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "#262626"),
                  panel.spacing = unit(0.1, "lines"))
        png(filename = paste0("plots/algae/algae-v-ASW/", ps_nr, ".png"),
            width = 18, height = 9, res = 300, units = "in")
        print(p)
        dev.off()
    }
    else if (ps_nr %in% w & !ps_nr %in% t){
        p <- ggplot(data = df1, aes(x=SampleGroup, y=sqrt(norm_intensity))) +
            geom_boxplot(aes(fill=SamplingTime), position = position_dodge(width = 0.8)) + 
            geom_point(pch=21, size=2,mapping =aes(fill=SamplingTime), 
                       position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("#FEC000","#2BB6AF"),
                              name = "Sampling time or sample type") +
            scale_y_continuous(name=expression(sqrt(Ribitol~normalised~intensity~(a.u.)))) +
            labs(x = "Sample group",
                 title = paste("Pseudospectra number ", ps_nr)) +
            add_pvalue(df3, xmin = "xmin", xmax = "xmax", 
                       y.position = sqrt(df3$y.position),
                       label = "p.adj.signif",
                       label.size = 5, fontfamily = "Arial",
                       step.increase = 0.1) +
            theme_bw() + 
            theme(text=element_text(size=12, family = "Arial", colour = "#262626"),
                  plot.title = element_text(size=14, family = "Arial", face = "bold",
                                            colour = "#262626", hjust = 0.5),
                  panel.grid = element_blank(),
                  legend.position = "top",
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "#262626"),
                  panel.spacing = unit(0.1, "lines"))
        png(filename = paste0("plots/algae/algae-v-ASW/", ps_nr, ".png"),
            width = 18, height = 9, res = 300, units = "in")
        print(p)
        dev.off()
    }
}
#remove ps_nr = 34 -> higher in blank
a <- a[a != 34]

#annotations based on standard compound table
algae.mlist2 <- list()
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
    algae.mlist2[[i]] <- df1 %>% dplyr::group_by(compound, ion) %>% 
        dplyr::summarise(n=n())
}
names(algae.mlist2) <- as.character(a)
algae.mlist2 <- bind_rows(algae.mlist2, .id = "ps_nr")

#annotations by comparison to mona db
algae.mona2 <- data.frame(ps_nr = NA, rt = NA,
                          annot1 = NA, dotprod1 = NA, 
                          annot2 = NA,dotprod2 = NA,
                          annot3 = NA,dotprod3 = NA, 
                          annot4 = NA, dotprod4 = NA,
                          annot5 = NA, dotprod5 = NA)
for (i in 1:length(a)){
    ps_nr = a[i]
    l <- lapply(mona.spec,
                function(x) compareSpectra(x,
                                           spclist.by.sample.spec.norm[[3]][[ps_nr]],
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
    algae.mona2 <- rbind(algae.mona2, tmp.df)
    
}
algae.mona2 <- na.omit(algae.mona2)

#
#17: Compare light and dark incubations of algae----
pspectra.q.met.alg <- pspectra.q.met %>% filter(Plant == "ALGAE")
pspectra.q.met.alg$SampleGroup4 <- paste0(pspectra.q.met.alg$SampleGroup,
                                          "_",
                                          pspectra.q.met.alg$SamplingTime)
#compare between sampling times within the same condition----


#get column names for table
i = 1
ps_nr = i
df1 <- pspectra.q.met.alg %>% filter(ps_nr == !!ps_nr)
cols <- df1 %>%
    rstatix::group_by(SampleGroup2) %>%
    rstatix::t_test(norm_intensity ~ SamplingTime) %>%
    rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
    rstatix::add_significance(p.col = "p.adj") %>% 
    rstatix::add_xy_position(x = "SamplingTime", dodge = 0.8)
cols <- cbind(ps_nr = NA, cols)
ttests.df.light <- cols[!any(is.na(cols)),]
rm(cols)

#do t-tests and p-value adjustment
for(i in 1:length(pspectra)){
    if(i == 23){next}
    else{
        ps_nr = i
        df1 <- pspectra.q.met.alg %>% filter(ps_nr == !!ps_nr)
        df1$intensity[is.na(df1$intensity)] <- 0
        df1$norm_intensity[is.na(df1$norm_intensity)] <- 0
        df2 <- df1 %>%
            rstatix::group_by(SampleGroup2) %>%
            rstatix::t_test(norm_intensity ~ SamplingTime) %>%
            rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
            rstatix::add_significance(p.col = "p.adj") %>% 
            rstatix::add_xy_position(x = "SamplingTime", dodge = 0.8)
        df2 <- cbind(ps_nr = i, df2)
        ttests.df.light <- rbind(ttests.df.light, df2)
    }
}

ttests.df.light.sig <- ttests.df.light %>% filter(p.adj <= 0.05)
ttests.df.light.sig$ps_nr %>% unique() %>% length() #n = 2

#set coordinates for plotting
ttests.df.light.sig$xmin[ttests.df.light.sig$SampleGroup2 == "WC_ALGAE_light"] <- 1.8
ttests.df.light.sig$xmax[ttests.df.light.sig$SampleGroup2 == "WC_ALGAE_light"] <- 2.2

a <- ttests.df.light.sig$ps_nr

#plot intensities
for (i in 1:length(a)){
    ps_nr = a[i]
    df1 <- pspectra.q.met.alg %>% filter(ps_nr == !!ps_nr)
    df2 <- ttests.df.light.sig%>% filter(ps_nr == !!ps_nr)
    p <- ggplot(data = df1, aes(x=Condition, y=sqrt(norm_intensity))) +
        geom_boxplot(aes(fill=SamplingTime), position = position_dodge(width = 0.8)) + 
        geom_point(pch=21, size=2,mapping =aes(fill=SamplingTime), 
                   position = position_dodge(width = 0.8)) +
        scale_fill_manual(values = c("#FEC000","#2BB6AF"),
                          name = "Sampling time or sample type") +
        scale_y_continuous(name=expression(sqrt(Ribitol~normalised~intensity~(a.u.)))) +
        labs(x = "Light condition", title = paste("Pseudospectra number ", ps_nr)) +
        add_pvalue(df2, xmin = "xmin", xmax = "xmax", 
                   y.position = sqrt(df2$y.position),
                   label = "p.adj.signif",
                   label.size = 5, fontfamily = "Arial") +
        theme_bw() + 
        theme(text=element_text(size=12, family = "Arial", colour = "#262626"),
              plot.title = element_text(size=14, family = "Arial", face = "bold",
                                        colour = "#262626", hjust = 0.5),
              panel.grid = element_blank(),
              legend.position = "top",
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "#262626"),
              panel.spacing = unit(0.1, "lines"))
    
    png(filename = paste0("plots/algae/light-conditions/", ps_nr, ".png"),
        width = 9, height = 4.5, res = 300, units = "in")
    print(p)
    dev.off()
}
#annotations based on standard compound table
algae.mlist3 <- list()
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
    algae.mlist3[[i]] <- df1 %>% dplyr::group_by(compound, ion) %>% 
        dplyr::summarise(n=n())
}
names(algae.mlist3) <- as.character(a)
algae.mlist3 <- bind_rows(algae.mlist3, .id = "ps_nr")
#annotations by comparison to mona db
algae.mona3 <- data.frame(ps_nr = NA, rt = NA,
                          annot1 = NA, dotprod1 = NA, 
                          annot2 = NA,dotprod2 = NA,
                          annot3 = NA,dotprod3 = NA, 
                          annot4 = NA, dotprod4 = NA,
                          annot5 = NA, dotprod5 = NA)
for (i in 1:length(a)){
    ps_nr = a[i]
    l <- lapply(mona.spec,
                function(x) compareSpectra(x,
                                           spclist.by.sample.spec.norm[[1]][[ps_nr]],
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
    algae.mona3 <- rbind(algae.mona3, tmp.df)
    
}
algae.mona3 <- na.omit(algae.mona3)




