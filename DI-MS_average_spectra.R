#This script can be used to average spectra recorded during direct infusion
#mass spectrometry within a defined time range. Spectra (MS1) will be averaged
#within each sample - this could also be changed to be within groups etc.
#if desired.

#Written by Margot Bligh, February 2020

#1: Install packages --------------------------------------------------------

library(xcms)
library(MSnbase)
library(mzR)
library(tidyverse)
library(data.table)
library(hablar)
library(ggplot2)

#2. Import MS1 data --------------------------------------------------------

setwd("~/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/1_standards/1C_GCC-SPE-yield")

#get file paths to mzML files
raw_files_path <- dir(path = "./mzML-files", all.files = FALSE, full.names = TRUE)
#only keep MS1 spectra
raw_files_path <- raw_files_path[-grep("MSMS", raw_files_path)]

#create phenodata data.frame
pd <- data.frame(name = sub(basename(raw_files_path), 
                                   pattern = ".mzML", 
                                   replacement = ""), 
                 concentration = sub(".mzML", "", basename(raw_files_path)) %>% 
                     sub("20201030_ESI_neg_MS31_GCCSPE-yield-test_", "", .) %>% 
                     sub("-.*", "", .) %>% 
                     sub("B\\d", 0, .) %>% 
                     sub("2$", 2.5, .),
                 fraction = sub(".mzML", "", basename(raw_files_path)) %>% 
                     sub("20201030_ESI_neg_MS31_GCCSPE-yield-test_", "", .) %>% 
                     gsub("\\d+-", "", .),
                 stringsAsFactors = FALSE)
#load the raw data
raw_data <- readMSData(files = raw_files_path, 
                       pdata = new("NAnnotatedDataFrame", 
                                   pd), 
                       mode = "onDisk")

#2. Get BPC (optional)------------------------------

bpis_all <- chromatogram(raw_data, aggregationFun = "max")
names(bpis_all) <- pd$name

#to plot all
#plot(bpis_all)

#to plot for a specific file, e.g.
#plot(bpis_all[[2]])

#3: Plot TIC (optional - airbubble was observed during experiment for a sample)-----

#extract tic
tic <- chromatogram(raw_data)

#format sample name
pd$fraction_fmt <- pd$fraction %>% 
    sub("E", "eluate ", .) %>% 
    sub("I", "initial", .)

#initialise empty data frame
tic.df <- data.frame(fraction = as.character,
                     intensity = as.numeric,
                     rt = as.numeric)

#fill in data frame
for (i in 1:length(pd$name)) {
    intensity = tic[[i]]@intensity
    rt = tic[[i]]@rtime
    fraction = rep(pd$fraction_fmt[i], length(rt))
    temp <- data.frame(fraction = fraction,
                       intensity = intensity,
                       rt = rt)
    tic.df <- rbind(tic.df,
                    temp)
}

#set variables to factors
tic.df$fraction <- factor(tic.df$fraction,
                          levels = unique(tic.df$fraction))
tic.df$fraction <- relevel(tic.df$fraction,
                           "initial")

#plot
svg("./analysis/tic_1ng.svg",
     height = 3,
     width = 12)

ggplot() +
    geom_line(mapping = aes(rt,
                            intensity),
              data = tic.df,
              colour = "black",
              lwd = 1.2) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          axis.text = element_text(size = 12)) +
    labs(x= "Infusion time (sec)",
         y = "Intensity (a.u.)") +
    scale_x_continuous(breaks = seq(0, 130, 10),
                       limits = c(0,130),
                       expand = c(0,0)) +
    facet_grid(rows = vars(fraction)) +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x = element_blank())
dev.off()


#4. Filter spectra by retention time  ------------------------------
#bubble in capillary from 1 min to 1 min 30 s for 1-E1
#filter time to 30 s to 60 s
#all others filter from 60 s to 90 s

#check number of spectra per file
table(fromFile(raw_data))

#filter
raw_data.rt30to60 <- raw_data %>%
    filterFile(1) %>% 
    filterRt(rt= c(30, 60))
raw_data.rt60to90 <- raw_data %>%
    filterRt(rt= c(60, 90)) %>% 
    filterFile(2:69)
raw_data.rtfilter <- c(raw_data.rt30to60, raw_data.rt60to90)

#check number of spectra per file has decreased and is approx equal across samples
lapply(raw_data.rtfilter, function(x) {table(fromFile(x))})

#5. Average all spectra within each sample -----------------------------------------
#average spectra by file
#this may give an error telling you something about vector is too large...
#fixed by just clearing environment and restarting, or restarting R...

combined <- lapply(raw_data.rtfilter, 
                   function(x) {
                       combineSpectra(x, 
                                      fcol = "fileIdx", 
                                      mzFun = base::mean)})


#each file should contain 1 spectrum
lapply(combined, function(x) {table(fromFile(x))})

#6. Plot examples to check -----------------------------------------
#EXAMPLE 1: sample ="20201030_ESI_neg_MS31_GCCSPE-yield-test_500-I"

x <- raw_data.rtfilter[[2]] %>% filterFile(57)
#check name of file
x$name 

#plot
par(mfrow = c(2, 1), mar = c(4.3, 4, 1, 1))
plot(mz(x[[20]]), intensity(x[[20]]), type = "h", col = "light blue")
points(mz(x[[1]]), intensity(x[[1]]), type = "h", col = "darkolivegreen4")
points(mz(x[[10]]), intensity(x[[10]]), type = "h", col = "thistle")
plot(mz(combined[[2]][[57]]), intensity(combined[[2]][[57]]), type = "h", col = "black")

#EXAMPLE 1: sample = "20201030_ESI_neg_MS31_GCCSPE-yield-test_2-5-E1"
x <- raw_data.rtfilter[[2]] %>% filterFile(18)

par(mfrow = c(2, 1), mar = c(4.3, 4, 1, 1))
plot(mz(x[[20]]), intensity(x[[20]]), type = "h", col = "light blue")
points(mz(x[[1]]), intensity(x[[1]]), type = "h", col = "darkolivegreen4")
points(mz(x[[10]]), intensity(x[[10]]), type = "h", col = "thistle")

plot(mz(combined[[2]][[18]]), intensity(combined[[2]][[18]]), type = "h", col = "black")

#7. Write averaged data to file -----------------------------------------
#create output directory if it doesn't exist
dir.create(file.path("mzML-files-averaged"), showWarnings = TRUE)
#create new file paths for averaged files
averaged_files_path <- raw_files_path %>% 
    sub("mzML-files", "mzML-files-averaged", .) %>%
    sub("\\.mzML", "-averaged.mzML", .)

#Feature data of OnDiskMSnExp contains a column named "fileIdx", feature data of MSnExp does not.
#Before averaging, object is class "OnDiskMSnExp"
#After averaging, onject is class "MSnExp"
#Need to ensure that column "fileIdx" is not copied during the combine step
#You will get the following error if an MSnExp object contains a fileIdx column in feature data:
#Error in validObject(x) : invalid class "MSnExp" object: Mismatch of files in assayData and processingData.
fData(combined[[1]])$fileIdx <- NULL 
fData(combined[[2]])$fileIdx <- NULL 

#write data to file
writeMSData(combined[[1]], 
            averaged_files_path[1])

writeMSData(combined[[2]], 
            averaged_files_path[2:69])



  