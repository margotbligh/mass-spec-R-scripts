# Analysis of SeaMet samples from Finland

#1 packages########################
### loading required packages ###

rm(list=ls())

#mass spectrometry (have to be installed using Bioconductor)
library(MSnbase)
library(msdata)
library(xcms)
library(CAMERA)

## data handling (can be installed using the Packages interface)
library(magrittr)
library(tidyverse)
library(reshape2)
library(stringi)
library(stringr)

## plotting and having fun with color (can be installed using the Packages interface)
library(ggplot2)
library(pcaMethods)
library(png)
library("viridis") 
library("paletteer")
library(scico)
library(pander)

source("https://raw.githubusercontent.com/jorainer/xcms-gnps-tools/master/customFunctions.R")

#2 working dir ########################
### setting up working directories ###

general_dir<-"C:/Users/admin/ownCloud/Assemble/Data/SeaMet/GCMS data analysis/"
data_dir<-"C:/Users/admin/ownCloud/Assemble/Data/SeaMet/GCMS data analysis/mzML data/"
results_dir<-"C:/Users/admin/ownCloud/Assemble/Data/SeaMet/GCMS data analysis/results/"
script_dir<-"C:/Users/admin/ownCloud/Assemble/Data/SeaMet/GCMS data analysis/scripts/"


#3 data import ####################

#3.1 metadata ######
### loading data and metadata ###
# Create a vector with the files
mzMLfiles = list.files(data_dir, recursive = TRUE, full.names = TRUE, pattern = ".mzML")

### create metadata ###
#vector with sample type
temp1_samplename<-unlist(strsplit(mzMLfiles,"/"))
temp2_samplename<-temp1_samplename[seq(10,760,10)]

#create empty data frame
raw_metadata<-data.frame(plant = character(), 
                         waterbody = character(),
                         samplingtime = character(),
                         condition = character(),
                         replicate = integer(),
                         seametdate = character(),
                         runinsequence = integer(),
                         enriched13C = character(),
                         stringsAsFactors = FALSE)

#create empty temporary vectors
temp_plant<-NULL
temp_waterbody<-NULL
temp_samplingtime<-NULL
temp_condition<-NULL
temp_replicate<-NULL
temp_seametdate<-NULL
temp_runinsequence<-NULL
temp_enriched13C<-NULL

#run for loop to extract information from sample names and fill into metadata data frame
for (i in seq(1,length(temp2_samplename),1)){
  if (grepl("BLK",temp2_samplename[i])){
    BLK_temp<-unlist(strsplit(temp2_samplename[i],"_"))
    temp_plant<-"NA"
    temp_waterbody<-"MQ"
    temp_samplingtime<-"NA"
    temp_condition<-"NA"
    temp_replicate<-"NA"
    temp_seametdate<-BLK_temp[2]
    RUN_temp<-unlist(strsplit(BLK_temp[3]," "))
    temp_runinsequence<-RUN_temp[1]
    temp_enriched13C<-"NO"
  }
  if (grepl("ASW",temp2_samplename[i])){
    ASW_temp<-unlist(strsplit(temp2_samplename[i],"_"))
    temp_plant<-"NA"
    temp_waterbody<-"ASW"
    temp_samplingtime<-"NA"
    temp_condition<-ASW_temp[2]
    temp_replicate<-"NA"
    temp_seametdate<-ASW_temp[3]
    RUN_temp<-unlist(strsplit(ASW_temp[4]," "))
    temp_runinsequence<-RUN_temp[1]
    temp_enriched13C<-"NO"
  }
  if (grepl("A_",temp2_samplename[i])){
    ALGAE_temp<-unlist(strsplit(temp2_samplename[i],"_"))
    temp_plant<-"ALGAE"
    EXP_temp<-gsub("[^a-zA-Z]", "_", ALGAE_temp[2])
    EXP_temp2<-unlist(strsplit(EXP_temp, "_"))
    temp_waterbody<-EXP_temp2[1]
    temp_samplingtime<-EXP_temp2[2]
    REP_COND_temp<-str_extract(ALGAE_temp[2],"[0-9]+")
    if (REP_COND_temp <= 3){
      temp_condition<-"light"
      temp_replicate<-REP_COND_temp
    }
    if (REP_COND_temp >=4){
      temp_condition<-"dark"
      temp_replicate<-as.numeric(REP_COND_temp)-3
    }
    temp_seametdate<-ALGAE_temp[3]
    RUN_temp<-unlist(strsplit(ALGAE_temp[4]," "))
    temp_runinsequence<-RUN_temp[1]
    temp_enriched13C<-"NO"
  }
  if (grepl("S_",temp2_samplename[i])){
    SEAGRASS_temp<-unlist(strsplit(temp2_samplename[i],"_"))
    temp_plant<-"SEAGRASS"
    EXP_temp<-gsub("[^a-zA-Z]", "_", SEAGRASS_temp[2])
    EXP_temp2<-unlist(strsplit(EXP_temp, "_"))
    temp_waterbody<-EXP_temp2[1]
    temp_samplingtime<-EXP_temp2[length(EXP_temp2)]
    REP_COND_temp<-str_extract(SEAGRASS_temp[2],"[0-9]+")
    if (as.numeric(REP_COND_temp) <= 3){
      temp_condition<-"light"
      temp_replicate<-REP_COND_temp
      temp_enriched13C<-"NO"
    }
    if (as.numeric(REP_COND_temp) >=4 && as.numeric(REP_COND_temp) <=6){
      temp_condition<-"dark"
      temp_replicate<-as.numeric(REP_COND_temp)-3
      temp_enriched13C<-"NO"
    }
    if (as.numeric(REP_COND_temp) >=7 && as.numeric(REP_COND_temp) <=9){
      temp_condition<-"light"
      temp_replicate<-as.numeric(REP_COND_temp)-6
      temp_enriched13C<-"YES"
    }
    if (as.numeric(REP_COND_temp) >=10){
      temp_condition<-"dark"
      temp_replicate<-as.numeric(REP_COND_temp)-9
      temp_enriched13C<-"YES"
    }
    temp_seametdate<-SEAGRASS_temp[3]
    RUN_temp<-unlist(strsplit(SEAGRASS_temp[4]," "))
    temp_runinsequence<-RUN_temp[1]
  }
  raw_metadata<-rbind(raw_metadata, c(temp_plant, 
                                      temp_waterbody, 
                                      temp_samplingtime,
                                      temp_condition,
                                      temp_replicate,
                                      temp_seametdate,
                                      temp_runinsequence,
                                      temp_enriched13C))
  colnames(raw_metadata)<-c("Plant", "Waterbody", "SamplingTime", "Condition", 
                            "Replicate", "SeametDate", "RunInSequence", "Enrichment13C")
  temp_plant<-NULL
  temp_waterbody<-NULL
  temp_samplingtime<-NULL
  temp_condition<-NULL
  temp_replicate<-NULL
  temp_seametdate<-NULL
  temp_runinsequence<-NULL
  temp_enriched13C<-NULL
}

raw_metadata$SampleGroup<-paste(raw_metadata$Waterbody,raw_metadata$Condition, 
                                raw_metadata$Plant, sep="_")


#3.2 load data ####
### import data with information ###

files<-mzMLfiles #input files

pd <- data.frame(sample_name = sub(basename(files), pattern = ".mzML",
                                   replacement = "", fixed = TRUE),
                 sample_group = raw_metadata$SampleGroup,
                 waterbody = raw_metadata$Waterbody,
                 plant = raw_metadata$Plant,
                 sampling_time= raw_metadata$SamplingTime,
                 condition = raw_metadata$Condition,
                 seamet_date = raw_metadata$SeametDate,
                 run_nr = raw_metadata$RunInSequence,
                 enriched = raw_metadata$Enrichment13C,
                 stringsAsFactors = FALSE)

#Load the raw data as an OnDiskMSnExp using MSnbase package
raw_data <- readMSData(files = files, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk", msLevel. = 1, centroided. = TRUE) 


#3.3 define colors #####
#Define a color palette for the variable "waterbody", which is the sample_group in the data set
group_color <- c(scico(5,end=0.9, palette = "lapaz", alpha = 0.5), 
                 scico(4, end=0.9, palette = "lajolla", alpha = 0.5))

#group_color[2]<-scico(1,palette="roma", alpha=0.2)
#group_color <- wesanderson::wes_palette("Darjeeling1", n=3)
names(group_color) <- c("WC_light_ALGAE", "WC_dark_ALGAE", "ASW_blank_NA", 
                        "ASW_metabolitestd_NA", "MQ_NA_NA",
                        "PW_dark_SEAGRASS", "PW_light_SEAGRASS",
                        "WC_dark_SEAGRASS", "WC_light_SEAGRASS")
leg <- c("WC_light_ALGAE", "WC_dark_ALGAE", "ASW_blank_NA", 
         "ASW_metabolitestd_NA", "MQ_NA_NA",
         "PW_dark_SEAGRASS", "PW_light_SEAGRASS",
         "WC_dark_SEAGRASS", "WC_light_SEAGRASS")
sample_colors<-c(rep(group_color[1],6),
                 rep(group_color[2],6),
                 rep(group_color[3],3),
                 rep(group_color[4],3),
                 rep(group_color[5],12),
                 rep(group_color[6],6),
                 rep(group_color[7],5),
                 rep(group_color[6],6),
                 rep(group_color[7],5),
                 rep(group_color[8],6),
                 rep(group_color[9],6),
                 rep(group_color[8],6),
                 rep(group_color[9],5))


#4 overview BPC and TIC#####
### assess the data with BPC and TIC ###

## Plotting BPC (base peak chromatogram)
BPC <- raw_data %>% chromatogram(.,rt=c(400,2000), aggregationFun = "max")

par(mar=c(5,4,3,15))
plot(BPC,col = group_color[raw_data$sample_group], main = "Base peak Chromatograms \nLB")
legend ("right", inset=c(0,0), xpd = TRUE, legend=leg, fill = group_color, cex = 0.6)


# Plot TIC per file 
#extract TIC and filenames
tic <- tic(raw_data, initial = FALSE) #Get tic
from <- fromFile(raw_data) #get sample ID
mydata <-data.frame(tic, from) #combine both vector in a dataframe
pd <- pd %>%
  rownames_to_column(var = "rowname")
mydata <- merge(mydata, pd, by.x = "from", by.y = "rowname")

#plot TIC by analysis day on GC-MS
ggplot(mydata, aes(x=as.integer(seamet_date), y=log2(tic), group=sample_name, fill=sample_group)) +
  #scale_fill_scico_d(4, palette="romaO", name= "Sample type", begin = 0, end = 0.8) +
  geom_boxplot() +
  coord_cartesian(ylim = c(10, 3e+01)) +
  labs(title = "TIC per file \nFinland samples", x="", y = "Intensity")+
  theme(legend.position = "right")

################################### ATTENTION, file removed ###############################
### based on very low TIC, S_PW9E_20201120_016 taken out of the analysis! ###


#5 peak picking ####

#5.1 test parameters ####
### testing parameters for peak picking ###

#parameters, empirical with ribitol peak at 217
mycentwaveparams<-CentWaveParam()
mycentwaveparams@peakwidth<-c(0.6,6)
mycentwaveparams@ppm<-300
mycentwaveparams@integrate<-1L
mycentwaveparams@noise<-0
mycentwaveparams@prefilter<-c(3,500)
mycentwaveparams@fitgauss<-T
mycentwaveparams@snthresh=6
mycentwaveparams@mzdiff=0.01
mycentwaveparams@mzCenterFun=c("apex")

#GNPS parameters
centW<-CentWaveParam()
centW@ppm = 100 #maximum m/z difference tolerated in ppm, 100 ppm = 0.1 difference
centW@peakwidth = c(5,10) # minimum, maximum peak width in sec
centW@snthresh = 6 # signal to noise ratio cut off 
centW@prefilter = c(3, 100) # 3 peaks with intensity equal or higher than 100 
centW@mzCenterFun = c("apex") # wMean on XMCS online
centW@integrate = 1L #peak limits are found through descenton the mexican hat filtered data
centW@mzdiff = 0.01 #minimum difference in m/z for peaks with overlapping retention times, can be negative to allow overlap
centW@fitgauss = TRUE
centW@noise = 0 #minimum intensity required for centroids to beconsidered in the first analysis step 
centW@verboseColumns = TRUE
#centW@sleep = 1 #number of seconds to pause between plotting peak finding cycles

# parameters that worked in terms of number of features for algae
algaeparams<-CentWaveParam()
algaeparams@peakwidth<-c(2,10)
algaeparams@ppm<-300


#5.2 ROI plotting ####

#ribitol test
rt_ribitol<-c(1000,1100)
mz_ribitol<-c(216, 219)

chr_raw<-chromatogram(raw_data, mz=mz_ribitol, rt=rt_ribitol)
plot(chr_raw, col = group_color[chr_raw$sample_group])
test_chrpeaks<-findChromPeaks(chr_raw, param = mycentwaveparams)
sample_colors <- group_color[test_chrpeaks$sample_group]
plot(test_chrpeaks, col = sample_colors,
     peakBg = sample_colors[chromPeaks(test_chrpeaks)[, "column"]])

#sugar test. Found Glucosamine at RT 1140
rt_sugars<-c(1050,1200)
mz_sugars_217<-c(217,220)
par(mfrow=c(2,3))
chr_raw2_initial_217<-chromatogram(raw_data, mz=mz_sugars_217, rt=rt_sugars)[1,c(2,4,6,8,10,12)]
plot(chr_raw2_initial_217, col= c(rep("red",3),rep("blue",3)), lwd=2, main="Algae initial 217-220")
legend("topright", xpd = TRUE, legend=leg, fill = group_color, inset=c(0,0))

chr_raw2_end_217<-chromatogram(raw_data, mz=mz_sugars_217, rt=rt_sugars)[1,c(1,3,5,9,11)]
plot(chr_raw2_end_217, col= c(rep("red",3),rep("blue",2)), lwd=2, main="Algae end 217-220")


mz_sugars_307<-c(307,310)
chr_raw2_initial_307<-chromatogram(raw_data, mz=mz_sugars_307, rt=rt_sugars)[1,c(2,4,6,8,10,12)]
plot(chr_raw2_initial_307, col= c(rep("red",3),rep("blue",3)), lwd=2, main="Algae initial 307-310")
#legend("topright", xpd = TRUE, legend=leg, fill = group_color, inset=c(0,0))

chr_raw2_end_307<-chromatogram(raw_data, mz=mz_sugars_307, rt=rt_sugars)[1,c(1,3,5,9,11)]
plot(chr_raw2_end_307, col= c(rep("red",3),rep("blue",2)), lwd=2, main="Algae end 307-310")

mz_sugars_319<-c(319,322)
chr_raw2_initial_319<-chromatogram(raw_data, mz=mz_sugars_319, rt=rt_sugars)[1,c(2,4,6,8,10,12)]
plot(chr_raw2_initial_319, col= c(rep("red",3),rep("blue",3)), lwd=2, main="Algae initial 319-322")
#legend("topright", xpd = TRUE, legend=leg, fill = group_color, inset=c(0,0))

chr_raw2_end_319<-chromatogram(raw_data, mz=mz_sugars_319, rt=rt_sugars)[1,c(1,3,5,9,11)]
plot(chr_raw2_end_319, col= c(rep("red",3),rep("blue",2)), lwd=2, main="Algae end 319-322")
#legend("topright", xpd = TRUE, legend=leg, fill = group_color, inset=c(0,0))


#sugar test. Moni
files
rt_sugars<-c(1400,1600)
mz_sugars_361<-c(361,363)
par(mfrow=c(3,2))
chr_raw2_initial_361<-chromatogram(raw_data, mz=mz_sugars_361, rt=rt_sugars)[37,c(38,39,40,41,42,43)]
plot(chr_raw2_initial_361, col= c(rep("red",3),rep("blue",3)), lwd=2, main="SG PW initial 361-363")
#legend("topright", xpd = TRUE, legend=leg, fill = group_color, inset=c(0,0))

chr_raw2_end_361<-chromatogram(raw_data, mz=mz_sugars_361, rt=rt_sugars)[37,c(38,39,40,41,42,43)]
plot(chr_raw2_end_361, col= c(rep("red",3),rep("blue",2)), lwd=2, main="Algae end 361-309")


mz_sugars_147<-c(147,149)
chr_raw2_initial_147<-chromatogram(raw_data, mz=mz_sugars_147, rt=rt_sugars)[37,c(38,39,40,41,42,43)]
plot(chr_raw2_initial_147, col= c(rep("red",3),rep("blue",3)), lwd=2, main="Algae initial 147-149")
#legend("topright", xpd = TRUE, legend=leg, fill = group_color, inset=c(0,0))

chr_raw2_end_147<-chromatogram(raw_data, mz=mz_sugars_147, rt=rt_sugars)[37,c(38,39,40,41,42,43)]
plot(chr_raw2_end_147, col= c(rep("red",3),rep("blue",2)), lwd=2, main="Algae end 147-149")

mz_sugars_217<-c(217,279)
chr_raw2_initial_217<-chromatogram(raw_data, mz=mz_sugars_217, rt=rt_sugars)[37,c(38,39,40,41,42,43)]
plot(chr_raw2_initial_217, col= c(rep("red",3),rep("blue",3)), lwd=2, main="Algae initial 217-219")
#legend("topright", xpd = TRUE, legend=leg, fill = group_color, inset=c(0,0))

chr_raw2_end_217<-chromatogram(raw_data, mz=mz_sugars_217, rt=rt_sugars)[37,c(38,39,40,41,42,43)]
plot(chr_raw2_end_217, col= c(rep("red",3),rep("blue",2)), lwd=2, main="Algae end 217-219")
#legend("topright", xpd = TRUE, legend=leg, fill = group_color, inset=c(0,0))


# Exploring an unknown

rt_unknown<-c(1100,1150)
mz_unknown<-c(203,204)
par(mfrow=c(1,2))
chr_unknown_203_initial<-chromatogram(raw_data, mz=mz_unknown, rt=rt_unknown)[1,c(59,61,63,65,67,69)]
plot(chr_unknown_203_initial, col= c(rep("red",3),rep("blue",3)), lwd=2, main="SG WC initial 203-204")
#legend("topright", xpd = TRUE, legend=leg, fill = group_color, inset=c(0,0))

chr_unknown_203_end<-chromatogram(raw_data, mz=mz_unknown, rt=rt_unknown)[1,c(60,62,64,66,68,70)]
plot(chr_unknown_203_end, col= c(rep("red",3),rep("blue",3)), lwd=2, main="SG WC end 203-204")


test_chrpeaks2<-findChromPeaks(chr_raw2, param = algaeparams)
plot(test_chrpeaks2, col=sample_colors,
     peakBg= sample_colors[chromPeaks(test_chrpeaks2)[,"column"]])

raw_data %>%
  filterMz(mz=mz_sugars) %>%
  filterRt(rt=rt_sugars) %>%
  filterFile(8)%>%
  plot(., type="XIC")

#sugar test 2
rt_sugars2<-c(900,1700)
mz_sugars2<-c(318,322)

chr_raw3<-chromatogram(raw_data, mz=mz_sugars2, rt=rt_sugars2)
plot(chr_raw3, col= group_color[chr_raw3$sample_group])

test_chrpeaks3<-findChromPeaks(chr_raw3, param = algaeparams)
plot(test_chrpeaks3, col=sample_colors,
     peakBg= sample_colors[chromPeaks(test_chrpeaks3)[,"column"]])


## Select one of your peak of  interest (in this case trehalose) in one file
tre_chr <- chromatogram(raw_data, rt = c(1500, 1600), mz = c(361.0, 361.4), aggregationFun = "max")[1,18]
## Plot the data
par(mfrow = c(1, 1), mar = c(4, 4.5, 1, 1))
plot(tre_chr)

## test run 
tre_peak <- findChromPeaks(tre_chr, param = cwp)

#5.3 ppm ####
### calculate and adapt ppm for dataset
## Restrict the data to signal from trehalose
tre <- raw_data %>%
  filterRt(rt = c(1530, 1565)) %>%
  filterMz(mz = c(360.0, 361.6)) %>%
  filterFile(c(13:18))
## Plot the data
plot(tre, type = "XIC") 
## Extract the trehalose data for one file as a data.frame
tre_df <- as(filterFile(tre, 1), "data.frame")
head(tre_df)
## calculate the difference in m/z values from consecutive scans in ppm
diff(tre_df$mz) * 1e6 / mean(tre_df$mz)

########################### NOTE peak shape ############################
### some peaks, in particular very high peaks, exhibit peak splitting 
### No parameters found that would prevent this. 
### Also merging neighboring peaks doesn't help, as smaller second peaks 
### are not picked at all. 

#mpp<-MergeNeighboringPeaksParam(expandRt = 4, ppm=300)
#test_chrpeaks3_merge<-refineChromPeaks(test_chrpeaks3, mpp)
#plot(test_chrpeaks3_merge, col=sample_colors,
#     peakBg= sample_colors[chromPeaks(test_chrpeaks3_merge)[,"column"]])


#5.4 peak picking #########################################
### peak picking based on empirically tested parameters ###


raw_peaks<-findChromPeaks(raw_data, param = centW)

save(raw_peaks, 
     file="C:/Users/admin/ownCloud/Assemble/Data/SeaMet/GCMS data analysis/results/raw_peaks_28Jan2021_centW.RData")
load("C:/Users/admin/ownCloud/Assemble/Data/SeaMet/GCMS data analysis/results/raw_peaks_28Jan2021_algaeparams.RData")


#5.5 assess peak picking ####

### ROIs of interest
test_pp<-chromatogram(raw_peaks, mz=mz_ribitol, rt=rt_ribitol)
plot(test_pp, col = group_color[test_pp$sample_group])
test_pp2<-chromatogram(raw_peaks, mz=mz_sugars, rt=rt_sugars)
plot(test_pp2, col = group_color[test_pp2$sample_group])
test_pp3<-chromatogram(raw_peaks, mz=mz_sugars2, rt=rt_sugars2)
plot(test_pp3, col = group_color[test_pp3$sample_group])


### peak distribution across RT
par(mfrow= c(1,1),mar = c(4, 20, 1, 1))
plotChromPeakImage(raw_peaks, binSize = 5, log = T)


### TIC
tic_peaks<-tic(raw_peaks, initial = FALSE)
from <- fromFile(raw_peaks) #get sample ID
mydata <-data.frame(tic_peaks, from) #combine both vector in a dataframe
pd <- pd %>%
  rownames_to_column(var = "rowname")
mydata <- merge(mydata, pd, by.x = "from", by.y = "rowname")

ggplot(mydata, aes(x=as.integer(seamet_date), y=log2(tic_peaks), group=sample_name, fill=sample_group)) +
  #scale_fill_scico_d(4, palette="romaO", name= "Sample type", begin = 0, end = 0.8) +
  geom_boxplot() +
  coord_cartesian(ylim = c(1e+1, 3e+01)) +
  labs(title = "TIC of picked peaks by run \nFinland SeaMet samples", x="", y = "Intensity")+
  theme(legend.position = "right")
# for mycentwaveparams, looks like sugars have been picked quite thoroughly,
# only mz range 318-322 looks like it has some un-picked peaks (in the std)

# for centW, very similar picture as with mycentwaveparams

# for algaeparams, the peak picking might actually look nicest -winner!

## Select one of your peak of  interest (in this case trehalose) in one file
algae_monosaccharide_chr <- chromatogram(raw_data, rt = c(1050, 1150), mz = c(319.0, 319.4),
                        aggregationFun = "max")[1,c(1,2)]
## Plot the data
par(mfrow = c(1, 1), mar = c(4, 4.5, 1, 1))
plot(algae_monosaccharide_chr, col = group_color[algae_monosaccharide_chr$sample_group])

legend ("topright", inset=c(0,0), xpd = TRUE, legend=leg, 
        fill = group_color[chromPeaks(algae_monosaccharide_chr)[,"column"]], cex = 0.8)


### summary table peaks

summary_fun <- function(z)
  c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))
T <- lapply(split.data.frame(
  chromPeaks(raw_peaks), f = chromPeaks(raw_peaks)[, "sample"]),
  FUN = summary_fun)
T <- do.call(rbind, T)
rownames(T) <- basename(fileNames(raw_peaks))
pandoc.table(
  T,
  caption = paste0("Summary statistics on identified chromatographic",
                   " peaks. Shown are number of identified peaks per",
                   " sample and widths/duration of chromatographic ",
                   "peaks."))


## Extract a list of per-sample peak intensities (in log2 scale)
ints <- split(log2(chromPeaks(raw_peaks)[, "into"]),
              f = chromPeaks(raw_peaks)[, "sample"])
boxplot(ints, varwidth = TRUE, col = group_color[raw_peaks$sample_group],
        ylab = expression(log[2]~intensity), main = "Peak intensities")
grid(nx = NA, ny = NULL)



#6 known compounds ####

#6.1 internal standards ####
#checking internal standards in all samples except BLK

# ribitol
mz_ribitol<-c(217,217.4)
rt_ribitol<-c(1000,1020)

chr_ribitol <- chromatogram(raw_peaks, mz = mz_ribitol, rt = rt_ribitol)
chromPeaks(chr_ribitol)

sample_colors <- group_color[chr_ribitol$sample_group]
plot(chr_ribitol, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_ribitol)[, "sample"]],
     peakBg = NA)
legend ("topright", inset=c(0,0), xpd = TRUE, legend=leg, 
        fill = group_color[chromPeaks(chr_ribitol)[,"column"]], cex = 0.8)

# cholestane

mz_cholestane<-c(372,373)
rt_cholestane<-c(1610,1630)

chr_cholestane <- chromatogram(raw_peaks, mz = mz_cholestane, rt = rt_cholestane)
chromPeaks(chr_cholestane)

sample_colors <- group_color[chr_cholestane$sample_group]
plot(chr_cholestane, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_cholestane)[, "sample"]],
     peakBg = NA)


#6.2 metabolite mix ####
# checking compounds in the metabolite mix

# Alanine

mz_alanine<-c(218,219)
rt_alanine<-c(930,950)

chr_alanine <- chromatogram(raw_peaks, mz = mz_alanine, rt = rt_alanine)
chromPeaks(chr_alanine)

sample_colors <- group_color[chr_alanine$sample_group]
plot(chr_alanine, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_alanine)[, "sample"]],
     peakBg = NA)


# Fructose

mz_fructose<-c(306,309)
rt_fructose<-c(1080,1150)

chr_fructose <- chromatogram(raw_peaks, mz = mz_fructose, rt = rt_fructose)
chromPeaks(chr_fructose)

sample_colors <- group_color[chr_fructose$sample_group]
plot(chr_fructose, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_fructose)[, "sample"]],
     peakBg = NA)


# Fumarate

mz_fumarate<-c(245,246)
rt_fumarate<-c(720,750)

chr_fumarate <- chromatogram(raw_peaks, mz = mz_fumarate, rt = rt_fumarate)
chromPeaks(chr_fumarate)

sample_colors <- group_color[chr_fumarate$sample_group]
plot(chr_fumarate, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_fumarate)[, "sample"]],
     peakBg = NA)


# Galactose

mz_galactose<-c(319,320)
rt_galactose<-c(1220,1240)

chr_galactose <- chromatogram(raw_peaks, mz = mz_galactose, rt = rt_galactose)
chromPeaks(chr_galactose)

sample_colors <- group_color[chr_galactose$sample_group]
plot(chr_galactose, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_galactose)[, "sample"]],
     peakBg = NA)


# Glucose

mz_glucose<-c(321,322)
rt_glucose<-c(1100,1150)

chr_glucose <- chromatogram(raw_peaks, mz = mz_glucose, rt = rt_glucose)
chromPeaks(chr_glucose)

sample_colors <- group_color[chr_glucose$sample_group]
plot(chr_glucose, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_glucose)[, "sample"]],
     peakBg = NA)


# Glutamate

mz_glutamate<-c(177,178)
rt_glutamate<-c(920,960)

chr_glutamate <- chromatogram(raw_peaks, mz = mz_glutamate, rt = rt_glutamate)
chromPeaks(chr_glutamate)

sample_colors <- group_color[chr_glutamate$sample_group]
plot(chr_glutamate, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_glutamate)[, "sample"]],
     peakBg = NA)

# glutamate peak not picked in sample 17, ASW metabolite STD

# Glycine

mz_glycine<-c(204,205)
rt_glycine<-c(630,660)

chr_glycine <- chromatogram(raw_peaks, mz = mz_glycine, rt = rt_glycine)
chromPeaks(chr_glycine)

sample_colors <- group_color[chr_glycine$sample_group]
plot(chr_glycine, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_glycine)[, "sample"]],
     peakBg = NA)


# Isoleucine and leucine?

mz_isoleucine<-c(232,233)
rt_isoleucine<-c(650, 700)

chr_isoleucine <- chromatogram(raw_peaks, mz = mz_isoleucine, rt = rt_isoleucine)
chromPeaks(chr_isoleucine)

sample_colors <- group_color[chr_isoleucine$sample_group]
plot(chr_isoleucine, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_isoleucine)[, "sample"]],
     peakBg = NA)


# lactate

mz_lactate<-c(193,194)
rt_lactate<-c(470,510)

chr_lactate <- chromatogram(raw_peaks, mz = mz_lactate, rt = rt_lactate)
chromPeaks(chr_lactate)

sample_colors <- group_color[chr_lactate$sample_group]
plot(chr_lactate, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_lactate)[, "sample"]],
     peakBg = NA)


# laurate - unclear RT, three peaks present

mz_laurate<-c(215,216)
rt_laurate<-c(920,1005)

chr_laurate <- chromatogram(raw_peaks, mz = mz_laurate, rt = rt_laurate)
chromPeaks(chr_laurate)

sample_colors <- group_color[chr_laurate$sample_group]
plot(chr_laurate, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_laurate)[, "sample"]],
     peakBg = NA)


# lysine probably 1130 sec peak, sister peak at 1070

mz_lysine<-c(174,175)
rt_lysine<-c(1050,1150)

chr_lysine <- chromatogram(raw_peaks, mz = mz_lysine, rt = rt_lysine)
chromPeaks(chr_lysine)

sample_colors <- group_color[chr_lysine$sample_group]
plot(chr_lysine, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_lysine)[, "sample"]],
     peakBg = NA)


# malate - finding two peaks, unsure which one's the right one

mz_malate<-c(233,234)
rt_malate<-c(820,880)

chr_malate <- chromatogram(raw_peaks, mz = mz_malate, rt = rt_malate)
chromPeaks(chr_malate)

sample_colors <- group_color[chr_malate$sample_group]
plot(chr_malate, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_malate)[, "sample"]],
     peakBg = NA)


# mannitol - completely confusing

mz_mannitol<-c(319,320)
rt_mannitol<-c(1100,1170)

chr_mannitol <- chromatogram(raw_peaks, mz = mz_mannitol, rt = rt_mannitol)
chromPeaks(chr_mannitol)

sample_colors <- group_color[chr_mannitol$sample_group]
plot(chr_mannitol, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_mannitol)[, "sample"]],
     peakBg = NA)


# mannose - not clear either

mz_mannose<-c(203,204)
rt_mannose<-c(1100,1170)

chr_mannose <- chromatogram(raw_peaks, mz = mz_mannose, rt = rt_mannose)
chromPeaks(chr_mannose)

sample_colors <- group_color[chr_mannose$sample_group]
plot(chr_mannose, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_mannose)[, "sample"]],
     peakBg = NA)


# myoinositol

mz_myoinositol<-c(305,306)
rt_myoinositol<-c(1210,1250)

chr_myoinositol <- chromatogram(raw_peaks, mz = mz_myoinositol, rt = rt_myoinositol)
chromPeaks(chr_myoinositol)

sample_colors <- group_color[chr_myoinositol$sample_group]
plot(chr_myoinositol, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_myoinositol)[, "sample"]],
     peakBg = NA)


# nacetylglucosamine - two peaks, unclear which one is NAG

mz_nacetylglucosamine<-c(319,320)
rt_nacetylglucosamine<-c(1200,1250)

chr_nacetylglucosamine <- chromatogram(raw_peaks, mz = mz_nacetylglucosamine, rt = rt_nacetylglucosamine)
chromPeaks(chr_nacetylglucosamine)

sample_colors <- group_color[chr_nacetylglucosamine$sample_group]
plot(chr_nacetylglucosamine, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_nacetylglucosamine)[, "sample"]],
     peakBg = NA)


# phenylalanine

mz_phenylalanine<-c(221,222)
rt_phenylalanine<-c(850,960)

chr_phenylalanine <- chromatogram(raw_peaks, mz = mz_phenylalanine, rt = rt_phenylalanine)
chromPeaks(chr_phenylalanine)

sample_colors <- group_color[chr_phenylalanine$sample_group]
plot(chr_phenylalanine, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_phenylalanine)[, "sample"]],
     peakBg = NA)


# proline

mz_proline<-c(216,217)
rt_proline<-c(660,770)

chr_proline <- chromatogram(raw_peaks, mz = mz_proline, rt = rt_proline)
chromPeaks(chr_proline)

sample_colors <- group_color[chr_proline$sample_group]
plot(chr_proline, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_proline)[, "sample"]],
     peakBg = NA)


# ribose

mz_ribose<-c(203,206)
rt_ribose<-c(900,1500)

chr_ribose <- chromatogram(raw_peaks, mz = mz_ribose, rt = rt_ribose)
chromPeaks(chr_ribose)

sample_colors <- group_color[chr_ribose$sample_group]
plot(chr_ribose, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_ribose)[, "sample"]],
     peakBg = NA)


# serine - peak picking seems quite shit

mz_serine<-c(204,205)
rt_serine<-c(730,770)

chr_serine <- chromatogram(raw_peaks, mz = mz_serine, rt = rt_serine)
chromPeaks(chr_serine)

sample_colors <- group_color[chr_serine$sample_group]
plot(chr_serine, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_serine)[, "sample"]],
     peakBg = NA)


# succinate

mz_succinate<-c(262,263)
rt_succinate<-c(650,760)

chr_succinate <- chromatogram(raw_peaks, mz = mz_succinate, rt = rt_succinate)
chromPeaks(chr_succinate)

sample_colors <- group_color[chr_succinate$sample_group]
plot(chr_succinate, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_succinate)[, "sample"]],
     peakBg = NA)


# sucrose includes maltose and trehalose peaks

mz_sucrose<-c(361,362)
rt_sucrose<-c(1470,1600)

chr_sucrose <- chromatogram(raw_peaks, mz = mz_sucrose, rt = rt_sucrose)
chromPeaks(chr_sucrose)

sample_colors <- group_color[chr_sucrose$sample_group]
plot(chr_sucrose, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_sucrose)[, "sample"]],
     peakBg = NA)


# threonine

mz_threonine<-c(291,292)
rt_threonine<-c(730,780)

chr_threonine <- chromatogram(raw_peaks, mz = mz_threonine, rt = rt_threonine)
chromPeaks(chr_threonine)

sample_colors <- group_color[chr_threonine$sample_group]
plot(chr_threonine, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_threonine)[, "sample"]],
     peakBg = NA)


# thymine

mz_thymine<-c(255,257)
rt_thymine<-c(700,800)

chr_thymine <- chromatogram(raw_peaks, mz = mz_thymine, rt = rt_thymine)
chromPeaks(chr_thymine)

sample_colors <- group_color[chr_thymine$sample_group]
plot(chr_thymine, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_thymine)[, "sample"]],
     peakBg = NA)


# valine

mz_valine<-c(188,190)
rt_valine<-c(550,700)

chr_valine <- chromatogram(raw_peaks, mz = mz_valine, rt = rt_valine)
chromPeaks(chr_valine)

sample_colors <- group_color[chr_valine$sample_group]
plot(chr_valine, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_valine)[, "sample"]],
     peakBg = NA)
#


#7 correspondence ####

#7.1 inital peak grouping ####
### initial correspondence analysis ###

## Initial peak grouping. Use sample_group as grouping variable
pdp <- PeakDensityParam(sampleGroups = raw_peaks$sample_group,
                             minFraction = 0.5, bw=2, binSize = 0.25)
grouped_peaks <- groupChromPeaks(raw_peaks, param = pdp)


#7.2 RT alignment ####
### Alignment ###
rtalign_params<-PeakGroupsParam(minFraction = 0.5, span = 0.4)

peaks_RTaligned <- adjustRtime(grouped_peaks, param = rtalign_params)


# ribitol
mz_ribitol<-c(217,217.4)
rt_ribitol<-c(1000,1020)

chr_ribitol <- chromatogram(peaks_RTaligned, mz = mz_ribitol, rt = rt_ribitol)
chromPeaks(chr_ribitol)

sample_colors <- group_color[chr_ribitol$sample_group]
plot(chr_ribitol, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_ribitol)[, "sample"]],
     peakBg = NA)
legend ("topright", inset=c(0,0), xpd = TRUE, legend=leg, 
        fill = group_color[chromPeaks(chr_ribitol)[,"column"]], cex = 0.8)

## Get the base peak chromatograms.
bpis <- chromatogram(raw_peaks, aggregationFun = "max", include = "none")
bpis_adj <- chromatogram(peaks_RTaligned, aggregationFun = "max", include = "none")
par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
plot(bpis, col = group_color[bpis$sample_group])
plot(bpis_adj, col = group_color[bpis_adj$sample_group])
## Plot also the difference of adjusted to raw retention time.
par(mfrow=c(1,1))
plotAdjustedRtime(peaks_RTaligned, col = group_color[peaks_RTaligned$sample_group])


# Comment on RT adjustment ------------------------------------------------

# It looks like the retention time adjustment actually worked really nicely!

# 7.3 test grouping ####
# Dry run

## Get BPC for ribitol and plot
mz_ribitol<-c(217,217.4)
rt_ribitol<-c(1000,1020)

bpc_ribitol <- chromatogram(peaks_RTaligned, mz = mz_ribitol, rt = rt_ribitol,
                            aggregationFun = "max")
par(mfrow = c(2, 1), mar = c(4, 4.3, 1, 0.5))
plot(bpc_ribitol)
highlightChromPeaks(peaks_RTaligned, mz = c(1540, 1550), whichPeaks = "apex_within")

pdp <- PeakDensityParam(sampleGroups = peaks_RTaligned$sample_group) # Get default parameter
plotChromPeakDensity(bpc_ribitol, param = pdp) # dry run to see if the parameters are good


## Get BPC for monosaccharides and plot
mz_monosaccharides<-c(319,322)
rt_monosaccharides<-c(1000,1200)

bpc_monosaccharides <- chromatogram(peaks_RTaligned, mz = mz_monosaccharides, rt = rt_monosaccharides,
                            aggregationFun = "max")
par(mfrow = c(2, 1), mar = c(4, 4.3, 1, 0.5))
plot(bpc_monosaccharides)
highlightChromPeaks(peaks_RTaligned, mz = c(1540, 1550), whichPeaks = "apex_within")

pdp <- PeakDensityParam(sampleGroups = peaks_RTaligned$sample_group) # Get default parameter
plotChromPeakDensity(bpc_monosaccharides, param = pdp) # dry run to see if the parameters are good



# 7.4 final grouping ####

my_peakdensityparams <- PeakDensityParam(sampleGroups = peaks_RTaligned$sample_group,
                        minFraction = 0.5, binSize = 0.25, bw = 2)
plotChromPeakDensity(bpc_monosaccharides, param = my_peakdensityparams)

grouped_peaks_RTaligned<-groupChromPeaks(peaks_RTaligned, param = my_peakdensityparams)


# 7.5 filling peaks

fmat <- featureValues(grouped_peaks_RTaligned, value = "into", method = "maxint")
sum(is.na(fmat)) # missing peaks before filling 

fpp <- FillChromPeaksParam(expandMz = 0.5, expandRt = 0.5, ppm = 20) # if you want to play with parameter

peaks_final <- fillChromPeaks(grouped_peaks_RTaligned, param = fpp)
sum(is.na(featureValues(peaks_final))) # missing peaks after rescuing

# filling reduced missing peaks from ~400,000 to 140,000 - WOW


#### plot final peaks of monosaccharides ####

setwd(results_dir)
tiff("BasePeakChromatogram_somesugars_216-219.tiff", res = 150, height = 12, width = 17.5, units = "cm")
mz_monosaccharides<-c(216,219)
rt_monosaccharides<-c(1050,1200)

bpc_monosaccharides <- chromatogram(peaks_final, mz = mz_monosaccharides, rt = rt_monosaccharides,
                                    aggregationFun = "max")
par(mfrow = c(2, 1), mar = c(4, 4.3, 1, 0.5))
#plot(bpc_monosaccharides)
highlightChromPeaks(peaks_RTaligned, mz = c(1540, 1550), whichPeaks = "apex_within")

pdp <- PeakDensityParam(sampleGroups = peaks_final$sample_group) # Get default parameter
plotChromPeakDensity(bpc_monosaccharides, param = pdp, col = group_color[peaks_final$sample_group]) 

dev.off()

setwd(results_dir)
tiff("BasePeakChromatogram_somesugars.tiff", res = 150, height = 12, width = 17.5, units = "cm")
mz_monosaccharides<-c(319,322)
rt_monosaccharides<-c(1000,1200)

bpc_monosaccharides <- chromatogram(peaks_final, mz = mz_monosaccharides, rt = rt_monosaccharides,
                                    aggregationFun = "max")
par(mfrow = c(2, 1), mar = c(4, 4.3, 1, 0.5))
#plot(bpc_monosaccharides)
highlightChromPeaks(peaks_RTaligned, mz = c(1540, 1550), whichPeaks = "apex_within")

pdp <- PeakDensityParam(sampleGroups = peaks_final$sample_group) # Get default parameter
plotChromPeakDensity(bpc_monosaccharides, param = pdp, col = group_color[peaks_final$sample_group]) 

dev.off()


mz_test<-c(146,149)
rt_test<-c(500,650)

bpc_test <- chromatogram(peaks_final, mz = mz_test, rt = rt_test,
                                    aggregationFun = "max")
plotChromPeakDensity(bpc_test, param = pdp, col = group_color[peaks_final$sample_group]) 
# Save result -------------------------------------------------------------

save(peaks_final, file="peaks_final_as XCMSnExp.RData")

# write MGF file for GNPS

writeMgfData(peaks_final, "peaks_final.mgf")
load("peaks_final.mgf")

# Extract intensity results
into_final <- featureValues(peaks_final, value = "into", method = "maxint")
#store the matrix to disk (avoid formating issues)
setwd(file.path(results_dir))
save(into_final, file = "df_exp_into.RData")

load("df_exp_into.RData")

# Extract mz and rt from the matrix 
mz <- featureValues(peaks_final, value = "mz")
head(mz)
rt <- featureValues(peaks_final, value = "rt")
head(rt)
# transform them
mz_long <- melt(mz)
head(mz_long)
colnames(mz_long) <- c("feature", "sample", "mz")
rt_long <- melt(rt)
head(rt_long)
colnames(rt_long) <- c("feature", "sample", "rt")
#combine
rt_mz <- merge(rt_long, mz_long, by = c("feature", "sample"))
head(rt_mz)
#export rt and mz info as tab delimited file
write.table(rt_mz, "exp_mzrt.csv", row.names = TRUE, sep="\t", dec = ".", col.names = TRUE)


# Quick Qualtity Control ####

par(mar=c(5,4,1,7))
boxplot(featureValues(peaks_final, value="into") +1, 
        col = group_color[peaks_final$sample_group], 
        log="y", las=2)
legend ("right", inset=c(-0.25,0), xpd = TRUE, legend=leg, fill = group_color)



# attempt GNPS ------------------------------------------------------------

filteredMs2Spectra <- featureSpectra(peaks_final, return.type = "MSpectra")
filteredMs2Spectra <- clean(filteredMs2Spectra, all = TRUE)
filteredMs2Spectra <- formatSpectraForGNPS(filteredMs2Spectra)




# Convert to xset ####

#Step 6: convert XCMSnExp back to xcmsSet (XCMSnExp are not compatible with all function and other packages)
xset <- as(peaks_final, "xcmsSet")
sampnames(xset) <- pData(peaks_final)$sample_name
sampclass(xset) <- pData(peaks_final)$sample_group
xset2 <- fillPeaks(xset) 
save(xset2, file = "SeaMet_final_asXset.RData")


# *** Savepoint *** ####
setwd("C:/Users/admin/ownCloud/Assemble/Data/SeaMet/GCMS data analysis/results/")
load("SeaMet_final_asXset.RData")

peaks_from_xset<-as(xset2, "xcmsnExp")
bpc_monosaccharides <- chromatogram(peaks_RTaligned, mz = mz_monosaccharides, rt = rt_monosaccharides,
                                    aggregationFun = "max")
par(mfrow = c(2, 1), mar = c(4, 4.3, 1, 0.5))
plot(bpc_monosaccharides)

#Plot mz deviation
plotQC(xset2, what="mzdevhist")

# plot rt deviation
plotQC(xset2, what="rtdevhist")

# mz deviation per sample
plotQC(xset2, what="mzdevsample")


# Quick PCA - but why?
sdThresh <- 4.0 ## Filter low-standard deviation rows for plot
data <- log(featureValues(peaks_final))+1
pca.result <- pca(data, nPcs=2)
pcap <- plotPcs(pca.result, type= "loadings",col = group_color[peaks_final$sample_group])
legend ("topright", inset=c(0,0), xpd = TRUE, legend=leg, fill = group_color)




# 8. Deconvolution --------------------------------------------------------

# create CAMERA object
an <- xsAnnotate(xset2, sample=seq(1,length(sampnames(xset2))))

#peak grouping by rt
an <- groupFWHM(an, perfwhm = 0.2)

# find isotopic peaks
an <- findIsotopes(an, mzabs = 0.01)

# verify peak groups
an <- groupCorr(an, graphMethod="lpc", calcIso = TRUE, calcCiS = TRUE, calcCaS = TRUE, cor_eic_th=0.6)

xsaFA <- findAdducts(an, polarity="positive")

peaklist<-getPeaklist(an, intval="into")

peaks_for_mgf<-writeMSData(xset2)
MGF_peaks<-writeMgfData(xset2)


write.csv(peaklist, file="peaklist_06March2021.csv")

# for GNPS, get pcgroup to the front and label as feature

peaks<-peaks[order(peaks$pcgroup),]
peaks[is.na(peaks)]<-0

head(peaks)


# write mgf with pseudospectra from xsannotate

print_COM<-"COM=Finland SeaMet pseudospectra exported from xsAnnotate"

print_start<-"BEGIN IONS"
print_end<-"END IONS"



plotEICs(an, pspec=c(13,18,23,24,25,33, 44, 227, 503), maxlabel = 3)

plotEICs(an, pspec=2, maxlabel=0)

plotPsSpectrum(an,pspec=1, maxlabel = 3)

an@psSamples

peaklist %>%
  subset(pcgroup == c(18)) %>%
  select(c(mz,intensity))
max(peaks$pcgroup)

psp.peaks <- getpspectra(xsaFA, 1) 

pspec2metfusion(object=xsaFA, filedir = "annotated_pseudospectra/")


## get feature definitions and intensities
featuresDef <- featureDefinitions(xsaFA)
featuresIntensities <- featureValues(processedData, value = "into")

## generate data table
dataTable <- merge(featuresDef, featuresIntensities, by = 0, all = TRUE)
dataTable <- dataTable[, !(colnames(dataTable) %in% c("peakidx"))]