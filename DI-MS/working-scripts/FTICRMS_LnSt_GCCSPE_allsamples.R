setwd("/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/4_FTICRMS/")
load("./analysis/RData/RData_2021_02_22_allsamples.RData")


# 1: Install packages --------------------------------------------------------
library(BiocStyle)
library(xcms)
library(faahKO)
library(pander)
library(RColorBrewer)
library(magrittr)
library(pheatmap)
library(MSnbase)
library(msdata)
library(png)
library(IPO)
library(tidyr)
library(detect)
library(devtools)
library(MetMSLine)
library(pcaMethods)
library(statTarget)
library(randomForest)
library(rlist)
library(purrr)
library(wesanderson)
library(ggplot2)
library(reshape2)
library(extrafont)
library(Rmisc)
library(edgeR)
library(limma)
library(mixOmics)
library(HTSFilter)
library(rstatix)
library(reshape2)
library(scales)
library(data.table)
library(remotes)
library(ggridges)
library(gridExtra)
library(tidyverse)
library(ggpubr)
library(viridis)
library(lemon)
library(gtable)
library(grid)
library(colorspace)
library(MassSpecWavelet)
library(SummarizedExperiment)

#2. Import and inspect data --------------------------------------------------------
#get file paths to mzML files
fp <- dir(path = "./data/MAHALE", 
          all.files = FALSE, 
          full.names = TRUE)
fp <- fp[grep("\\.mzML", fp)]

#create phenodata data.frame

name <- basename(fp) %>%
    sub("MGC21006\\d{2}_LSt_MB_RP_FTMS_", "", .) %>% 
    sub("MGC21006\\d{2}_LSt_RPLSt_FTMS_", "", .) %>% 
    sub(".mzML", "", .)

sample_type <- name %>% 
    sub(".*Blank.*", "blank", .) %>% 
    sub("EC_STD.*", "standard mix extraction",.) %>% 
    sub("STD_Mix_all.*", "standard mix all",.) %>% 
    sub("STD_Mix_pH2.*", "standard mix all pH2",.) %>% 
    sub("STD_Mix_1.*", "vitamin B1", .) %>% 
    sub("STD_Mix_2.*", "vitamin B2", .) %>% 
    sub("STD_Mix_3.*", "N metabolites", .) %>% 
    sub("STD_Mix_4.*", "neutral oligos", .) %>% 
    sub("STD_Mix_6.*", "neutral monomers", .) %>% 
    sub("STD_mann_acid_1.*", "tri-mannuronic acid", .) %>% 
    sub("Wash.*", "wash", .)

medium <- name %>% 
    sub(".*_1A_.*", "MQ", .) %>% 
    sub(".*_2A_.*", "ASW", .) %>% 
    sub(".*\\d+.*", "NA", .)

concentration <- name %>% 
    sub(".*Blank.*", 0, .) %>% 
    sub(".*STD_5g.*", 5, .) %>% 
    sub(".*STD_1g.*", 1, .) %>% 
    sub(".*STD_0_1g.*", 0.1, .) %>% 
    sub("[[:alpha:]]+.*", "NA", .)

k = name %>% sub("_\\d{5}$", "", .) %>% 
    unique() %>% 
    length()

n = round(length(fp) / k, 0)


rep <- rep(c(1:n), k)
rep <- c(rep, 4, 5)


pd <- data.frame(name = name,
                 sample_type = sample_type,
                 medium = medium,
                 concentration = concentration,
                 replicate = rep,
                 stringsAsFactors = FALSE)

#load data
data <- readMSData(files = fp,
                    pdata = new("NAnnotatedDataFrame",
                                pd),
                    mode = "onDisk")

#set ms level to 1... for some reason all have ms level set to 2??
data@featureData@data$msLevel <- 1
data@featureData@data$isolationWindowLowerOffset <- NA
data@featureData@data$isolationWindowUpperOffset <- NA
data@featureData@data$isolationWindowTargetMZ <- NA

#3: Create initial output directories -------------------------------------
dir.create("./analysis",
           showWarnings = FALSE)
dir.create("./analysis/RData",
           showWarnings = FALSE)

#4: Pick peaks (mass spec wavelet) -------------------------------------
#msw default param:
#   snthresh: 3 
#   verboseColumns: FALSE 
#   scales: 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,36,40,44,48,52,56,60,64 
#   nearbyPeak: TRUE 
#   peakScaleRange: 5 
#   ampTh: 0.01 
#   minNoiseLevel: 0.003333333 
#   ridgeLength: 24 
#   peakThr:  
#   tuneIn: FALSE 

msw <- MSWParam(nearbyPeak = FALSE, 
                snthresh = 0.01,
                ampTh= 0.001,
                minNoiseLevel= 0.01,
                peakScaleRange = 1,
                SNR.method = "data.mean",
                scales = c(1),
)

#msw <- MSWParam()

data <- findChromPeaks(data, param = msw)
peaks <- chromPeaks(data)

#5: Group peaks -------------------------------------

## Using default settings but define sample group assignment
mzc <- MzClustParam(sampleGroups = data$sample_type,
                    minFraction = 0.01)
data <- groupChromPeaks(data, param = mzc)


#6: Fill missing peaks -------------------------------------

data <- fillChromPeaks(data, 
                       param = FillChromPeaksParam())

featVal <- featureValues(data)
featDef <- featureDefinitions(data)


#7: Isotope and adduct detection ----
#set as xcmsSet object
xset <- as(data, "xcmsSet")
sampnames(xset) <- pData(data)$name
sampclass(xset) <- pData(data)$sample_type
##create xsannotate object
an <- xsAnnotate(xset)
##Annotate isotope peaks
an <- findIsotopes(an, 
                   mzabs=0.01)
##Find adducts
an <- findAdducts(an, 
                  polarity="negative")
##Get peak list
pl <-getPeaklist(an)

#8: Collapse isotopes ----
setDT(pl)
#split out features without an isotope detected
noiso <- pl[pl$isotopes=="",]
iso <- pl[!pl$isotopes=="",]
#make column for the isotope group
iso$isotope_group <- iso$isotopes %>% 
    sub("\\[M.*", "", .)
#order isotopes within each group correctly
iso$isotope_number <- iso$isotopes %>% 
    sub(".*\\[M\\].*", "0", .) %>% 
    sub(".*\\[M\\+", "", .) %>% 
    sub("\\].*", "", .) %>% 
    as.numeric()
iso <- iso[order(isotope_group, 
                 isotope_number),]
#get concatenated list of isotopes per group
iso_concat <- iso[, 
                  list(isotopes = paste(isotopes, 
                                        collapse = ', ')), 
                  by = isotope_group]
#remove duplicates within each isotope group (will keep [M] isotope)
#because of ordering
iso <- unique(iso, 
              by = "isotope_group")
#merge to get concatenated isotope lists
iso <- merge(iso,
             iso_concat,
             by = "isotope_group")
#clean up df
iso <- iso %>% 
    select(-c("isotope_group",
              "isotope_number",
              "isotopes.x"))
names(iso)[names(iso) == 'isotopes.y'] <- 'isotopes'

#merge features with and without isotopes
pl_old <- pl 
pl <- rbind.fill(iso,
                 noiso)
rm(iso,
   noiso,
   iso_concat)

#9: Annotate with predicted oligosaccharides  ----

#predicted with: 
#sugarMassesPredict.py -dp 1 6 -p 1 -m sulphate carboxyl deoxy anhydrobridge nacetyl -n 1 -ld D -oh 0 -b 0 -i neg -s 92 2000 -o script-prediction_1.txt

pred1 <- fread("script-prediction_1.txt")
pred1$index <- NULL

good_names <- c("glucose",
                "fucose",
                "xylose",
                "laminaribiose",
                "mannotriose",
                "cellotetraose",
                "laminarihexaose",
                "mannuronic acid",
                "di-mannuronic acid",
                "tri-mannuronic acid",
                "penta-mannuronic acid",
                "tetra-guluronic acid",
                "hexa-guluronic acid",
                "k-carrageenan DP4",
                "N-acetylglucosamine")
pred_names <- c("hex-1",
                "hex-1-deoxy-1",
                "pent-1",
                "hex-2",
                "hex-3",
                "hex-4",
                "hex-6",
                "hex-1-carboxyl-1",
                "hex-2-carboxyl-2",
                "hex-3-carboxyl-3",
                "hex-5-carboxyl-5",
                "hex-4-carboxyl-4",
                "hex-6-carboxyl-6",
                "hex-4-sulphate-2-anhydrobridge-2",
                "hex-1-nacetyl-1")
names(pred_names) <- good_names

#filter to only keep standards
pred2 <- pred1 %>% 
    filter(name %in% pred_names)

#add common names
pred2 <- cbind(common_name = names(pred_names[order(match(pred_names,pred2$name))]),
               pred2)

#make long format
pred2_wide <- pred2
pred2 <- gather(pred2_wide, 
                key = "ion", 
                value = "mz", 
                -common_name,
                -name,
                -dp,
                -mass,
                -formula
)
pred2 <- na.omit(pred2)

pred2$name <- NULL

#match with overlap (+-1 ppm)
pred2$mzmin <- pred2$mz - (pred2$mz / 1000000)
pred2$mzmax <- pred2$mz + (pred2$mz / 1000000)
setDT(pl)
setDT(pred2)
setkey(pred2, mzmin, mzmax)
pl_annot <- foverlaps(pl,
                      pred2)


#10 : add vitamins / others not predicted in script ----
H_diff <- 1.007277
Cl_diff <- 34.9694


extra_common_name = c("glucosamine", 
                      "GABA",
                      "thiamin",
                      "riboflavin",
                      "pantothenic acid",
                      "pyridoxine",
                      "niacin",
                      "biotin",
                      "folate")


pred_extra <- data.frame(common_name = extra_common_name,
                         dp = c(1,1,rep(NA, length(extra_common_name)-2)),
                         mass = c(179.079373,
                                  103.063329,
                                  300.081161,
                                  376.138286,
                                  219.110674,
                                  169.073894,
                                  123.032029,
                                  244.088165,
                                  441.139683),
                         formula = c("C6H13NO5",
                                     "C4H9NO2",
                                     "C12H17ClN4OS",
                                     "C17H20N4O6",
                                     "C9H17NO5",
                                     "C8H11NO3",
                                     "C6H5NO2",
                                     "C10H16N2O3S",
                                     "C19H19N7O6"))
pred_extra$`[M-H]` <- pred_extra$mass - H_diff
pred_extra$`[M-2H]` <- (pred_extra$mass - 2*H_diff) / 2
pred_extra$`[M+Cl]` <- pred_extra$mass + Cl_diff

pred_extra_wide <- pred_extra
pred_extra <- gather(pred_extra_wide, 
                key = "ion", 
                value = "mz", 
                -common_name,
                -dp,
                -mass,
                -formula
)

#match with overlap
pred_extra$mzmin <- pred_extra$mz - (pred_extra$mz/1000000)
pred_extra$mzmax <- pred_extra$mz + (pred_extra$mz/1000000)
setDT(pl_annot)
setDT(pred_extra)
setkey(pred_extra, mzmin, mzmax)
pl_annot2 <- foverlaps(pl,
                       pred_extra)


#merge
pl_annot2$mzmin <- NULL
pl_annot2$mzmax <- NULL
names(pl_annot2)[names(pl_annot2) == "mz"] <- "theoretical_mz"
names(pl_annot2)[names(pl_annot2) == "i.mz"] <- "mz"
names(pl_annot2)[names(pl_annot2) == "i.mzmin"] <- "mzmin"
names(pl_annot2)[names(pl_annot2) == "i.mzmax"] <- "mzmax"


pl_annot$mzmin <- NULL
pl_annot$mzmax <- NULL
names(pl_annot)[names(pl_annot) == "mz"] <- "theoretical_mz"
names(pl_annot)[names(pl_annot) == "i.mz"] <- "mz"
names(pl_annot)[names(pl_annot) == "i.mzmin"] <- "mzmin"
names(pl_annot)[names(pl_annot) == "i.mzmax"] <- "mzmax"

pl_annot.merge <- rbind(pl_annot,
                        pl_annot2)

setDT(pl_annot.merge)
pl_annot.merge <- pl_annot.merge[order(common_name),]

pl_annot.merge <- unique(pl_annot.merge, 
              by = "mz")

#11 : adducts ----
adducts <- pl_annot.merge[pl_annot.merge$adduct != ""]
adducts <- adducts %>% 
    select(c(mz, adduct))
adduct_masses <- adducts$adduct %>% 
    gsub("\\[\\d{0,1}M.*?\\]\\d{0,1}-\\s{0,1}", "", .) %>% 
    sub("^\\s", "", .) %>% 
    strsplit(.,  "  ") %>% 
    unlist %>% 
    as.numeric()

adduct_masses.df <- data.frame(mass = adduct_masses,
                               massmin = adduct_masses - 0.001,
                               massmax = adduct_masses + 0.001)  

setDT(adduct_masses.df)  
adduct_masses.df <- unique(adduct_masses.df,
                           by = "mass")

    
pred_adduct <- rbind.fill(pred2_wide, pred_extra_wide)
pred_adduct <- pred_adduct %>% 
    select(c("common_name",
             "dp",
             "mass",
             "formula"))

pred_adduct$massmin <- pred_adduct$mass - 0.001
pred_adduct$massmax <- pred_adduct$mass + 0.001

setDT(pred_adduct)    

setkey(adduct_masses.df, massmin, massmax)    
pred_adduct_matched <- foverlaps(pred_adduct,
                                 adduct_masses.df)
#only matched metabolites already annotated

dir.create("./analysis/analysis_tables",
           showWarnings = FALSE)


#12: write peak list tables and feature summaries to file ----

pl_annot.merge$pcgroup <- NULL

fwrite(pl_annot.merge,
       "./analysis/analysis_tables/peaklist_annotated_allsamples_20210222.txt",
       sep = "\t")

pl_annot.merge_matchedonly <- pl_annot.merge[pl_annot.merge$common_name != "NA",]

fwrite(pl_annot.merge_matchedonly,
       "./analysis/analysis_tables/peaklist_annotated_allsamples_matchedonly_20210222.txt",
       sep = "\t")


#get per sample counts
allFeatSum <- featureSummary(data,
                          group = pd$sample_type,
                          perSampleCounts = TRUE)

allFeatMz <- featDef@listData$mzmed
allFeatNames <- featDef@rownames
allFeatDef.df <- data.frame(featnames = allFeatNames,
                            mz = allFeatMz)

matchedFeatMz <- pl_annot.merge_matchedonly$mz
matchedFeatNames <- allFeatDef.df$allFeatDef.df$featnames[allFeatDef.df$mz %in% matchedFeatMz]
matchedFeatDef.df <- data.frame(featnames = matchedFeatNames,
                                mz = matchedFeatMz)

matchedFeatSum <- allFeatSum[rownames(allFeatSum) %in% matchedFeatNames,]

matchedFeatSum.df <- as.data.frame(matchedFeatSum)
matchedFeatSum.df$featnames <- rownames(matchedFeatSum)


setDT(pl_annot.merge_matchedonly); setDT(matchedFeatDef.df)
setkey(pl_annot.merge_matchedonly, mz)
setkey(matchedFeatDef.df, mz)
pl_annot.merge_matchedonly2 <- merge(pl_annot.merge_matchedonly,
                                     matchedFeatDef.df)

setDT(matchedFeatSum.df)
setkey(matchedFeatSum.df, featnames)
setkey(pl_annot.merge_matchedonly2, featnames)
matchedFeatSumAnnot.df <- merge(matchedFeatSum.df,
                                pl_annot.merge_matchedonly2) 

extraCol <- names(matchedFeatSumAnnot.df)[
    grepl("multi_count|multi_perc|_rsd",
          names(matchedFeatSumAnnot.df))]
matchedFeatSumAnnot.df <- matchedFeatSumAnnot.df %>% 
    select(-extraCol)

extraCol <- setdiff(setdiff(names(pl_annot.merge_matchedonly2), 
                            names(pred2)),
                    c("isotopes",
                      "featnames"))

matchedFeatSumAnnot.df <- matchedFeatSumAnnot.df %>% 
    select(-extraCol)

fwrite(matchedFeatSumAnnot.df,
       "./analysis/analysis_tables/featSum_annotated_allsamples_matchedonly_20210222.txt",
       sep = "\t")

matchedFeatSumAnnotSamples.df <- matchedFeatSumAnnot.df
matchedFeatSumAnnotSamples.df <- matchedFeatSumAnnotSamples.df %>% 
    select(c("common_name", "ion", basename(fp)))
matchedFeatSumAnnotSamples.df <- matchedFeatSumAnnotSamples.df %>% 
    gather(key = "sample", value = "count", -common_name, -ion)

matchedFeatSumAnnotSamples.df$sample_type <- matchedFeatSumAnnotSamples.df$sample %>% 
    sub("MGC21006\\d{2}_LSt_RPLSt_FTMS_", "", .) %>% 
    sub(".*Blank.*", "blank", .) %>% 
    sub("EC_STD.*", "standard mix extraction",.) %>% 
    sub("STD_Mix_all.*", "standard mix all",.) %>% 
    sub("STD_Mix_pH2.*", "standard mix all pH2",.) %>% 
    sub("STD_Mix_1.*", "vitamin B1", .) %>% 
    sub("STD_Mix_2.*", "vitamin B2", .) %>% 
    sub("STD_Mix_3.*", "N metabolites", .) %>% 
    sub("STD_Mix_4.*", "neutral oligos", .) %>% 
    sub("STD_Mix_6.*", "neutral monomers", .) %>% 
    sub("STD_mann_acid_1.*", "tri-mannuronic acid", .) %>% 
    sub("Wash.*", "wash", .) %>% 
    sub("MGC21006\\d{2}_LSt_MB_RP_FTMS_", "", .)

matchedFeatSumAnnotSamples.df$medium <- matchedFeatSumAnnotSamples.df$sample %>% 
    sub("MGC21006\\d{2}_LSt_RPLSt_FTMS_", "", .) %>% 
    sub("MGC21006\\d{2}_LSt_MB_RP_FTMS_", "", .) %>% 
    sub(".*_1A_.*", "MQ", .) %>% 
    sub(".*_2A_.*", "ASW", .) %>% 
    sub(".*\\d+.*", "NA", .)

matchedFeatSumAnnotSamples.df$concentration <- matchedFeatSumAnnotSamples.df$sample %>% 
    sub("MGC21006\\d{2}_LSt_RPLSt_FTMS_", "", .) %>% 
    sub("MGC21006\\d{2}_LSt_MB_RP_FTMS_", "", .) %>% 
    sub(".*Blank.*", 0, .) %>% 
    sub(".*STD_5g.*", 5, .) %>% 
    sub(".*STD_1g.*", 1, .) %>% 
    sub(".*STD_0_1g.*", 0.1, .) %>% 
    sub("[[:alpha:]]+.*", "NA", .)

setDF(matchedFeatSumAnnotSamples.df)
matchedFeatSumAnnotSamplesSum.df <- matchedFeatSumAnnotSamples.df %>% 
    group_by(common_name,
             ion,
             sample_type, 
             medium,
             concentration) %>% 
    dplyr::summarise(count_sum=sum(count))

matchedFeatSumAnnotSamplesSum.df <- matchedFeatSumAnnotSamplesSum.df %>% 
    filter(count_sum > 0)


fwrite(matchedFeatSumAnnotSamplesSum.df,
       "./analysis/analysis_tables/featSum_annotated_allsamples_summary_matchedonly_20210222.txt",
       sep = "\t")


#13: plotting representation of results summary ----

metabolites <- matchedFeatSumAnnotSamplesSum.df$common_name[
    matchedFeatSumAnnotSamplesSum.df$sample_type == "standard mix all" &
        matchedFeatSumAnnotSamplesSum.df$count_sum >= 2 |
        matchedFeatSumAnnotSamplesSum.df$sample_type == "neutral oligos" &
        matchedFeatSumAnnotSamplesSum.df$count_sum >= 2] %>% 
    unique()

notMetabolites <- matchedFeatSumAnnotSamplesSum.df %>% 
    filter(sample_type != "standard mix extraction" & count_sum >1) %>% 
    group_by(common_name, sample_type) %>% 
    dplyr::summarise(n = n()) %>% 
    group_by(common_name) %>% 
    dplyr::summarise(n = n()) %>% 
    filter(n > 4) %>% 
    select(common_name)

notMetabolites <- notMetabolites$common_name

metabolites <- metabolites[!metabolites %in% notMetabolites]

res <- matchedFeatSumAnnotSamplesSum.df %>% 
    filter(common_name %in% metabolites)

x <- data.frame(common_name = c("k-carrageenan DP4", 
                                "mannuronic acid"),
                sample_type = c("standard mix extraction",
                                "standard mix extraction"),
                concentration = c(0.1, 0.1),
                medium = c("MQ", "MQ"))

res <- rbind.fill(res,
                  x)

res$common_name <- factor(res$common_name,
                          levels= c("glucose",
                                    "laminaribiose",
                                    "mannotriose",
                                    "cellotetraose",
                                    "laminarihexaose",
                                    "mannuronic acid",
                                    "k-carrageenan DP4",
                                    "N-acetylglucosamine",
                                    "pantothenic acid",
                                    "pyridoxine",
                                    "biotin"))

res$common_name <- factor(res$common_name,
                          levels = rev(levels(res$common_name)))

res$concentration <- as.character(res$concentration)
res$concentration <- factor(res$concentration,
                            levels = c("0.1",
                                       "1",
                                       "5"))

res$medium <- factor(res$medium,
                     levels = c("MQ", "ASW"))


tiff("./analysis/plots/peakcount_summary_v1.tiff",
     units = "in",
     res = 300,
     width = 12,
     height = 6)

svg("./analysis/plots/peakcount_summary_v1.svg",
     width = 12,
     height = 6)

ggplot(res[res$sample_type == "standard mix extraction",],
       aes(x = concentration,
           y = common_name,
           colour = medium,
           size = count_sum,
           group = medium)) +
    scale_colour_manual(values = c("#7DB1DE",
                                   "#0B4150"),
                        name = "Medium",
                        labels = c("MilliQ-water",
                                   "Artificial sea water")) +
    geom_point(stat = "identity",
               position = position_dodge(width = 0.3)) +
    theme_classic() +
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 12)) +
    labs(x= "Initial amount (ug)",
         y = "Metabolite")
dev.off()    


#PLOT RANDOM SPECTRUM FOR DEFENCE (REPRESENT FT-ICR-MS)    

tiff("/Users/margotbligh/Google_Drive/MPI_Masters/MSc_thesis/Thesis/Defence/Figures/ft-icr-ms.tiff",
     res = 600,
     units = "in",
     width = 4.01,
     height = 2.41)
ggplot(y, aes(x = mz, y = intensity)) + 
    geom_segment(aes(x=mz, xend=mz, y=0, yend=intensity)) + 
    scale_x_continuous(name = expression(paste(italic("m/z"))),
                       limits = c(50,750),
                       expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(text = element_text(family = "Avenir", size = 16),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.line = element_blank(),
          axis.ticks.y = element_blank())
dev.off()









