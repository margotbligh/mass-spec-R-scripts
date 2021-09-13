setwd("/Users/margotbligh/Google_Drive/MPI_PhD/Lab-things/alpha-mannan/Polaribacter_Hel1-33-78_enzymes/202109_GH99_PLx_GH92")
load("./analysis/RData/RData_20210910.RData")

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
library(ggsci)

#2. Import and inspect MS data --------------------------------------------------------
#get file paths to mzML files
fp <- dir(path = "MS31_20210907", 
          all.files = FALSE, 
          full.names = TRUE)

#create phenodata data.frame
#each sample must have a unique name!
pd <- data.frame(name = basename(fp) %>%
                     sub("MS31_20210907_", "", .) %>%
                     sub("_\\d\\d.mzML", "", .),
                 sample_type = basename(fp) %>% 
                     sub(".*_[2345][ab]_.*", "digest", .) %>% 
                     sub(".*_1[ab]_.*", "neg", .) %>% 
                     sub(".*blank.*", "solvent blank", .) %>% 
                     sub(".*std.*", "standard", .),
                 enzyme = basename(fp) %>% 
                     sub(".*_2[ab]_.*", "PLx", .) %>% 
                     sub(".*_3[ab]_.*", "GH99-sulphatase", .) %>% 
                     sub(".*_5[ab]_.*", "GH92 + GH99-sulphatase", .) %>%
                     sub("MS31.*", "NA", .),
                 spiked = basename(fp) %>% 
                     sub(".*_[1235]a_.*", "y", .) %>% 
                     sub(".*_[1235]b_.*", "n", .) %>% 
                     sub("MS31_.*", "NA", .),
                 group = basename(fp) %>% 
                     sub(".*_1[ab]_.*", "no digest", .) %>% 
                     sub(".*_2[ab]_.*", "PLx digest", .) %>% 
                     sub(".*_3[ab]_.*", "GH99-sulphatase digest", .) %>% 
                     sub(".*_5[ab]_.*", "GH92 + GH99-sulphatase digest", .) %>% 
                     sub(".*blank.*", "solvent blank", .) %>% 
                     sub(".*std.*", "carrageenan standard", .),
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
dir.create("./analysis/processing_plots",
           showWarnings = FALSE)

#4: Plot TIC ----
dir.create("./analysis/processing_plots/tic",
           showWarnings = FALSE)

pal_group <- pal_npg("nrc", alpha = 1)(length(unique(pd$group)))
names(pal_group) <- unique(pd$group)

#plot the tic as boxplot
tc <- split(tic(all_data), 
            f = fromFile(all_data))
cairo_pdf("./analysis/processing_plots/tic/tic_boxplot.pdf",
          family = "Avenir",
          width = 12,
          height = 9)
par(mar=c(9,5,1,1))
boxplot(tc, 
        col = pal_group[all_data$group],
        ylab = "intensity", 
        main = "total ion current",
        names = all_data$name,
        las=2,
        cex.axis = 0.8)
dev.off()


#5: Plot kappa-carrageenan EIC ----
k2.mh <- 403.055189579909
k4.m2h <- 394.04990707990896
k6.m3h <- 391.04814624657564
kik6.m4h <- 313.0234948299089

error <- 0.001

chr.k2 <- chromatogram(all_data,
                       mz = c(k2.mh-error,
                              k2.mh+error))
chr.k4 <- chromatogram(all_data,
                       mz = c(k4.m2h-error,
                              k4.m2h+error))
chr.k6 <- chromatogram(all_data,
                       mz = c(k6.m3h-error,
                              k6.m3h+error))
chr.kik6 <- chromatogram(all_data,
                         mz = c(kik6.m4h-error,
                                kik6.m4h+error))

chr_kappa <- list(chr.k2, chr.k4, chr.k6, chr.kik6)
chr_kappa.df <- data.frame(sample = as.character(),
                           group = as.character(),
                           ion = as.character(),
                           rt = as.numeric(),
                           intensity = as.numeric())
kappa.ions <- c(k2.mh, k4.m2h, k6.m3h,kik6.m4h)

names(kappa.ions) <- c("kappa-2 [M-H]-",
                       "kappa-4 [M-2H]2-",
                       "kappa-6 [M-3H]3-",
                       "kappa-iota-kappa-6 [M-4H]4-")

for (i in 1:length(kappa.ions)){
    ion.tmp = names(kappa.ions)[i]
    mz.tmp = kappa.ions[i]
    for (j in 1:length(pd$name)){
        rt = chr_kappa[[i]][[j]]@rtime/60
        intensity = chr_kappa[[i]][[j]]@intensity
        mz = rep(mz.tmp, length(rt))
        sample = rep(pd$name[j], length(rt))
        group = rep(pd$group[j], length(rt))
        ion = rep(ion.tmp, length(rt))
        temp <- data.frame(sample = sample,
                           group = group,
                           ion = ion,
                           rt = rt,
                           intensity = intensity)
        chr_kappa.df <- rbind(temp,
                              chr_kappa.df)
    }
    rm(temp, sample, group, rt, intensity, ion)
}

#set NA to 0
chr_kappa.df[is.na(chr_kappa.df)] <- 0

#make directory
dir.create("./analysis/analysis_plots/kappa-carrageenan-eic",
           showWarnings = FALSE)

#make y axis labelling function
scientific_function <- function(x) {
    text <- gsub("E0", "", gsub("e\\+0", "E", scales::scientific_format()(x)))
    text
}

#make strip text labelling function
ion_names <- list(
    'kappa-2 [M-H]-'= expression(paste(kappa, 
                                       "-carrageenan DP2 [M-H]"^-{})),
    "kappa-4 [M-2H]2-"=expression(paste(kappa, 
                                        "-carrageenan DP4 [M-2H]"^-{2})),
    "kappa-6 [M-3H]3-"=expression(paste(kappa, 
                                        "-carrageenan DP6 [M-3H]"^-{3})),
    "kappa-iota-kappa-6 [M-4H]4-"=expression(paste(kappa,"-",iota, "-", kappa, 
                                                   "-carrageenan DP6 [M-4H]"^-{4}))
)

ion_labeller <- function(variable,value){
    return(ion_names[value])
}

#rename standard group variable
chr_kappa.df$groupchr_kappa.df$group %>% 
    sub("standard", "oligosaccharide standard", .)

#plot eic
tiff("./analysis/analysis_plots/kappa-carrageenan-eic/facet-by-ion_v1.tiff",
     width = 9,
     height = 6,
     res = 600,
     units = "in")
ggplot() +
    geom_line(mapping = aes(rt,
                            intensity,
                            colour = group,
                            group = sample),
              data = chr_kappa.df,
              lwd = 1) +
    scale_color_npg(name = "") +
    labs(x= "Retention time (min)",
         y = "Intensity (a.u.)") +
    facet_grid(rows = vars(ion), scales = "free_y", labeller = ion_labeller) +
    xlim(0, 30) +
    scale_y_continuous(breaks = breaks_extended(n=3),
                       labels = scientific_function) +
    theme_classic() +
    theme(strip.text.y = element_text(size = 8, 
                                      family = "Avenir LT 65 Medium",
                                      angle = 360,
                                      hjust = 0),
          strip.background = element_blank(),
          axis.text = element_text(size = 12,
                                   family = "Avenir"),
          axis.title = element_text(size = 12, 
                                    family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "#848587",
                                      size = 0.5,
                                      fill = NA),
          legend.position = "bottom",
          legend.text = element_text(size = 10, 
                                     family = "Avenir"),
          axis.line = element_blank())
dev.off()

predicted <- fread("predicted_sugars.txt")
data_digest <- filterFile(data,
                          file = which(pd$sample_type == "digest"))


x <- predicted$`[M-3H]-3`[predicted$name == "hex-2-sulphate-3"]






