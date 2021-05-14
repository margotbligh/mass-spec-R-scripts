# 1: Install packages --------------------------------------------------------

library(xcms)
library(MSnbase)
library(mzR)
library(tidyverse)
library(data.table)
library(ggplot2)
library(bit64)
library(RColorBrewer)

#2. Import MS1 data --------------------------------------------------------



setwd("~/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_FITDOG/2a_neutral_standards/MS-direct-injection_20201111")
#get file paths to mzML files
raw_files_path <- dir(path = "./mzML-files", all.files = FALSE, full.names = TRUE)

#create phenodata data.frame
sample = basename(raw_files_path) %>% 
    sub(".mzML", "", .) %>% 
    sub("20201111_neg_DI\\d+_", "", .) %>% 
    sub("_100000xdilute", "", .) %>% 
    sub("-\\d", "", .) %>% 
    sub("-digested|-undigested", "", .) %>% 
    sub("lam", "laminarin", .) %>% 
    sub("yeastman", "yeast-mannan", .)

treatment = basename(raw_files_path) %>% 
    sub(".mzML", "", .) %>% 
    sub("20201111_neg_DI\\d+_", "", .) %>% 
    sub("_100000xdilute", "", .) %>% 
    sub("-\\d", "", .) %>% 
    sub("lam|yeastman", "", .) %>% 
    gsub("methanol|reaction|-blank|-", "", .)

replicate = basename(raw_files_path) %>% 
    sub(".mzML", "", .) %>% 
    sub("20201111_neg_DI\\d+_", "", .) %>% 
    sub("_100000xdilute", "", .) %>% 
    gsub("lam|yeastman|methanol|reaction|-digested|-undigested|-blank|-", "", .)


pd <- data.frame(name = sub(basename(raw_files_path), 
                                   pattern = ".mzML", 
                                   replacement = ""), 
                 sample = sample, 
                 treatment = treatment,
                 replicate = replicate,
                 stringsAsFactors = FALSE)

pd <- apply(pd, 2, function(x) gsub("^$|^ $", NA, x)) #replace empty cells with NA
pd <- as.data.frame(pd)

pd$short_name <- paste0(pd$sample, "-", 
                        pd$treatment, "-", 
                        pd$replicate) %>% 
    gsub("-NA", "", .)

#load the raw data
raw_data <- readMSData(files = raw_files_path, 
                       pdata = new("NAnnotatedDataFrame", 
                                   pd), 
                       mode = "onDisk")

#2. Plot chromatograms for specific ions ---------------------------
#import mz values
mz_values.df <- fread("~/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/2_FITDOG/2a_neutral_standards/masses_predicted/masses_predicted_dp2-13_neg_withoutpent_scan150-2500.txt")

#error
error = 0.001

#create output directory
dir.create(file.path("analysis/plots"), showWarnings = FALSE)

# colours
sample_colours <- paste0(wes_palette(n = 4, name = "GrandBudapest2"))
names(sample_colours) <- unique(pd$sample)

sample_name_colours <- brewer.pal(9, "Set3")
names(sample_name_colours) <- pd$short_name

# extract chromatogram and plot for each dp
dp = mz_values.df$dp

for (i in 1:length(dp)){
    mz <- mz_values.df %>% 
        filter(dp == dp[i]) %>% 
        select('[M-H]') %>% 
        as.double()
    mzr <- c(mz - error, 
             mz + error)
    chr_raw <- chromatogram(raw_data, 
                            mz = mzr)
    
    cairo_pdf(paste0("analysis/plots/chromatogram_DP",
                     dp[i],
                     ".pdf"),
              family = "Avenir",
              pointsize = 25,
              width = 12,
              height = 9)
    
    par(mar=c(5,4,2,7))
    
    plot(chr_raw, 
         col = sample_name_colours[chr_raw$short_name],
         main = "",
         xlab = "",
         ylab = "")
    
    title(main = paste0("DP",
                        dp[i],
                        " [M-H] chromatrogam: mz = ",
                        mz,
                        " +- 0.001"),
          cex.main = 0.8,
          font.main = 2,
          xlab = "injection time (s)",
          ylab = "signal intensity",
          cex.lab = 0.8,
          font.lab = 2)
         
    
    legend("topleft",
           legend = pd$short_name,
           fill = sample_name_colours,
           pt.cex = 0.3,
           cex = 0.5,
           bty = "n",
           ncol = 1,
           x.intersp = 0.2,
           inset=c(1,0),
           xpd=TRUE)
    
    dev.off()
    
    
}



  