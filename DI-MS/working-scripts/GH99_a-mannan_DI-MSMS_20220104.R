setwd("/Users/margotbligh/ownCloud/marglyco/Data/LC-MS_data/alpha-mannan/20211221_Direct_injection")

library(xcms) ; library(data.table) ; library(scales) ; library(tidyverse)
library(ggplot2)


data <- readMSData(files = "20211221_MS31_alphamannan_DI2_GH99_1to100dil.mzML", 
                   mode = "onDisk")
data.ms2 <- data[data@featureData@data$msLevel == 2]

precursor_mz <- 301.02349557990897
data.ms2.301 <- data.ms2[data.ms2@featureData@data$precursorMZ > 300 &
                             data.ms2@featureData@data$precursorMZ < 302]

ms2.301.mean <- combineSpectra(data.ms2.301, fcol = "fileIdx", mzFun = base::mean)

plot(mz(ms2.301.mean[[1]]), intensity(ms2.301.mean[[1]]), type = "h", col = "light blue")

df1 <- data.frame(mz = mz(ms2.301.mean[[1]]), intensity = intensity(ms2.301.mean[[1]]))
df1$norm_intensity <- df1$intensity/max(df1$intensity)

df2 <- df1 %>% filter(norm_intensity > 0.01)
df2$mz_round <- round(df2$mz, 4)

p <- ggplot(data = df2) + 
    geom_segment(aes(x=mz_round, xend=mz_round,  y=0, yend=norm_intensity),
                 lwd = 1.3) +
    #precursor arrow
    geom_point(aes(x = 301.0229,y = 0.02),shape = 25, size = 3,
               fill = "#FEC000", colour = "black") +
    #fragment ion m/z values
    geom_text(aes(x = mz_round, y = norm_intensity, 
                  label = as.character(mz_round)),
              angle = 90, vjust = 0.5, nudge_y =0.05,
              hjust = 0)+ 
    ylim(0,1.5) +
    labs(x = expression(italic(m/z)), y = "Normalised intensity (a.u.)")+
    theme_classic() + 
    theme(axis.text = element_text(size = 10,family = "Avenir"),
          axis.title = element_text(size = 10,family = "Avenir LT 65 Medium"),
          panel.border = element_rect(colour = "black",size = 0.5,fill = NA),
          plot.title = element_text(size = 10,family = "Avenir LT 65 Medium",
                                    hjust = 0.5),
          axis.line = element_blank())

png(filename = "301-ms2_plot_v2.png", width = 7, height = 4, units = "in",
    res = 400)
print(p)
dev.off()


svg(filename = "301-ms2_plot_v2.svg", width = 7, height = 4)
print(p)
dev.off()

lc.data <- readMSData(files = "../MS78_20211216_GH99-alpha-mannan_25.mzML", 
                   mode = "onDisk")

chr_301 <- chromatogram(lc.data, mz = c(301.02349557990897-0.001,
                                        301.02349557990897+0.001))

chr_301.df <- data.frame(rt = chr_301[,1]@rtime/60,
                         intensity = chr_301[,1]@intensity)

scientific_function <- function(x) {
    text <- gsub("E0", "", gsub("e\\+0", "E", scales::scientific_format()(x)))
    text
}

p <- ggplot(data = chr_301.df) +
    geom_line(aes(x = rt, y = intensity)) +
    xlim(0,25) +
    scale_y_continuous(labels = scientific_function) +
    labs(x = "Retention time (min)", y = "Intensity (a.u.)") +
    theme_classic() + 
    theme(text = element_text(family = "Avenir"),
          panel.border = element_rect(colour = "#848587", size = 0.5, fill = NA),
          axis.line = element_blank())
png(filename = "301_chromatogram_v1.png",
    height = 3, width = 8, units = "in", res = 300)
print(p)
dev.off()

chr_247 <- chromatogram(lc.data, mz = c(247.0057-0.001,
                                        247.0057+0.001))

plot(chr_301)


