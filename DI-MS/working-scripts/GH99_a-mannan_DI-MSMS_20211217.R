setwd("/Users/margotbligh/ownCloud/marglyco/Data/LC-MS_data/alpha-mannan/20211215_Direct_injection")

library(xcms) ; library(data.table) ; library(scales) ; library(tidyverse)
library(ggplot2)


data <- readMSData(files = "20211216_MS31_alphamannan_DI2_GH99_1to100dil.mzML", 
                   mode = "onDisk")
data.ms2 <- data[data@featureData@data$msLevel == 2]
data.ms2.rt <- data.ms2 %>% filterRt(rt = c(330,400))
table(fromFile(data))
table(fromFile(data.ms2))
table(fromFile(data.ms2.rt))


ms2.mean <- combineSpectra(data.ms2.rt, fcol = "fileIdx", mzFun = base::mean)

plot(mz(ms2.mean[[1]]), intensity(ms2.mean[[1]]), type = "h", col = "light blue")

df1 <- data.frame(mz = mz(ms2.mean[[1]]), intensity = intensity(ms2.mean[[1]]))
df1$norm_intensity <- df1$intensity/max(df1$intensity)

df2 <- df1 %>% filter(norm_intensity > 0.1)
df2$mz_round <- round(df2$mz, 4)

df2$annotation[df2$mz_round == 247.049] <- "tri-sulphated mannotriose [M-3H]-3"
df2$annotation[df2$mz_round == 97.0457] <- "HSO4-"
df2$annotation[df2$mz_round == 97.0457] <- "HSO4-"

p <- ggplot(data = df2) + 
    geom_segment(aes(x=mz_round, xend=mz_round,  y=0, yend=norm_intensity),
                 lwd = 1.3) +
    #precursor arrow
    geom_point(aes(x = 247.0490,y = 0.02),shape = 25, size = 3,
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

svg(filename = "ms2_plot_v2.svg", width = 7, height = 4)
print(p)
dev.off()


