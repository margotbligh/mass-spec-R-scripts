#setwd("/Users/margotbligh/ownCloud/marglyco/Data/LC-MS_data/alpha-mannan")
setwd("/Users/mbligh/ownCloud/marglyco/Data/LC-MS_data/alpha-mannan")

load("20220211_RData.RData")

#change directory names for new Macbook
x <- dirname(data.di.ms2) 
x <- x %>% sub("margotbligh", "mbligh", .)
dirname(data.di.ms2) <- x

#palettes
cold <- hcl.colors(n = 5, palette = "Cold")
warm <- hcl.colors(n = 5, palette = "Warm")
pretty3 <- c(warm[1], cold[1:2])
#which runs do I want to plot?
#03.02 - full digests
#08.02 - full digest June, tSIM
#09.02 - fractions 10-11, 13-18, full digest June tSIM
#10.02 - fractions 10-11, 13-18, full digest June (tSIM, fully deprotonated); fractions 19-22 and 23-26 tSIM

#Load libraries
library(tidyr) ; library(xcms) ; library(data.table) ; library(ggplot2)
library(scales) ; library(ggsci) ; library(ggridges) ; library(ggpubr)
library(reticulate) ; library(ggrepel) ; library(MSnbase)

#1: Load LC-MS data----
#get filepaths
fp <- dir(".", pattern = ".mzML", recursive = T)
fp <- fp[!grepl("Direct|solvent|20220207|20211216", fp)]
fp <- fp[c(5,6,9:17)]

#create phenodata frame
pd <- data.frame(name = basename(fp) %>% 
                     sub(".*MS31_", "", .) %>% 
                     sub(".mzML", "", .),
                 sampletype = basename(fp) %>% 
                     sub(".*_fractions.*", "pooled fractions", .) %>% 
                     sub(".*full.*", "full digest", .),
                 fractions = basename(fp) %>% 
                     sub(".*10-11.*", "10-11", .) %>% 
                     sub(".*13-18.*", "13-18", .) %>% 
                     sub(".*19-22.*", "19-22", .) %>% 
                     sub(".*23-26.*", "23-26", .) %>% 
                     sub(".*full.*", "", .),
                 date = basename(fp) %>% 
                     sub("MS31_", "", .) %>% 
                     sub("_.*", "", .))
pd$runtype <- "tSIM"
pd$runtype[pd$date == "20220203"] <- "full MS"
pd$targetions <- "most abundant"
pd$targetions[pd$date == "20220210" & pd$fractions != "19-22" &
                  pd$fractions != "23-26"] <- "fully deprotonated"
pd$targetions[pd$runtype == "full MS"] <- ""

#load data
data <- readMSData(files = fp, pdata = new("NAnnotatedDataFrame", pd), 
                   mode = "onDisk")

data.ms1 <- data[data@featureData@data$msLevel == 1]

#2: Extract chromatograms-------
mz.for.chr <- c(259.0131, 421.0658, 250.0078, 377.0853, 583.1192, 331.0337,
                452.039, 573.0444, 694.0494, 543.0338, 301.023, 286.0182,
                277.015, 271.0129)
names(mz.for.chr) <- c("hex-1-sulphate-1: [M-H]-", "hex-2-sulphate-1: [M-H]-",
                       "hex-2-sulphate-2: [M-2H]-2", "hex-2: [M+Cl]-",
                       "hex-3-sulphate-1: [M-H]-", "hex-3-sulphate-2: [M-2H]-2",
                       "hex-4-sulphate-3: [M-2H]-2", "hex-5-sulphate-4: [M-2H]-2",
                       "hex-6-sulphate-5: [M-2H]-2", "hex-7-sulphate-6: [M-3H]-3",
                       "hex-4-sulphate-3: [M-3H]-3", "hex-5-sulphate-4: [M-4H]-4",
                       "hex-6-sulphate-5: [M-5H]-5", "hex-7-sulphate-6: [M-6H]-6")

chr.list <- list()
for(i in 1:length(mz.for.chr)){
    chr.list[[i]] <- chromatogram(data.ms1, mz = c(mz.for.chr[i] - 0.001,
                                                   mz.for.chr[i] + 0.001))
}

#3: Create dataframe with chromatograms-----

#build dataframe for plotting with ggplot
chr.df <- data.frame(name = as.character(),
                     sampletype = as.character(),
                     fractions = as.character(),
                     date = as.character(), runtype = as.character(),
                     targetions = as.character(),
                     mz = as.numeric(), ion = as.character(),
                     rt = as.numeric(), intensity = as.numeric())

for(i in 1:length(mz.for.chr)){
    chr <- chr.list[[i]]
    mz <- mz.for.chr[i]
    ion <- names(mz.for.chr)[i]
    
    lengths <- unlist(lapply(chr, function(x) length(x)))
    
    name <- rep(chr$name, times = lengths)
    sampletype <- rep(chr$sampletype, times = lengths)
    fractions <- rep(chr$fractions, times = lengths)
    date <- rep(chr$date, times = lengths)
    runtype <- rep(chr$runtype, times = lengths)
    targetions <- rep(chr$targetions, times = lengths)
    mz <- rep(mz, sum(lengths))
    ion <- rep(ion, sum(lengths))
    rt <- unlist(lapply(chr, function(x) x@rtime/60))
    intensity <- unlist(lapply(chr, function(x) x@intensity))
    
    
    temp <- data.frame(name = name, sampletype = sampletype,
                       fractions = fractions, date = date,
                       runtype = runtype, targetions = targetions,
                       mz = mz, ion = ion,
                       rt = rt, intensity = intensity)
    
    chr.df <- rbind(chr.df, temp)
}


#4: Chromatogram plots ----
scientific_function <- function(x){
    text <- gsub("E0", "", gsub("e\\+0", "E", scales::scientific_format()(x)))
    text
}

chr.df$group <- paste0(chr.df$sampletype, ": ", chr.df$fractions) %>% 
    sub(":\\s$", "", .) 

#basic chromatograms in tSIM mode
chr.df.tSIM.abundant <- chr.df %>% 
    dplyr::filter(runtype == "tSIM") %>% 
    dplyr::filter(targetions == "most abundant")

chr.df.tSIM.abundant.noNA <- chr.df.tSIM.abundant[complete.cases(chr.df.tSIM.abundant),]
chr.df.tSIM.abundant.noNA2 <- chr.df.tSIM.abundant.noNA %>% 
    dplyr::filter(name != "20220208_fulldigest_june_100xdilute_05")

chr.df.tSIM.abundant.noNA2$ion <- factor(chr.df.tSIM.abundant.noNA2$ion,
                                         levels = c("hex-1-sulphate-1: [M-H]-",
                                                    "hex-2: [M+Cl]-",  
                                                    "hex-2-sulphate-1: [M-H]-",
                                                    "hex-2-sulphate-2: [M-2H]-2", 
                                                    "hex-3-sulphate-1: [M-H]-",
                                                    "hex-3-sulphate-2: [M-2H]-2",
                                                    "hex-4-sulphate-3: [M-2H]-2", 
                                                    "hex-5-sulphate-4: [M-2H]-2",
                                                    "hex-6-sulphate-5: [M-2H]-2",
                                                    "hex-7-sulphate-6: [M-3H]-3"),
                                         labels = c("sulphated mannose: [M-H]-",
                                                    "mannobiose: [M+Cl]-",  
                                                    "sulphated mannobiose: [M-H]-",
                                                    "disulphated mannobiose: [M-2H]-2", 
                                                    "sulphated mannotriose: [M-H]-",
                                                    "disulphated mannotriose: [M-2H]-2",
                                                    "trisulphated mannotetraose: [M-2H]-2", 
                                                    "tetrasulphated mannopentaose: [M-2H]-2",
                                                    "pentasulphated mannohexaose: [M-2H]-2",
                                                    "hexasulphated mannoheptaose: [M-3H]-3"))

p1 <- ggplot(data = chr.df.tSIM.abundant.noNA2 %>% 
                 dplyr::filter(ion != "disulphated mannobiose: [M-2H]-2"), 
       aes(x = rt, y = sqrt(intensity), height = sqrt(intensity), group = ion, 
           colour = ion)) +
    geom_ridgeline(lwd = 0.8, fill = NA) +
    scale_colour_npg(name = "") +
    facet_grid(rows = vars(group), scales = "free") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function) +
    labs(x = "Retention time (min)", y= expression(bold(sqrt(Intensity~"(a.u.)")))) +
    scale_x_continuous(breaks = seq(0,25, by = 2)) +
    theme_classic() +
    theme(text = element_text(family = "Arial"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 12, face = "bold"),
          axis.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 12),
          legend.position = "top",
          strip.background = element_rect(colour = NA)) +
    guides(colour=guide_legend(nrow=4,byrow=TRUE))

png("GH99_plots_202202/chromatograms_all_tSIM_v1.png",
    res = 450, width = 12, height = 7, units = "in")
p1
dev.off()


#chromatogram in tSIM mode for fully deprotonated ions
chr.df.tSIM.deprot <- chr.df %>% 
    dplyr::filter(runtype == "tSIM") %>% 
    dplyr::filter(targetions == "fully deprotonated")

chr.df.tSIM.deprot.noNA <- chr.df.tSIM.deprot[complete.cases(chr.df.tSIM.deprot),]

chr.df.tSIM.deprot.noNA$ion <- factor(chr.df.tSIM.deprot.noNA$ion,
                                         levels = c("hex-1-sulphate-1: [M-H]-",
                                                    "hex-2: [M+Cl]-",  
                                                    "hex-2-sulphate-1: [M-H]-",
                                                    "hex-2-sulphate-2: [M-2H]-2", 
                                                    "hex-3-sulphate-1: [M-H]-",
                                                    "hex-3-sulphate-2: [M-2H]-2",
                                                    "hex-4-sulphate-3: [M-3H]-3", 
                                                    "hex-5-sulphate-4: [M-4H]-4",
                                                    "hex-6-sulphate-5: [M-5H]-5",
                                                    "hex-7-sulphate-6: [M-6H]-6"),
                                         labels = c("sulphated mannose: [M-H]-",
                                                    "mannobiose: [M+Cl]-",  
                                                    "sulphated mannobiose: [M-H]-",
                                                    "disulphated mannobiose: [M-2H]-2", 
                                                    "sulphated mannotriose: [M-H]-",
                                                    "disulphated mannotriose: [M-2H]-2",
                                                    "trisulphated mannotetraose: [M-3H]-3", 
                                                    "tetrasulphated mannopentaose: [M-4H]-4",
                                                    "pentasulphated mannohexaose: [M-5H]-5",
                                                    "hexasulphated mannoheptaose: [M-6H]-6"))

p2 <- ggplot(data = chr.df.tSIM.deprot.noNA %>% 
           dplyr::filter(ion != "disulphated mannobiose: [M-2H]-2") %>% 
           dplyr::filter(name == "20220210_fulldigest_june_100xdilute_04"), 
             aes(x = rt, y = sqrt(intensity), height = sqrt(intensity), group = ion, 
                 colour = ion)) +
    geom_ridgeline(lwd = 0.8, fill = NA) +
    scale_colour_npg(name = "") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function) +
    labs(x = "Retention time (min)", y= expression(bold(sqrt(Intensity~"(a.u.)")))) +
    scale_x_continuous(breaks = seq(0,25, by = 2)) +
    theme_classic() +
    theme(text = element_text(family = "Arial"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 12, face = "bold"),
          axis.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 12),
          legend.position = "top",
          strip.background = element_rect(colour = NA)) +
    guides(colour=guide_legend(nrow=4,byrow=TRUE))
p3 <- ggplot(data = chr.df.tSIM.abundant.noNA2 %>% 
                 dplyr::filter(ion != "disulphated mannobiose: [M-2H]-2") %>% 
                 dplyr::filter(group == "full digest"), 
             aes(x = rt, y = sqrt(intensity), height = sqrt(intensity), group = ion, 
                 colour = ion)) +
    geom_ridgeline(lwd = 0.8, fill = NA) +
    scale_colour_npg(name = "") +
    scale_y_continuous(breaks = breaks_pretty(3), labels = scientific_function) +
    labs(x = "Retention time (min)", y= expression(bold(sqrt(Intensity~"(a.u.)")))) +
    scale_x_continuous(breaks = seq(0,25, by = 2)) +
    theme_classic() +
    theme(text = element_text(family = "Arial"),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 12, face = "bold"),
          axis.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 12),
          legend.position = "top",
          strip.background = element_rect(colour = NA)) +
    guides(colour=guide_legend(nrow=4,byrow=TRUE))

png("GH99_plots_202202/chromatograms_tSIM_full_tSIM.png",
    res = 450, width = 12, height = 6, units = "in")
ggarrange(p3, p2, nrow =2)
dev.off()

#5: Load direct injection MS/MS data----
fp.di <- dir(path = "./20220210_Direct_injection", pattern = ".mzML",
             full.names = T)
fp.di <- fp.di[grep("MSMS", fp.di)]

data.di <- readMSData(files = fp.di, mode = "onDisk")
data.di.ms2 <- data.di[data.di@featureData@data$msLevel == 2]

#average data
data.di.ms2.mean <- combineSpectra(data.di.ms2, fcol = "fileIdx", 
                                   intensityFun = max, mzd = 0.005)

#6: Annotate and normalise ms2 spectra-----
#predict sugars
predictSugars <- function(dp1,dp2,ESI_mode, scan_range1, scan_range2,
                          pent_option=NULL,modifications=NULL,
                          nmod_max=NULL,double_sulphate=NULL,label=NULL){
    py_install("pandas", "numpy")
    source_python("../sugarMassesPredict-r.py")
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

predicted_sugars <- predictSugars(dp1 = 1, dp2 = 7, ESI_mode = 'neg', 
                                  scan_range1 = 50, scan_range2 = 1000,
                                  pent_option = 0, modifications = c("sulphate"),
                                  double_sulphate = 1, nmod_max = 2)

predicted <- predicted_sugars %>% 
    pivot_longer(cols = starts_with("[M"), values_to = "mz",
                 names_to = "ion")
predicted <- predicted[complete.cases(predicted),]
predicted <- predicted[!grepl("CHOO|Cl", predicted$ion),]

#add triple sulphated hexose and h2so4
predicted <- rbind(predicted,
                   c(1, "hex-1-sulphate-3", 419.933841, "C6H12O15S3", 
                     "[M-3H]-3", 138.9701))
predicted <- rbind(predicted,
                   c(0, "HSO4", 96.959557, "HSO4", 
                     "", 96.959557))

#extra spectra as dataframe
msms.df <- data.frame(precursorMz = as.numeric(),
                      mz = as.numeric(),
                      intensity = as.numeric())
for(i in 1:length(fp.di)){
    mz <- data.di.ms2.mean[[i]]@mz
    intensity <- data.di.ms2.mean[[i]]@intensity
    precursorMz <- rep(data.di.ms2.mean[[i]]@precursorMz, length(mz))
    tmp <- data.frame(precursorMz = precursorMz,mz =mz,intensity = intensity)
    msms.df <- rbind(msms.df, tmp)
}

msms.df$precursorMz <- round(msms.df$precursorMz, 2)

msms.df$precursorIon[msms.df$precursorMz == 250.00] <- "disulphated mannobiose: [M-2H]-2"
msms.df$precursorIon[msms.df$precursorMz == 271.00] <- "hexasulphated mannoheptaose: [M-6H]-6"
msms.df$precursorIon[msms.df$precursorMz == 277.00] <- "pentasulphated mannohexaose: [M-5H]-5"
msms.df$precursorIon[msms.df$precursorMz == 286.00] <- "tetrasulphated mannopentaose: [M-4H]-4"
msms.df$precursorIon[msms.df$precursorMz == 301.00] <- "trisulphated mannotetraose: [M-3H]-3"
msms.df$precursorIon[msms.df$precursorMz == 331.05] <- "disulphated mannotriose: [M-2H]-2"
msms.df$precursorIon[msms.df$precursorMz == 452.05] <- "trisulphated mannotetraose: [M-2H]-2"
msms.df$precursorIon[msms.df$precursorMz == 573.05] <- "tetrasulphated mannopentaose: [M-2H]-2"

#overlap annotation
predicted$mzmin <- as.numeric(predicted$mz) - 0.001
predicted$mzmax <- as.numeric(predicted$mz) + 0.001

msms.df$mzmin <- as.numeric(msms.df$mz) - 0.001
msms.df$mzmax <- as.numeric(msms.df$mz) + 0.001

setDT(predicted); setDT(msms.df)
setkey(predicted, mzmin, mzmax)
msms.df_nopred <- msms.df
msms.df <- foverlaps(msms.df_nopred,predicted)

#change NA values created during matching (features with no match) to be blank
#remove extra columns
msms.df$mz <- NULL
msms.df$mzmin <- NULL
msms.df$mzmax <- NULL

names(msms.df) <- names(msms.df) %>% 
    sub("^i\\.", "", .)

msms.df$mzmin <- NULL
msms.df$mzmax <- NULL

msms.df <- msms.df %>% 
    replace_na(list("dp"="","name"="", "mass"= "", "formula" = "", "ion"= ""))

#normalise intensities
msms.df <- msms.df %>% 
    dplyr::group_by(precursorIon) %>% 
    dplyr::mutate(normalised_intensity = intensity/max(intensity))

#7: Plots-----
    #disulphated mannotriose: [M-2H]-2 ----
h3s2.df <- msms.df %>% 
    dplyr::filter(precursorIon == "disulphated mannotriose: [M-2H]-2") %>% 
    dplyr::filter(normalised_intensity > 0.005)

h3s2.df <- h3s2.df[h3s2.df$name != "hex-6-sulphate-11",]

h3s2.df$label <- paste0(h3s2.df$name, "\n", h3s2.df$ion) %>% 
    sub("\\\n$", "", .)
h3s2.df$label[round(h3s2.df$mz,2)== 331.06] <- "disulphated mannotriose\n[M-2H]-2"
h3s2.df$label[round(h3s2.df$mz,2)== 96.96] <- "HSO4"

png(filename = "GH99_plots_202202/disulphated_mannotriose-ms2_plot_v2.png", 
    width = 6, height = 3, units = "in", res = 400)
ggplot(data = h3s2.df) + 
    geom_segment(aes(x=mz, xend=mz,  y=0, yend=normalised_intensity),
                 lwd = 1) +
    #precursor arrow
    geom_point(aes(x = precursorMz, y = 0.02),shape = 25, size = 3,
               fill = "#FEC000", colour = "black") +
    #labels of fragments
    coord_cartesian(clip = "off") +
    geom_label_repel(data = h3s2.df[h3s2.df$label != "",],
                     aes(x = mz, y = normalised_intensity, label = label),
                     min.segment.length = 0, box.padding = 0.5, nudge_y = 0.1,
                     max.overlaps = Inf,  xlim = c(0, 450), ylim = c(0.1, Inf),
                     size = 3, segment.linetype = "dashed", 
                     fill = alpha(c("white"),0.7),segment.colour = "#2BB6AF") +
    ylim(0, 1.5) +
    labs(x = expression(italic(m/z)), y = "Normalised intensity (a.u.)",
         tag = "A")+
    scale_x_continuous(expand = expansion(c(0.3,0.2))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10,family = "Avenir"),
          axis.title = element_text(size = 10,family = "Avenir LT 65 Medium"),
          #panel.border = element_rect(colour = "black",size = 0.5,fill = NA),
          plot.title = element_text(size = 10,family = "Avenir LT 65 Medium",
                                    hjust = 0.5))
dev.off()


    #trisulphated mannotetraose: [M-3H]-3-----
h4s3.df <- msms.df %>% 
    dplyr::filter(precursorIon == "trisulphated mannotetraose: [M-3H]-3") %>% 
    dplyr::filter(normalised_intensity > 0.001)

h4s3.df <- h4s3.df[h4s3.df$name != "hex-6-sulphate-11",]

h4s3.df$label <- paste0(h4s3.df$name, "\n", h4s3.df$ion) %>% 
    sub("\\\n$", "", .)
h4s3.df$label[round(h4s3.df$mz,2)== 286.01] <- "trisulphated mannotetraose\n[M-3H]-3"
h4s3.df$label[round(h4s3.df$mz,2)== 96.96] <- "HSO4"

png(filename = "GH99_plots_202202/trisulphated_mannotetraose-ms2_plot_v2.png", 
    width = 6, height = 3, units = "in", res = 400)
ggplot(data = h4s3.df) + 
    geom_segment(aes(x=mz, xend=mz,  y=0, yend=normalised_intensity),
                 lwd = 1) +
    #precursor arrow
    geom_point(aes(x = precursorMz, y = 0.02),shape = 25, size = 3,
               fill = "#FEC000", colour = "black") +
    #labels of fragments
    coord_cartesian(clip = "off") +
    geom_label_repel(data = h4s3.df[h4s3.df$label != "",],
                     aes(x = mz, y = normalised_intensity, label = label),
                     min.segment.length = 0, box.padding = 0.5, nudge_y = 0.1,
                     max.overlaps = Inf,  xlim = c(0, 450), ylim = c(0.1, Inf),
                     size = 3, segment.linetype = "dashed", 
                     fill = alpha(c("white"),0.7),segment.colour = "#2BB6AF") +
    ylim(0, 1.5) +
    labs(x = expression(italic(m/z)), y = "Normalised intensity (a.u.)",
         tag = "B")+
    scale_x_continuous(expand = expansion(c(0.3,0.2))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10,family = "Avenir"),
          axis.title = element_text(size = 10,family = "Avenir LT 65 Medium"),
          #panel.border = element_rect(colour = "black",size = 0.5,fill = NA),
          plot.title = element_text(size = 10,family = "Avenir LT 65 Medium",
                                    hjust = 0.5))
dev.off()

    #tetrasulphated mannopentaose: [M-4H]-4-----
h5s4.df <- msms.df %>% 
    dplyr::filter(precursorIon == "tetrasulphated mannopentaose: [M-4H]-4") %>% 
    dplyr::filter(normalised_intensity > 0.005)
h5s4.df <- h5s4.df[h5s4.df$name != "hex-6-sulphate-11",]

h5s4.df$label <- paste0(h5s4.df$name, "\n", h5s4.df$ion) %>% 
    sub("\\\n$", "", .)
h5s4.df$label[round(h5s4.df$mz,2)== 286.01] <- "tetrasulphated mannopentaose\n[M-4H]-4"

png(filename = "GH99_plots_202202/tetrasulphated_mannopentaose-ms2_plot_v2.png", 
    width = 6, height = 3, units = "in", res = 400)
ggplot(data = h5s4.df) + 
    geom_segment(aes(x=mz, xend=mz,  y=0, yend=normalised_intensity),
                 lwd = 1) +
    #precursor arrow
    geom_point(aes(x = precursorMz, y = 0.02),shape = 25, size = 3,
               fill = "#FEC000", colour = "black") +
    #labels of fragments
    coord_cartesian(clip = "off") +
    geom_label_repel(data = h5s4.df[h5s4.df$label != "",],
                     aes(x = mz, y = normalised_intensity, label = label),
                     min.segment.length = 0, box.padding = 0.5, nudge_y = 0.1,
                     max.overlaps = Inf,  xlim = c(0, 450), ylim = c(0.1, Inf),
                     size = 3, segment.linetype = "dashed", 
                     fill = alpha(c("white"),0.7),segment.colour = "#2BB6AF") +
    ylim(0, 1.5) +
    labs(tag = "C") +
    labs(x = expression(italic(m/z)), y = "Normalised intensity (a.u.)")+
    scale_x_continuous(expand = expansion(c(0.3,0.2))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10,family = "Avenir"),
          axis.title = element_text(size = 10,family = "Avenir LT 65 Medium"),
          #panel.border = element_rect(colour = "black",size = 0.5,fill = NA),
          plot.title = element_text(size = 10,family = "Avenir LT 65 Medium",
                                    hjust = 0.5))
dev.off()


    #pentasulphated mannohexaose: [M-5H]-5-----
h6s5.df <- msms.df %>% 
    dplyr::filter(precursorIon == "pentasulphated mannohexaose: [M-5H]-5") %>% 
    dplyr::filter(normalised_intensity > 0.005)
h6s5.df <- h6s5.df[h6s5.df$name != "hex-6-sulphate-11",]

h6s5.df$label <- paste0(h6s5.df$name, "\n", h6s5.df$ion) %>% 
    sub("\\\n$", "", .)
h6s5.df$label[round(h6s5.df$mz,2)== 277.08] <- "pentasulphated mannohexaose\n[M-5H]-5"

png(filename = "GH99_plots_202202/pentasulphated_mannohexaose-ms2_plot_v2.png", 
    width = 6, height = 3, units = "in", res = 400)
ggplot(data = h6s5.df) + 
    geom_segment(aes(x=mz, xend=mz,  y=0, yend=normalised_intensity),
                 lwd = 1) +
    #precursor arrow
    geom_point(aes(x = precursorMz, y = 0.02),shape = 25, size = 3,
               fill = "#FEC000", colour = "black") +
    #labels of fragments
    coord_cartesian(clip = "off") +
    geom_label_repel(data = h6s5.df[h6s5.df$label != "",],
                     aes(x = mz, y = normalised_intensity, label = label),
                     min.segment.length = 0, box.padding = 0.5, nudge_y = 0.1,
                     max.overlaps = Inf,  xlim = c(0, 450), ylim = c(0.1, Inf),
                     size = 3, segment.linetype = "dashed", 
                     fill = alpha(c("white"),0.7),segment.colour = "#2BB6AF") +
    ylim(0, 1.5) +
    labs(tag = "D") +
    labs(x = expression(italic(m/z)), y = "Normalised intensity (a.u.)")+
    scale_x_continuous(expand = expansion(c(0.3,0.2))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10,family = "Avenir"),
          axis.title = element_text(size = 10,family = "Avenir LT 65 Medium"),
          #panel.border = element_rect(colour = "black",size = 0.5,fill = NA),
          plot.title = element_text(size = 10,family = "Avenir LT 65 Medium",
                                    hjust = 0.5))
dev.off()


#8: MSMS plots for mass spec group meeting----
    #hexasulphated mannoheptaose: [M-6H]-6-----
h7s6.df <- msms.df %>% 
    dplyr::filter(precursorIon == "hexasulphated mannoheptaose: [M-6H]-6") %>% 
    dplyr::filter(normalised_intensity > 0.005)
h7s6.df <- h7s6.df[h7s6.df$name != "hex-6-sulphate-11",]

h7s6.df$label <- paste0(h7s6.df$name, "\n", h7s6.df$ion) %>% 
    sub("\\\n$", "", .)
h7s6.df$label[round(h7s6.df$mz,2)== 271.08] <- "hex-7-sulphate-6\n[M-6H]-6"

svg(filename = "GH99_plots_202202/hexasulphated_mannoheptaose-ms2_plot_v2.svg", 
    width = 8, height = 6)
cairo_pdf(filename = "GH99_plots_202202/hexasulphated_mannoheptaose-ms2_plot_v2.pdf", 
    width = 8, height = 6)
ggplot(data = h7s6.df) + 
    geom_segment(aes(x=mz, xend=mz,  y=0, yend=normalised_intensity),
                 lwd = 1) +
    #precursor arrow
    geom_point(aes(x = precursorMz, y = 0.02),shape = 25, size = 3,
               fill = "#FEC000", colour = "black") +
    #labels of fragments
    coord_cartesian(clip = "off") +
    geom_label_repel(data = h7s6.df[h7s6.df$label != "",],
                     aes(x = mz, y = normalised_intensity, label = label),
                     min.segment.length = 0, box.padding = 0.5, nudge_y = 0.1,
                     max.overlaps = Inf,  xlim = c(0, 450), ylim = c(0.1, Inf),
                     size = 2.5, segment.linetype = "dashed", 
                     fill = alpha(c("white"),0.7),segment.colour = "#2BB6AF",
                     family = "Avenir") +
    #label mz values
    geom_text(aes(x = mz, y = normalised_intensity, label = round(mz, 3)),
              angle = 90, nudge_y = 0.08, size = 2,
              family = "Avenir") +
    ylim(0, 1.5) +
    labs(x = expression(italic(m/z)), y = "Normalised intensity (a.u.)")+
    scale_x_continuous(expand = expansion(c(0.3,0.2))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10,family = "Avenir"),
          axis.title = element_text(size = 10,family = "Avenir LT 65 Medium"),
          #panel.border = element_rect(colour = "black",size = 0.5,fill = NA),
          plot.title = element_text(size = 10,family = "Avenir LT 65 Medium",
                                    hjust = 0.5))
dev.off()


    #pentasulphated mannohexaose: [M-5H]-5-----
cairo_pdf(filename = "GH99_plots_202202/pentasulphated_mannohexaose-ms2_plot_v2.pdf", 
          width = 8, height = 6)
ggplot(data = h6s5.df) + 
    geom_segment(aes(x=mz, xend=mz,  y=0, yend=normalised_intensity),
                 lwd = 1) +
    #precursor arrow
    geom_point(aes(x = precursorMz, y = 0.02),shape = 25, size = 3,
               fill = "#FEC000", colour = "black") +
    #labels of fragments
    coord_cartesian(clip = "off") +
    geom_label_repel(data = h6s5.df[h6s5.df$label != "",],
                     aes(x = mz, y = normalised_intensity+0.1, label = label),
                     min.segment.length = 0, box.padding = 0.5, nudge_y = 0.1,
                     max.overlaps = Inf,  xlim = c(0, 450), ylim = c(0.1, Inf),
                     size = 2.5, segment.linetype = "dashed", 
                     fill = alpha(c("white"),0.7),segment.colour = "#2BB6AF",
                     family = "Avenir") +
    #label mz values
    geom_text(aes(x = mz, y = normalised_intensity, label = round(mz, 3)),
              angle = 90, nudge_y = 0.06, size = 2,
              family = "Avenir") +
    ylim(0, 1.5) +
    labs(tag = "D") +
    labs(x = expression(italic(m/z)), y = "Normalised intensity (a.u.)")+
    scale_x_continuous(expand = expansion(c(0.3,0.2))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10,family = "Avenir"),
          axis.title = element_text(size = 10,family = "Avenir LT 65 Medium"),
          #panel.border = element_rect(colour = "black",size = 0.5,fill = NA),
          plot.title = element_text(size = 10,family = "Avenir LT 65 Medium",
                                    hjust = 0.5))
dev.off()



    #disulphated mannotriose: [M-2H]-2-----

h3s2.df$label[round(h3s2.df$mz,2)== 331.06] <- "hex-3-sulphate-2\n[M-2H]-2"
cairo_pdf(filename = "GH99_plots_202202/disulphated_mannotriose-ms2_plot_v2.pdf", 
          width = 8, height = 6)
ggplot(data = h3s2.df) + 
    geom_segment(aes(x=mz, xend=mz,  y=0, yend=normalised_intensity),
                 lwd = 1) +
    #precursor arrow
    geom_point(aes(x = precursorMz, y = 0.02),shape = 25, size = 3,
               fill = "#FEC000", colour = "black") +
    #labels of fragments
    coord_cartesian(clip = "off") +
    geom_label_repel(data = h3s2.df[h3s2.df$label != "",],
                     aes(x = mz, y = normalised_intensity+0.1, label = label),
                     min.segment.length = 0, box.padding = 0.5, nudge_y = 0.1,
                     max.overlaps = Inf,  xlim = c(0, 450), ylim = c(0.1, Inf),
                     size = 2.5, segment.linetype = "dashed", 
                     fill = alpha(c("white"),0.7),segment.colour = "#2BB6AF",
                     family = "Avenir") +
    #label mz values
    geom_text(aes(x = mz, y = normalised_intensity, label = round(mz, 3)),
              angle = 90, nudge_y = 0.06, size = 2,
              family = "Avenir") +
    ylim(0, 1.5) +
    labs(tag = "D") +
    labs(x = expression(italic(m/z)), y = "Normalised intensity (a.u.)")+
    scale_x_continuous(expand = expansion(c(0.3,0.2))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10,family = "Avenir"),
          axis.title = element_text(size = 10,family = "Avenir LT 65 Medium"),
          #panel.border = element_rect(colour = "black",size = 0.5,fill = NA),
          plot.title = element_text(size = 10,family = "Avenir LT 65 Medium",
                                    hjust = 0.5))
dev.off()



#9: Average ion intensities by collision energy-----
#see which columns can be use to combine on
data.di.ms2.f <- fData(data.di.ms2)
#combine on filterString - contains precursor and ionisation enegery
data.di.ms2.mean2 <- combineSpectra(data.di.ms2, fcol = "filterString", 
                                    intensityFun = max, mzd = 0.005)

#10: Extract averaged ms2 data as df----
msms2.df <- data.frame(precursorMz = as.numeric(),
                       nce = as.numeric(),
                       mz = as.numeric(),
                       intensity = as.numeric())
for(i in 1:length(data.di.ms2.mean2)){
    mz <- data.di.ms2.mean2[[i]]@mz
    intensity <- data.di.ms2.mean2[[i]]@intensity
    precursorMz <- rep(data.di.ms2.mean2[[i]]@precursorMz, length(mz))
    nce <- rep(data.di.ms2.mean2[[i]]@collisionEnergy, length(mz))
    
    tmp <- data.frame(precursorMz = precursorMz, nce = nce,
                      mz =mz,intensity = intensity)
    msms2.df <- rbind(msms2.df, tmp)
}

msms2.df$precursorMz <- round(msms2.df$precursorMz, 2)

msms2.df$precursorIon[msms2.df$precursorMz == 250.00] <- "disulphated mannobiose: [M-2H]-2"
msms2.df$precursorIon[msms2.df$precursorMz == 271.00] <- "hexasulphated mannoheptaose: [M-6H]-6"
msms2.df$precursorIon[msms2.df$precursorMz == 277.00] <- "pentasulphated mannohexaose: [M-5H]-5"
msms2.df$precursorIon[msms2.df$precursorMz == 286.00] <- "tetrasulphated mannopentaose: [M-4H]-4"
msms2.df$precursorIon[msms2.df$precursorMz == 301.00] <- "trisulphated mannotetraose: [M-3H]-3"
msms2.df$precursorIon[msms2.df$precursorMz == 331.05] <- "disulphated mannotriose: [M-2H]-2"
msms2.df$precursorIon[msms2.df$precursorMz == 452.05] <- "trisulphated mannotetraose: [M-2H]-2"
msms2.df$precursorIon[msms2.df$precursorMz == 573.05] <- "tetrasulphated mannopentaose: [M-2H]-2"

#11: Annotate ions-----
predicted$mzmin <- as.numeric(predicted$mz) - 0.001
predicted$mzmax <- as.numeric(predicted$mz) + 0.001
names(predicted)[names(predicted)=="name"] <- "sugar"

msms2.df_nopred <- msms2.df
msms2.df$mzmin <- as.numeric(msms2.df$mz) - 0.001
msms2.df$mzmax <- as.numeric(msms2.df$mz) + 0.001

setDT(predicted); setDT(msms2.df)
setkey(predicted, mzmin, mzmax)
msms2.df <- foverlaps(msms2.df,predicted)

#change NA values created during matching (features with no match) to be blank
#remove extra columns
msms2.df$mz <- NULL
msms2.df$mzmin <- NULL
msms2.df$mzmax <- NULL

names(msms2.df) <- names(msms2.df) %>% 
    sub("^i\\.", "", .)

msms2.df$mzmin <- NULL
msms2.df$mzmax <- NULL

msms2.df <- msms2.df %>% 
    replace_na(list("dp"="","sugar"="", "mass"= "", "formula" = "", "ion"= ""))

#12: Normalise intensities----
msms2.df <- msms2.df %>% 
    dplyr::group_by(precursorIon, nce) %>% 
    dplyr::mutate(normalised_intensity = intensity/max(intensity))
msms2.df.filt <- msms2.df %>% dplyr::filter(normalised_intensity > 0.005)

#13: disulphated mannotriose------
h3s2.df2 <- msms2.df %>% 
    dplyr::filter(precursorIon == "disulphated mannotriose: [M-2H]-2")
h3s2.df2$mz_round <- round(h3s2.df2$mz, 3)

#precursor
h3s2.df2$sugar[h3s2.df2$mz_round == 331.073 |
                   h3s2.df2$mz_round == 331.052 |
                   h3s2.df2$mz_round == 331.039] <- "hex-3-sulphate-2"
h3s2.df2$ion[h3s2.df2$mz_round == 331.073 |
                 h3s2.df2$mz_round == 331.052 |
                 h3s2.df2$mz_round == 331.039] <- "[M-2H]-2"
#HSO4
h3s2.df2$sugar[h3s2.df2$mz_round == 96.963 ] <- "HSO4"

h3s2.df2$label <- paste0(h3s2.df2$sugar, ": ", h3s2.df2$ion) %>% 
    sub(":\\s$", "", .)

sug <- c("HSO4", "hex-2-sulphate-1: [M-H]-", "hex-3-sulphate-2: [M-2H]-2",
         "hex-1-sulphate-1: [M-H]-", "hex-1-sulphate-2: [M-2H]-2", 
         "hex-1-sulphate-3: [M-3H]-3",
         "hex-4-sulphate-3: [M-3H]-3", "hex-7-sulphate-6: [M-6H]-6")
h3s2.df2$class <- ""
h3s2.df2$class[h3s2.df2$label == "hex-3-sulphate-2: [M-2H]-2"] <- "precursor"
h3s2.df2$class[h3s2.df2$label %in% c("HSO4", "hex-2-sulphate-1: [M-H]-",
                                     "hex-1-sulphate-1: [M-H]-", 
                                     "hex-1-sulphate-2: [M-2H]-2")] <- "expected"
h3s2.df2$class[h3s2.df2$class == "" &
                 h3s2.df2$label %in% sug] <- "not expected"
h3s2.df2$class <- factor(h3s2.df2$class, levels = c("precursor",
                                                    "expected", "not expected"))


png(filename = "GH99_plots_202202/NCE_vs_intensity_disulphated-mannotriose_v1.png",
    height = 6, width = 8, units = "in", res = 300)
ggplot(data = h3s2.df2 %>% 
         dplyr::filter(label %in% sug)) +
  geom_line(aes(x = nce, y = intensity, colour = class,
                group = label), lwd = 0.5) +
  geom_point(aes(x = nce, y = intensity, colour = class),
             shape = 21, fill = "white", size = 2) +
  facet_wrap(~label, scales = "free_y") +
  scale_colour_manual(values = pretty3, name = "Class") +
  theme_classic() +
  theme(text = element_text(family = "Avenir"),
        legend.position = "top")
dev.off()

#13: trisulphated mannotetraose------
h4s3.df2 <- msms2.df %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(precursorIon == "trisulphated mannotetraose: [M-3H]-3" &
                  !is.na(mz)) %>% 
  mutate(mz_round =  round(mz, 3))

#precursor
h4s3.df2$sugar[h4s3.df2$mz_round == 301.029 |
                   h4s3.df2$mz_round == 301.031 |
                   h4s3.df2$mz_round == 301.055 |
                   h4s3.df2$mz_round == 301.026] <- "hex-4-sulphate-3"
h4s3.df2$ion[h4s3.df2$mz_round == 301.029 |
                 h4s3.df2$mz_round == 301.031 |
                 h4s3.df2$mz_round == 301.055 |
                 h4s3.df2$mz_round == 301.026] <- "[M-3H]-3"
#HSO4
h4s3.df2$sugar[h4s3.df2$mz_round == 96.954 ] <- "HSO4"

h4s3.df2$label <- paste0(h4s3.df2$sugar, ": ", h4s3.df2$ion) %>% 
    sub(":\\s$", "", .)

#others
h4s3.df2$label[h4s3.df2$mz_round == 331.037] <- "hex-3-sulphate-2: [M-2H]-2"
h4s3.df2$label[h4s3.df2$mz_round == 247.002 |
                 h4s3.df2$mz_round == 247.004] <- "hex-3-sulphate-3: [M-3H]-3"


sug <- c("hex-4-sulphate-3: [M-3H]-3", "HSO4", "hex-3-sulphate-2: [M-2H]-2",
         "hex-1-sulphate-2: [M-2H]-2", "hex-3-sulphate-3: [M-3H]-3",
         "hex-5-sulphate-4: [M-4H]-4", "hex-2-sulphate-1: [M-H]-",
         "hex-1-sulphate-1: [M-H]-", 
         "hex-1-sulphate-3: [M-3H]-3",
         "hex-4-sulphate-3: [M-3H]-3", "hex-7-sulphate-6: [M-6H]-6",
         "hex-7-sulphate-4: [M-4H]-4", "hex-6-sulphate-5: [M-5H]-5")

h4s3.df2$class <- ""
h4s3.df2$class[h4s3.df2$label == "hex-4-sulphate-3: [M-3H]-3"] <- "precursor"
h4s3.df2$class[h4s3.df2$label %in% c("HSO4", "hex-2-sulphate-1: [M-H]-",
                                     "hex-1-sulphate-1: [M-H]-", 
                                     "hex-1-sulphate-2: [M-2H]-2",
                                     "hex-1-sulphate-3: [M-3H]-3",
                                     "hex-3-sulphate-2: [M-2H]-2",
                                     "hex-3-sulphate-3: [M-3H]-3")] <- "expected"
h4s3.df2$class[h4s3.df2$class == "" &
                 h4s3.df2$label %in% sug] <- "not expected"
h4s3.df2$class <- factor(h4s3.df2$class, levels = c("precursor",
                                                    "expected", "not expected"))


png(filename = "GH99_plots_202202/NCE_vs_intensity_trisulphated-mannotetraose_v1.png",
    height = 6, width = 10, units = "in", res = 300)
ggplot(data = h4s3.df2 %>% 
         dplyr::filter(label %in% sug)) +
  geom_line(aes(x = nce, y = intensity, colour = class,
                group = label), lwd = 0.5) +
  geom_point(aes(x = nce, y = intensity, colour = class),
             shape = 21, fill = "white", size = 2) +
  facet_wrap(~label, scales = "free_y") +
  scale_colour_manual(values = pretty3, name = "Class") +
  theme_classic() +
  theme(text = element_text(family = "Avenir"),
        legend.position = "top")
dev.off()

