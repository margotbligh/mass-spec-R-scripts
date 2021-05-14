setwd("~/Google_Drive/MPI_Masters/MSc_thesis/Lab_things/Experiments/1_standards/1C_GCC-SPE-yield")
load("analysis/RData/RData_final.RData")

#1: Install packages --------------------------------------------------------

library(data.table)
library(tidyverse)
library(ggplot2)
library(lsmeans)

library(wesanderson)
library(hablar)
library(patchwork)
library(RColorBrewer)
library(gplots)
library(grid)
library(gridExtra)
library(scales)
library(cowplot)
library(gtable)
library(ggpubr)
library(ggpmisc)
library(grid)
library(car)
library(compute.es)
library(effects)
library(multcomp)
library(pastecs)
library(WRSS)

library(rstatix)

#2: Import data --------------------------------------------------------

mmass_peaks.filepath <- dir(path = "./mmass-peak_lists", 
                            all.files = FALSE, 
                            full.names = TRUE)

#create sample list
samples.df <- data.frame(sampleID = basename(mmass_peaks.filepath) %>% 
                             sub(".txt", "", .),
                         concentration = basename(mmass_peaks.filepath) %>% 
                             sub("-[[:upper:]]+\\d*.txt", "", .) %>% 
                             sub("-", ".", .) %>%  
                             as.numeric(),
                         fraction = basename(mmass_peaks.filepath) %>% 
                             sub(".txt", "", .) %>% 
                             gsub("\\d+-", "", .))

#import data
files_list <- lapply(mmass_peaks.filepath, fread)

#create vector of sample IDs
sampleID.vec <- basename(mmass_peaks.filepath) %>% 
    sub(".txt", "", .)

#add column to each dataframe in list with sample ID    
files_list <- mapply(`[<-`, 
                     files_list, 
                     'sampleID', 
                     value = sampleID.vec, 
                     SIMPLIFY = FALSE)

#change column names in each dataframe
colnames <- c("mz", "intensity", "sampleID")
files_list <- lapply(files_list, setNames, colnames)

#3: Filter data for known ions-----------------
#read in table of masses
masses <- fread("GCC-SPE_yield_masses.csv")
masses$`m/z` <- as.numeric(masses$`m/z`)

#remove rows where m/z is 0
masses <- masses[-which(masses$`m/z`==0), ]

#combine all samples into one dataframe
master.df <- bind_rows(files_list)

#initialise empty columns
master.df$sugar <- NA
master.df$ion <- NA
master.df$theoretical_mz <- NA
master.df$mz_difference <- NA
master.df$type<- NA

#match known m/z values in table based on defined mass tolerance
error = 0.0025

for (i in 1:nrow(masses)){
    upper = masses$`m/z`[i] + error
    lower = masses$`m/z`[i] - error
    sugar = masses$ID[i]
    ion = masses$Ion[i]
    type = masses$Type[i]
    index <- which(master.df$mz >= lower & master.df$mz <= upper )
    master.df$sugar[index] <- sugar
    master.df$ion[index] <- ion
    master.df$type[index] <- type
    master.df$theoretical_mz[index] <- masses$`m/z`[i]
    master.df$mz_difference[index] <- abs(masses$`m/z`[i] - master.df$mz[index]) 
}


#build dataframe for ions
#cannot just filter, as need missing values to be NA (rather than just missing)
#when plotting
sugarID.vec.all <- masses$ID
ion.vec.all <- masses$Ion
type.vec.all <- masses$Type
theoretical_mz.vec.all <- masses$`m/z`
conc.vec.all <- samples.df$concentration
fraction.vec.all <- samples.df$fraction

ions <- data.frame(sampleID = rep(sampleID.vec, 
                                  each = length(ion.vec.all)),
                   fraction = rep(fraction.vec.all,
                                  each = length(ion.vec.all)),
                   concentration = rep(conc.vec.all,
                                       each = length(ion.vec.all)),
                   sugar = rep(sugarID.vec.all,
                               length(sampleID.vec)),
                   ion = rep(ion.vec.all,
                             length(sampleID.vec)),
                   theoretical_mz = rep(theoretical_mz.vec.all,
                                        length(sampleID.vec)),
                   type = rep(type.vec.all,
                              length(sampleID.vec)))

#use data.table merge function to fill in values
setDT(ions); setDT(master.df)
identified_ions <- merge(ions, 
                         master.df,
                         by = c("sampleID",
                                "sugar",
                                "ion",
                                "theoretical_mz",
                                "type"),
                         all.x = TRUE)

#check class of variables
#lapply(identified_ions,class)

#see how many signals are duplicated for each sample
#should have no duplicates
duplicates <- identified_ions[duplicated(identified_ions, 
                                         by = c("sugar", 
                                                "ion", 
                                                "fraction", 
                                                "concentration")),]

if (nrow(duplicates) == 0) {
    print("no duplicates") ; rm(duplicates)
} else if (nrow(duplicates) != 0) {
        print("duplicates found!")
    }

#account for dilution factors
x <- c(1, 2.5, 5, 10, 25)
identified_ions$intensity[identified_ions$concentration %in% x] <- identified_ions$intensity[identified_ions$concentration %in% x] * 0.1
identified_ions$intensity[identified_ions$concentration == 2500] <- identified_ions$intensity[identified_ions$concentration == 2500] * 10

#format ion names
identified_ions$ion_fmt <- identified_ions$ion %>%
    sub("M", "[M", .) %>% 
    sub("-H", "-H]-", .) %>% 
    sub("-2H", "-2H]2-", .) %>% 
    sub("-3H", "-3H]3-", .) %>% 
    sub("-4H", "-4H]4-", .) %>% 
    sub("-5H", "-5H]5-", .) %>% 
    sub("\\+Cl", "+Cl]-", .) %>%
    sub("\\+CH00", "+CHOO]-", .)

#format sugar names
identified_ions$sugar_fmt <- identified_ions$sugar %>% 
    sub("BM3","b-mannotriose", .) %>% 
    sub("C6","cellohexaose", .) %>% 
    sub("G1", "glucose", .) %>% 
    sub("L2", "laminaribiose", .) %>% 
    sub("L4", "laminaritetraose", .) %>% 
    sub("KI", "k-i-k-carrageenan DP", .) %>% 
    sub("K", "k-carrageenan DP", .) %>% 
    sub("MA1", "mannuronic acid", .) %>% 
    sub("MA3", "tri-mannuronic acid", .) %>% 
    sub("MA5", "penta-mannuronic acid", .)

masses$sugar_fmt <- masses$ID %>% 
    sub("BM3","b-mannotriose", .) %>% 
    sub("C6","cellohexaose", .) %>% 
    sub("G1", "glucose", .) %>% 
    sub("L2", "laminaribiose", .) %>% 
    sub("L4", "laminaritetraose", .) %>% 
    sub("KI", "k-i-k-carrageenan DP", .) %>% 
    sub("K", "k-carrageenan DP", .) %>% 
    sub("MA1", "mannuronic acid", .) %>% 
    sub("MA3", "tri-mannuronic acid", .) %>% 
    sub("MA5", "penta-mannuronic acid", .)

#sum all ions for a each sugar in each sample
#need to first change NA to 0, then set back to NA
identified_ions[is.na(identified_ions)] <- 0

identified_ions.sums <- identified_ions %>% 
    group_by(sampleID, sugar_fmt, concentration, fraction, type) %>% 
    summarise(intensity = sum(intensity)) %>% 
    ungroup()

identified_ions.sums$intensity[identified_ions.sums$intensity == 0] <- NA

#get concatenated list of which ions were summed for each sugar within each sample
ions_temp <- identified_ions[identified_ions$intensity!= 0,]
ions_temp$ion_mz <- paste0(ions_temp$ion_fmt,
                           " mz=",
                           sprintf("%.3f",ions_temp$mz)) 

ions_temp_2 <- ions_temp[, 
                         list(ion = paste(ion_mz, 
                                          collapse = ' + ')), 
                         by = c("sampleID", 
                                "sugar_fmt", 
                                "concentration", 
                                "fraction",
                                "type")]
identified_ions.sums_noIons <- identified_ions.sums #keep to be safe

identified_ions.sums <- merge(identified_ions.sums,
                              ions_temp_2,
                              by = c("sampleID", 
                                     "sugar_fmt", 
                                     "concentration", 
                                     "fraction",
                                     "type"),
                              all.x = TRUE)

rm(ions_temp,
   ions_temp_2)

#format fraction names
identified_ions.sums$fraction_fmt <- identified_ions.sums$fraction %>% 
    sub("I", "initial", .) %>% 
    sub("SF", "sample flow-through", .) %>% 
    sub("W", "wash", .) %>% 
    sub("E", "eluate ", .)

#split into blanks and without blanks
noblanks <- identified_ions.sums %>% 
    filter(!grepl("B", fraction))

blanks <- identified_ions.sums %>% 
    filter(grepl("B", fraction))

#4: Regression ----
#a: Check for linearity between log2(concentration) and log2(intensity) ----
sugar.vec <- masses$sugar_fmt %>% 
    unique()

fraction_levels <- c("initial", 
                     "sample flow-through", 
                     "wash", 
                     "eluate 1", 
                     "eluate 2", 
                     "eluate 3")

lm.df <- data.frame(sugar = rep(sugar.vec, 
                                each = length(fraction_levels)),
                    fraction = rep(fraction_levels, 
                                   length(sugar.vec)),
                    c = rep("", 
                            length(sugar.vec)*length(fraction_levels)),
                    m = rep("", 
                            length(sugar.vec)*length(fraction_levels)),
                    r2 = rep("", 
                             length(sugar.vec)*length(fraction_levels)),
                    pval = rep("",
                               length(sugar.vec)*length(fraction_levels)))
#linear regression
#only do if >= three data points
for (i in 1:length(sugar.vec)){
    sug = sugar.vec[i]
    for (j in 1:length(fraction_levels)){
        frac = fraction_levels[j]
        x <- noblanks %>% 
            filter(sugar_fmt == !!sug,
                   fraction_fmt == !!frac)
        if(10 - length(which(is.na(x$intensity))) < 3) {next}
        else {
            fit <- lm(log2(intensity) ~ log2(concentration), 
                      data = x)
            c = summary(fit)$coefficients[1,1]
            m = summary(fit)$coefficients[2,1]
            pval = summary(fit)$coefficients[2,4]
            r2 = summary(fit)$adj.r.squared
            lm.df$c[lm.df$sugar == sug & 
                        lm.df$fraction == frac] <- c
            lm.df$m[lm.df$sugar == sug & 
                        lm.df$fraction == frac] <- m
            lm.df$r2[lm.df$sugar == sug & 
                         lm.df$fraction == frac] <- r2
            lm.df$pval[lm.df$sugar == sug & 
                           lm.df$fraction == frac] <- pval
        }
    }
}

lm.df$c <- as.numeric(lm.df$c)
lm.df$m <- as.numeric(lm.df$m)
lm.df$r2 <- as.numeric(lm.df$r2)
lm.df$pval <- as.numeric(lm.df$pval)

#see which sugars have significant relationships in initial
signif_sug <- lm.df %>% 
    filter(fraction == "initial") %>% 
    filter(pval <= 0.05) %>% 
    filter(r2 >= 0.65 ) 
signif_sug.vec <- signif_sug$sugar

#keep only these sugars, and only fractions with significant regressions
lm.df.sigSug <- lm.df %>% 
    filter(sugar %in% signif_sug.vec) %>% 
    filter(pval <= 0.05) %>% 
    filter(r2 >= 0.65 ) 
    

#b: Check for homogeneity of slopes, within each sugar ----
#create output directories
dir.create(file.path("analysis"), showWarnings = TRUE)
dir.create(file.path("analysis/linear_models"), showWarnings = TRUE)
dir.create(file.path("analysis/linear_models/assumption_plots"), showWarnings = TRUE)


lm.df.sigSug.aov <- data.frame(sugar = signif_sug.vec,
                               interaction = rep(NA, 
                                                 length(signif_sug.vec)),
                               interaction.pval = rep(NA, 
                                                      length(signif_sug.vec)),
                               differences = rep(NA, 
                                                 length(signif_sug.vec)))

for (i in 1:length(signif_sug.vec)){
    sug = signif_sug.vec[i]
    
    #check that there is more than one fraction
    d <- noblanks %>% 
        filter(sugar_fmt == !!sug) %>%
        filter(fraction_fmt %in% 
                   lm.df.sigSug$fraction[lm.df.sigSug$sugar == !!sug])
    
    if (length(unique(d$fraction_fmt)) <= 1 ) {
        next
    } else { 
        #get linear model
        x <- d %>% 
            lm(log2(intensity)~log2(concentration)*fraction_fmt, .)
        
        #check assumptions
        tiff(paste0("./analysis/linear_models/assumption_plots/",
                    sug,
                    "_density-residuals.tiff"),
             res = 300,
             units  = "in",
             width = 3,
             height = 3)
        plot(density(x$residuals))
        dev.off()
        tiff(paste0("./analysis/linear_models/assumption_plots/",
                    sug,
                    "_qqplot.tiff"),
             res = 300,
             units  = "in",
             width = 3,
             height = 3)
        qqnorm(x$residuals)
        qqline(x$residuals, datax = FALSE, 
               distribution = qnorm, 
               probs = c(0.25, 0.75))
        dev.off()
        
        #check for interaction
        aov <- anova(x)
        rownames(aov) <- c("log2(concentration)", 
                           "fraction_fmt",
                           "interaction",
                           "residuals")
        pval.int <- aov["interaction", "Pr(>F)"]
        if (pval.int > 0.05) {
            lm.df.sigSug.aov$interaction[
                lm.df.sigSug.aov$sugar == sug] <- "no"
            lm.df.sigSug.aov$interaction.pval[
                lm.df.sigSug.aov$sugar == sug] <- pval.int
            next
        } else { 
            lm.df.sigSug.aov$interaction[
                lm.df.sigSug.aov$sugar == sug] <- "yes"
            lm.df.sigSug.aov$interaction.pval[
                lm.df.sigSug.aov$sugar == sug] <- pval.int
        
            #obtain slopes
            m.lst <- lstrends(x, "fraction_fmt", var="log2(concentration)")
            #compare slopes
            m.comp <- pairs(m.lst)
            m.sum <- summary(m.comp)
            pairs.diff <- paste(m.sum$contrast[m.sum$p.value<0.05], 
                                collapse = ', ')
            lm.df.sigSug.aov$differences[
                lm.df.sigSug.aov$sugar == sug] <- pairs.diff
        }
    }
}


#5:Save RData ----
#create directory
dir.create(file.path("analysis/RData"), showWarnings = TRUE)

#save RData
save.image("analysis/RData/RData_final.RData")







