#THIS IS A "CLEAN" WORKSPACE AND SCRIPT FOR MARGOT TO WORK ON GC-MS DATA
#objects are imported from workspace generated from previous script

#starting point:
#pseudospectra from xsAnnotate object with >= 3 associated features (n = 389)
#each pseudospectra = one pcgroup from CAMERA grouping of features
#I already know that ribitol is number 23 in the list of psuedospectra, and
#the ion with max intensity is 

setwd("/Users/margotbligh/Google_Drive/MPI_Masters/Assemble")

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
library(MSnbase)

#2: Load data objects----
##peaklist
load("./RData/peaklist.RData")
load("./RData/pksInt.RData")
##pd
load("./RData/pd.RData")
##mona spec
load("./RData/mona.spec.RData")
load("./RData/mona.names.RData")
##an
load("./RData/an.RData")
##pspectra
load("./RData/pspectra.RData")
load("./RData/spclist.by.sample.RData")
load("./RData/spclist.by.sample.spec.norm.RData")
load("./RData/rt.spclist.RData")

#3: Normalise all ions by ribitol intensity (per sample)----
for(i in 1:length(spclist.by.sample)){
    x <- spclist.by.sample[[i]] #subset
    x <- sapply(x, function(x) as.data.frame(x)) #make df
    #get those with only 1 row (weirdly form 3x1 df)
    n <- which(sapply(spclist.by.sample[[i]], function(y) is.null(nrow(y))))
    for(j in 1:length(n)){
        index = n[j]
        x[[index]] <- t(x[[index]]) #transpose
    }
    x <- lapply(x, transform, norm_intensity = intensity/r) #normalise
    spclist.by.sample[[i]] <- x #add back in
}

#4: Get quantification ions----
#set up df
pspectra.q.ions <- data.frame(ps_nr = seq(1, length(pspectra),1),
                              pcgroup = names(spclist.by.sample[[1]]),
                              quant_mz = NA, max_sample = NA)
pspectra.q <- data.frame(sample_name = rep(names(pksInt), length(pspectra)),
                         sample_nr = rep(1:nrow(pd), length(pspectra)),
                         ps_nr = rep(1:length(pspectra), each = nrow(pd)),
                         quant_mz = NA, intensity = NA, norm_intensity = NA)
pname <- names(pspectra.q)

for(i in 1:length(pspectra)){
    #subset list and bind rows
    x <- bind_rows(lapply(spclist.by.sample, "[[", i), .id = "sample_nr")
    #get max ion
    q <- x$mz[which(x$norm_intensity == max(x$norm_intensity))]
    pspectra.q.ions$quant_mz[i] <- q
    pspectra.q.ions$max_sample[i] <- x$sample_nr[
        which(x$norm_intensity == max(x$norm_intensity))]
    pspectra.q$quant_mz[pspectra.q$ps_nr == i] <- q
    
    #add into table
    x <- x[x$mz == q,]
    x$quant_mz <- q
    x$ps_nr <- i
    x <- x[,c("sample_nr","ps_nr","quant_mz", "intensity", "norm_intensity")]
    x$sample_nr <- as.numeric(x$sample_nr)
    pspectra.q <- left_join(pspectra.q, x, by = c("sample_nr", "quant_mz", "ps_nr")) %>% 
        mutate(intensity = coalesce(intensity.x, intensity.y),
               norm_intensity = coalesce(norm_intensity.x, norm_intensity.y)) %>% 
        dplyr::select(all_of(pname))
}

#format name
pd$sample_name_original <- pd$sample_name
pd$sample_name <- pd$sample_name_original %>% sub("\\s-\\s", "...", .)
pspectra.q.met <- left_join(pspectra.q, pd, by = "sample_name")
#remove MQ (normalisation gives crazy values because no ribitol added)
pspectra.q.met.wMQ <- pspectra.q.met
pspectra.q.met <- pspectra.q.met %>% filter(sample_group != "MQ_NA_NA")


#5: Significance testing----
#T-tests withtin groups----
#set up grouping variables
pspectra.q.met$sample_group2 <- paste0(pspectra.q.met$sample_group, "_",
                                       pspectra.q.met$sampling_time) %>% 
    gsub("_dark|_light|_NA", "", .)
pspectra.q.met$sample_group <- pspectra.q.met$sample_group2 %>% 
    sub("_I|_E|_blank|_metab.*", "", .) 
pspectra.q.met$sampling_time[pspectra.q.met$sample_group2=="ASW_blank"] <- "I"
pspectra.q.met$sampling_time[pspectra.q.met$sample_group2=="ASW_metabolitestd"] <- "E"

#get column names for table
i = 1
ps_nr = i
df1 <- pspectra.q.met %>% filter(ps_nr == !!ps_nr)
cols <- df1 %>%
    rstatix::group_by(sample_group) %>%
    rstatix::t_test(norm_intensity ~ sampling_time) %>%
    rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
    rstatix::add_significance(p.col = "p.adj") %>% 
    rstatix::add_xy_position(x = "sampling_time", dodge = 0.8)
cols <- cbind(ps_nr = NA, cols)
ttests.df <- cols[!any(is.na(cols)),]
rm(cols)

#do t-tests and p-value adjustment
for(i in 1:length(pspectra)){
    ps_nr = i
    df1 <- pspectra.q.met %>% filter(ps_nr == !!ps_nr)
    df1$intensity[is.na(df1$intensity)] <- 0
    df1$norm_intensity[is.na(df1$norm_intensity)] <- 0
    df2 <- df1 %>%
        rstatix::group_by(sample_group) %>%
        rstatix::t_test(norm_intensity ~ sampling_time) %>%
        rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
        rstatix::add_significance(p.col = "p.adj") %>% 
        rstatix::add_xy_position(x = "sampling_time", dodge = 0.8)
    df2 <- cbind(ps_nr = i, df2)
    ttests.df <- rbind(ttests.df, df2)
}

ttests.df.sig <- ttests.df %>% filter(p.adj <= 0.05)
ttests.df.sig$ps_nr %>% unique() %>% length() #n = 101

#Wilcox tests between groups and ASW blank ----
#get column names for table
i = 1
ps_nr = i
df1 <- pspectra.q.met %>% filter(ps_nr == !!ps_nr)
df1 <- df1 %>% filter(sample_group2 != "ASW_metabolitestd")
df1$norm_intensity[is.na(df1$norm_intensity)] <- 0
df1 <- df1 %>% 
    rstatix::wilcox_test(norm_intensity ~ sample_group, ref.group = "ASW") %>% 
    rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
    rstatix::add_significance(p.col = "p.adj") %>% 
    rstatix::add_xy_position()
cols <- cbind(ps_nr = NA, df1)
wilcox.df <- cols[!any(is.na(cols)),]
rm(cols)

#do t-tests and p-value adjustment
for(i in 1:length(pspectra)){
    ps_nr = i
    df1 <- pspectra.q.met %>% filter(ps_nr == !!ps_nr)
    df1 <- df1 %>% filter(sample_group2 != "ASW_metabolitestd")
    df1$norm_intensity[is.na(df1$norm_intensity)] <- 0
    k <- df1 %>% group_by(sample_group) %>% summarise(sum = sum(norm_intensity))
    if (nrow(k[k$sum == 0,]) >= 2){next}
    df2 <- df1 %>% 
        rstatix::wilcox_test(norm_intensity ~ sample_group, ref.group = "ASW") %>% 
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

#6: Plot pseudospec with significant differences----
pspectra.q.met.sig <- pspectra.q.met %>% filter(ps_nr %in% ttests.df.sig$ps_nr |
                                                    ps_nr %in% wilcox.df.sig$ps_nr)
a <- pspectra.q.met.sig$ps_nr %>% unique() #n = 163

pspectra.q.met.sig$sampling_time <- factor(pspectra.q.met.sig$sampling_time,
                                           levels = c("I", "E"))

ttests.df.sig$xmin[ttests.df.sig$sample_group == "ASW"] <- 0.8
ttests.df.sig$xmax[ttests.df.sig$sample_group == "ASW"] <- 1.2
ttests.df.sig$xmin[ttests.df.sig$sample_group == "PW_SEAGRASS"] <- 1.8
ttests.df.sig$xmax[ttests.df.sig$sample_group == "PW_SEAGRASS"] <- 2.2
ttests.df.sig$xmin[ttests.df.sig$sample_group == "WC_ALGAE"] <- 2.8
ttests.df.sig$xmax[ttests.df.sig$sample_group == "WC_ALGAE"] <- 3.2
ttests.df.sig$xmin[ttests.df.sig$sample_group == "WC_SEAGRASS"] <- 3.8
ttests.df.sig$xmax[ttests.df.sig$sample_group == "WC_SEAGRASS"] <- 4.2

wilcox.df.sig$xmin <- 0.8


w <- wilcox.df.sig$ps_nr %>% unique()
t <- ttests.df.sig$ps_nr %>% unique()


for (i in 1:length(a)){
    ps_nr = a[i]
    df1 <- pspectra.q.met.sig %>% filter(ps_nr == !!ps_nr)
    df2 <- ttests.df.sig%>% filter(ps_nr == !!ps_nr)
    df3 <- wilcox.df.sig%>% filter(ps_nr == !!ps_nr)
    if (ps_nr %in% w & ps_nr %in% t){
        p <- ggplot(data = df1, aes(x=sample_group, y=sqrt(norm_intensity))) +
            geom_boxplot(aes(fill=sampling_time), position = position_dodge(width = 0.8)) + 
            geom_point(pch=21, size=2,mapping =aes(fill=sampling_time), 
                       position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("#FEC000","#2BB6AF"),
                              name = "Sampling time or sample type") +
            scale_y_continuous(name=expression(sqrt(Ribitol~normalised~intensity~(a.u.)))) +
            labs(x = "Sample group") +
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
                  panel.grid = element_blank(),
                  legend.position = "top",
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "#262626"),
                  panel.spacing = unit(0.1, "lines"))
        png(filename = paste0("test-plots-3/both/", ps_nr, ".png"),
            width = 18, height = 9, res = 300, units = "in")
        print(p)
        dev.off()
    }
    else if (ps_nr %in% w & !ps_nr %in% t){
        p <- ggplot(data = df1, aes(x=sample_group, y=sqrt(norm_intensity))) +
            geom_boxplot(aes(fill=sampling_time), position = position_dodge(width = 0.8)) + 
            geom_point(pch=21, size=2,mapping =aes(fill=sampling_time), 
                       position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("#FEC000","#2BB6AF"),
                              name = "Sampling time or sample type") +
            scale_y_continuous(name=expression(sqrt(Ribitol~normalised~intensity~(a.u.)))) +
            labs(x = "Sample group") +
            add_pvalue(df3, xmin = "xmin", xmax = "xmax", 
                       y.position = sqrt(df3$y.position),
                       label = "p.adj.signif",
                       label.size = 5, fontfamily = "Arial",
                       step.increase = 0.1) +
            theme_bw() + 
            theme(text=element_text(size=12, family = "Arial", colour = "#262626"),
                  panel.grid = element_blank(),
                  legend.position = "top",
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "#262626"),
                  panel.spacing = unit(0.1, "lines"))
        png(filename = paste0("test-plots-3/between-ASW-and-groups-only/", 
                              ps_nr, ".png"),
            width = 18, height = 9, res = 300, units = "in")
        print(p)
        dev.off()
    }
    else if (!ps_nr %in% w & ps_nr %in% t){
        p <- ggplot(data = df1, aes(x=sample_group, y=sqrt(norm_intensity))) +
            geom_boxplot(aes(fill=sampling_time), position = position_dodge(width = 0.8)) + 
            geom_point(pch=21, size=2,mapping =aes(fill=sampling_time), 
                       position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("#FEC000","#2BB6AF"),
                              name = "Sampling time or sample type") +
            scale_y_continuous(name=expression(sqrt(Ribitol~normalised~intensity~(a.u.)))) +
            labs(x = "Sample group") +
            add_pvalue(df2, xmin = "xmin", xmax = "xmax", 
                       y.position = sqrt(df2$y.position),
                       label = "p.adj.signif",
                       label.size = 5, fontfamily = "Arial") +
            theme_bw() + 
            theme(text=element_text(size=12, family = "Arial", colour = "#262626"),
                  panel.grid = element_blank(),
                  legend.position = "top",
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "#262626"),
                  panel.spacing = unit(0.1, "lines"))
        
        if (nrow(df2) == 1 & df2$sample_group == "ASW"){
            png(filename = paste0("test-plots-3/within-groups-only/std-only/", 
                                  ps_nr, ".png"),
                width = 18, height = 9, res = 300, units = "in")
            print(p)
            dev.off()
        }
        else{
            png(filename = paste0("test-plots-3/within-groups-only/samples/", 
                                  ps_nr, ".png"),
                width = 18, height = 9, res = 300, units = "in")
            print(p)
            dev.off()
        }
    }
}


#7: Get annotation for interesting pseudospectra----
#make spectra centroided (required to for plot comparison)
spclist.by.sample.spec.norm.uncentroided <- spclist.by.sample.spec.norm
for (i in 1:length(spclist.by.sample.spec.norm.uncentroided)){
    x <- spclist.by.sample.spec.norm.uncentroided[[i]]
    x <- lapply(x, function(x) new("Spectrum1", 
                                   mz = x@mz,
                                   intensity = x@intensity,
                                   centroided = TRUE))
    spclist.by.sample.spec.norm[[i]] <- x
}


#start with those that are both higher in samples than ASW, and have a within group significance
a <- dir(path = "test-plots-3/both", all.files = FALSE)
a <- a %>% sub(".png", "", .) %>% as.numeric()
pspectra.sigdiff1 <- data.frame(ps_nr = NA, annot1 = NA, dotprod1 = NA,
                               annot2 = NA,dotprod2 = NA,
                               annot3 = NA,dotprod3 = NA,
                               annot4 = NA, dotprod4 = NA,
                               annot5 = NA, dotprod5 = NA)

for (i in 1:length(a)){
    ps_nr = a[i]
    sample_nr = as.numeric(pspectra.q.ions$max_sample[
        pspectra.q.ions$ps_nr == ps_nr])
    l <- lapply(mona.spec,
                function(x)
                    compareSpectra(x,
                                   spclist.by.sample.spec.norm[[sample_nr]][[ps_nr]],
                                   fun = "dot"))
    l <- unlist(l)
    names(l) <- 1:length(l)
    l <- sort(l, decreasing = TRUE)
    tmp.df <- data.frame(ps_nr = ps_nr,
                         annot1 = mona.names[as.numeric(names(l)[1])], dotprod1 = l[1],
                         annot2 = mona.names[as.numeric(names(l)[2])], dotprod2 = l[2],
                         annot3 = mona.names[as.numeric(names(l)[3])], dotprod3 = l[3],
                         annot4 = mona.names[as.numeric(names(l)[4])], dotprod4 = l[4],
                         annot5 = mona.names[as.numeric(names(l)[5])], dotprod5 = l[5])
    pspectra.sigdiff1 <- rbind(pspectra.sigdiff1, tmp.df)
    
}

x <- pspectra.sigdiff1
pspectra.sigdiff1 <- pspectra.sigdiff1[complete.cases(pspectra.sigdiff1),]
pspectra.sigdiff1 <- pspectra.sigdiff1[order(pspectra.sigdiff1$ps_nr),]
fwrite(pspectra.sigdiff1, file = "results/pspectra.sigdiff1.txt", sep = " ")


#those that have a within group significance in samples
a <- dir(path = "test-plots-3/within-groups-only/samples", all.files = FALSE)
a <- a %>% sub(".png", "", .) %>% as.numeric()
pspectra.sigdiff2 <- data.frame(ps_nr = NA, annot1 = NA, dotprod1 = NA,
                                annot2 = NA,dotprod2 = NA,
                                annot3 = NA,dotprod3 = NA,
                                annot4 = NA, dotprod4 = NA,
                                annot5 = NA, dotprod5 = NA)

for (i in 1:length(a)){
    ps_nr = a[i]
    sample_nr = as.numeric(pspectra.q.ions$max_sample[
        pspectra.q.ions$ps_nr == ps_nr])
    l <- lapply(mona.spec,
                function(x)
                    compareSpectra(x,
                                   spclist.by.sample.spec.norm[[sample_nr]][[ps_nr]],
                                   fun = "dot"))
    l <- unlist(l)
    names(l) <- 1:length(l)
    l <- sort(l, decreasing = TRUE)
    tmp.df <- data.frame(ps_nr = ps_nr,
                         annot1 = mona.names[as.numeric(names(l)[1])], dotprod1 = l[1],
                         annot2 = mona.names[as.numeric(names(l)[2])], dotprod2 = l[2],
                         annot3 = mona.names[as.numeric(names(l)[3])], dotprod3 = l[3],
                         annot4 = mona.names[as.numeric(names(l)[4])], dotprod4 = l[4],
                         annot5 = mona.names[as.numeric(names(l)[5])], dotprod5 = l[5])
    pspectra.sigdiff2 <- rbind(pspectra.sigdiff2, tmp.df)
    
}

x <- pspectra.sigdiff2
pspectra.sigdiff2 <- pspectra.sigdiff2[complete.cases(pspectra.sigdiff2),]
pspectra.sigdiff2 <- pspectra.sigdiff2[order(pspectra.sigdiff2$ps_nr),]
fwrite(pspectra.sigdiff2, file = "results/pspectra.sigdiff2.txt", sep = " ")


#those that have a between group significance 
a <- dir(path = "test-plots-3/between-ASW-and-groups-only", all.files = FALSE)
a <- a %>% sub(".png", "", .) %>% as.numeric()
pspectra.sigdiff3 <- data.frame(ps_nr = NA, annot1 = NA, dotprod1 = NA,
                                annot2 = NA,dotprod2 = NA,
                                annot3 = NA,dotprod3 = NA,
                                annot4 = NA, dotprod4 = NA,
                                annot5 = NA, dotprod5 = NA)

for (i in 1:length(a)){
    ps_nr = a[i]
    sample_nr = as.numeric(pspectra.q.ions$max_sample[
        pspectra.q.ions$ps_nr == ps_nr])
    l <- lapply(mona.spec,
                function(x)
                    compareSpectra(x,
                                   spclist.by.sample.spec.norm[[sample_nr]][[ps_nr]],
                                   fun = "dot"))
    l <- unlist(l)
    names(l) <- 1:length(l)
    l <- sort(l, decreasing = TRUE)
    tmp.df <- data.frame(ps_nr = ps_nr,
                         annot1 = mona.names[as.numeric(names(l)[1])], dotprod1 = l[1],
                         annot2 = mona.names[as.numeric(names(l)[2])], dotprod2 = l[2],
                         annot3 = mona.names[as.numeric(names(l)[3])], dotprod3 = l[3],
                         annot4 = mona.names[as.numeric(names(l)[4])], dotprod4 = l[4],
                         annot5 = mona.names[as.numeric(names(l)[5])], dotprod5 = l[5])
    pspectra.sigdiff3 <- rbind(pspectra.sigdiff3, tmp.df)
    
}

x <- pspectra.sigdiff3
pspectra.sigdiff3 <- pspectra.sigdiff3[complete.cases(pspectra.sigdiff3),]
pspectra.sigdiff3 <- pspectra.sigdiff3[order(pspectra.sigdiff3$ps_nr),]
fwrite(pspectra.sigdiff3, file = "results/pspectra.sigdiff3.txt", sep = " ")

