
# Deconvolution --------------------------------------------------------
#start with xset from GCMS-data-analysis-script_ASSEMBLE.r

# checking with Margot 14. May 2021

setwd("C:/Users/admin/ownCloud/Assemble/Data/SeaMet/GCMS data analysis/results/")
load(file = "RData_metaMStesting_20210511.RData")


plot(cholestane.sp)
plot(ribitol.sp)




# create CAMERA object
an <- xsAnnotate(xset)
#peak grouping by rt
an <- groupFWHM(an)
# find isotopic peaks
an <- findIsotopes(an, mzabs=0.01)
# verify peak groups
an <- groupCorr(an, 
                cor_eic_th=0.75)
# Convert pseudospectra to msp format-----
#metaMS has function "to.msp" but gives error:
#Error in allpks[, intensity] : subscript out of bounds
#"to.msp" written for list of xset objects with each object for one sample
#I tried to adapt function code as best I could to generate same result

#define settings
object <- an
ndigit = 0
minfeat = 5

class(an) == "xsAnnotate" #should be true

#get intensities of samples as dataframe
allpks <- object@groupInfo
allpks <- as.data.frame(allpks)
allpksInt <- select(allpks, contains("Copy"))

#filter pseudospectra to have minimum 5 associated features
pspectra <- object@pspectra
npeaks <- sapply(pspectra, length)
pspectra <- pspectra[npeaks >= minfeat]


#this makes a nested list as is the output from the metaMS wrappers
#each sample is a single element - within that each pseudospetra is an element
spclist.by.sample <- vector("list", length(pd$sample_name))
for (i in 1:length(pd$sample_name)) {
    spclist.by.sample[[i]] <- lapply(pspectra, 
                                   function(x)
                                       cbind(mz = round(allpks[x,"mz"], 
                                                        digits = ndigit),
                                             intensity = allpksInt[x, i],
                                             rt = allpks$rt[x]))
}

#construct msp objects (one per sample)
list1 <- vector("list", 264)
extra.info.list <- list()
extra.info.list[[1]] <- list1 #empty list, strcucture needed for construct.msp()

msp.list <- vector("list", 75)
for (i in 1:length(pd$sample_name)){
    msp.list[[i]] <- metaMS::construct.msp(spclist.by.sample[[i]],
                                           extra.info = extra.info.list)
}

#write msp objects to file
dir.create("./results/msp-by-sample",
           showWarnings = FALSE)

for (i in 1:length(pd$sample_name)){
    metaMS::write.msp(msp.list[[i]], 
                      file = paste0("./results/msp-by-sample/",
                                    pd$sample_name[i],
                                    ".msp"),
                      newFile = TRUE)
}

# Checking standards -----

rt.spclist <- sapply(spclist, "[", 1, 77) # subset from ps list, row 1, column 77 (contains rt of ps)

#ribitol: rt 950-1020, mz_ribitol 217-217.4
ribitol.msp <- read.msp("./msp-mona-singlecompounds/ribitol_MoNA002421.msp")

ribitol.sp <- mona.spec[which(mona.names=="Ribitol")][[1]]

ribitol.sp <- new("Spectrum1", 
                  mz = ribitol.msp[[1]]$pspectrum[,1], 
                  intensity = ribitol.msp[[1]]$pspectrum[,2])

index <- which(rt.spclist > 950 & rt.spclist < 1020)
for (i in 1:length(index)){
    n = index[i]
    if(217 %in% spclist.by.sample[[17]][[n]][,1]){ # in spclist.by.sample[[sample]][[ps]]
        print(n)
    }
}

#alt 1 

dot_ribitol = compareSpectra(ribitol.sp,
                             spclist.by.sample.spec.norm[[20]][[22]],
                             fun = "dot")
plot(ribitol.sp,
     spclist.by.sample.spec.norm[[20]][[22]],     
     main = paste0("Ribitol: dot product = ", dot_ribitol))

#alt 2


samp13_ribitol <- as.data.frame(spclist.by.sample[[13]][[22]])
samp13_ribitol$intensity <- samp13_ribitol$intensity / 
    max(samp13_ribitol$intensity)*100

samp13_ribitol.sp <- new("Spectrum1", 
                  mz = samp13_ribitol$mz, 
                  intensity = samp13_ribitol$intensity)

centroided(ribitol.sp) <- TRUE
centroided(samp13_ribitol.sp) <- TRUE

dot_ribitol = compareSpectra(ribitol.sp,
                             samp13_ribitol.sp,
                             fun = "dot")

plot(ribitol.sp,
     samp13_ribitol.sp,     
     main = paste0("Ribitol: dot product = ", dot_ribitol))

#cholestane: rt 1610-1630, mz_ribitol 373
cholestane.msp <- read.msp("./msp-mona-singlecompounds/cholestane_JP004122.msp")
cholestane.sp <- new("Spectrum1", 
                  mz = cholestane.msp[[1]]$pspectrum[,1], 
                  intensity = cholestane.msp[[1]]$pspectrum[,2])


index <- which(rt.spclist > 1610 & rt.spclist < 1630)
for (i in 1:length(index)){
    n = index[i]
    if(373 %in% spclist.by.sample[[13]][[n]][,1]){
        print(n)
    }
}

samp13_cholestane <- as.data.frame(spclist.by.sample[[13]][[17]])
samp13_cholestane <- na.omit(samp13_cholestane)
samp13_cholestane$intensity <- samp13_cholestane$intensity / 
    max(samp13_cholestane$intensity)*100

samp13_cholestane.sp <- new("Spectrum1", 
                         mz = samp13_cholestane$mz, 
                         intensity = samp13_cholestane$intensity)

centroided(cholestane.sp) <- TRUE
centroided(samp13_cholestane.sp) <- TRUE

dot_cholestane = compareSpectra(cholestane.sp,
                                samp13_cholestane.sp,
                                fun = "dot")
plot(samp13_cholestane.sp,
     cholestane.sp,     
     main = paste0("Cholestane: dot product = ", dot_cholestane))


# Compare pseudospectra to MoNA database ---------------
#read in external database 
mona.msp <- metaMS::read.msp(file = "MoNA-export-GC-MS_Spectra.msp")
mona.names <- unlist(sapply(mona.msp, "[", 1))

#example how to filter database

finding_cholestane<-mona.msp[which(mona.names == "Ribitol")]

which(mona.names == "Ribitol")


#create spectrum objects for each db pseudospectrum
mona.spec <- lapply(mona.msp,
                    function(x) new("Spectrum1", 
                                    mz = x$pspectrum[,1], 
                                    intensity = x$pspectrum[,2],
                                    centroided = TRUE))

#compare spectra within database

compareSpectra(mona.spec[[12125]], mona.spec[[15061]], fun="dot")



#change NA to 0 in sample spectra
spclist.by.sample.noNA <- vector("list", length(pd$sample_name))
for (i in 1:length(pd$sample_name)){
    spclist.by.sample.noNA[[i]] <- lapply(spclist.by.sample[[i]],
                                          function(x) replace_na(x,0))
}
#create spectrum objects for each sample pseudospectrum
spclist.by.sample.spec <-vector("list", length(pd$sample_name))
for (i in 1:length(pd$sample_name)) {
        spclist.by.sample.spec[[i]] <- lapply(spclist.by.sample.noNA[[i]],
                                              function(x) new("Spectrum1", 
                                                              mz = x[,1], 
                                                              intensity = x[,2],
                                                              centroided = TRUE))
}
#normalise
spclist.by.sample.spec.norm <-vector("list", length(pd$sample_name))
for (i in 1:length(pd$sample_name)) {
    for (j in 1:length(pspectra)){
        y <- spclist.by.sample.spec[[i]][[j]]
        if(sum(y@intensity) == 0) {
            next
        } 
        else if (sum(y@intensity) > 0){
            spclist.by.sample.spec.norm[[i]][[j]] <- 
                normalise(spclist.by.sample.spec[[i]][[j]],"max")
            spclist.by.sample.spec.norm[[i]][[j]]@intensity <- 
                spclist.by.sample.spec.norm[[i]][[j]]@intensity * 100
        }
    }
}
#check which samples have all psuedospectra (i.e. none with zero intensity for all)
for (i in 1:length(pd$sample_name)){
   y = which(lapply(spclist.by.sample.spec.norm[[i]], function(x) is.null(x)) == TRUE)
   if (length(y) == 0) {print(i)}
   else {next}
}

#returns samples number 16, 17, 18 which are the metabolite standard mixes --> promising???
comp.samples.test <- c(16,17,18)
names(comp.samples.test) <- pd$sample_name[comp.samples.test]

#initialise dataframe
spectra.comparison <- data.frame(sample.num = as.numeric(),
                                 sample.name = as.character(),
                                 ps.num = as.numeric(),
                                 mona.num = as.numeric(),
                                 mona.name = as.character(),
                                 dot.product = as.numeric())



#compare all spectra
#takes approx 8-9 hours to finish...
for(i in 1:length(comp.samples.test)){
    sample.num = comp.samples.test[i]
    sample.name = names(comp.samples.test)[i]
    sample.spec <- spclist.by.sample.spec.norm[[sample.num]]
    for(j in 1:length(mona.spec)){
        dot.temp <- lapply(sample.spec, function(x)
            compareSpectra(x, mona.spec[[j]], fun ="dot"))
        dot.temp<- unlist(dot.temp)
        temp.df <- data.frame(sample.num = rep(sample.num, length(dot.temp)),
                              sample.name = rep(sample.name, length(dot.temp)),
                              ps.num = 1:length(dot.temp),
                              mona.num = rep(j, length(dot.temp)),
                              mona.name = rep(mona.names[j], length(dot.temp)),
                              dot.product = dot.temp)
        temp.df2 <- temp.df %>% filter(dot.product >= 0.5)
        spectra.comparison <- rbind(spectra.comparison,
                                    temp.df2)
        rm(temp.df, temp.df2, dot.temp)
    }
}


# compare spectra database vs. our data
pd$sample_name

# sample information about ps:
spclist.by.sample[[16]][[27]]

compareSpectra(mona.spec[[14299]], spclist.by.sample.spec.norm[[16]][[22]], fun="dot")

mona.msp[[14299]]


# Hagen tryouts

db_entries<-which(mona.names == "Ribitol")
dot_product<-data.frame(db_entry=numeric(),
                        dot_product=numeric(),
                        sample_number=numeric(),
                        pseudospectrum=numeric())

nr_of_pseudospectra<-264

sample_nrs<-c(16, 17, 18)

for (ps in 1:nr_of_pseudospectra){
    
    for (sample in sample_nrs){
        for (entry in db_entries){
        dot_product_temp<-compareSpectra(mona.spec[[entry]], 
                                     spclist.by.sample.spec.norm[[sample]][[ps]], 
                                     fun="dot")
        dot_product<-rbind(dot_product, c(entry, dot_product_temp, sample, ps))
        colnames(dot_product)<-c("db_entry","dot_product", "sample_number", "pseudospectrum")
        }
    }
}



# go by pseudospectra, find best matching DB entry

db_entries<-mona.names[1:10]
dot_product<-data.frame(db_entry=numeric(),
                        dot_product=numeric(),
                        sample_number=numeric(),
                        pseudospectrum=numeric())

nr_of_pseudospectra<-264

sample_nrs<-c(16, 17, 18)

for (ps in 1:nr_of_pseudospectra){
    
    for (sample in sample_nrs){
        for (entry in db_entries){
            dot_product_temp<-compareSpectra(mona.spec[[entry]], 
                                             spclist.by.sample.spec.norm[[sample]][[ps]], 
                                             fun="dot")
            list_

        }
    }
    dot_product<-rbind(dot_product, c(entry, dot_product_temp, sample, ps))
    colnames(dot_product)<-c("db_entry","dot_product", "sample_number", "pseudospectrum")
}


# Finding ps corresponding to compounds in metmix ---------------------------------------------

rt.spclist <- sapply(spclist, "[", 1, 77) # subset from ps list, row 1, column 77 (contains rt of ps)

# Ribitol -----------------------------------------------------------------

#ribitol: rt 950-1020, mz_ribitol 217-217.4
sample_nr<-16
ribitol_rtmin<-950
ribitol_rtmax<-1020
ribitol_mz<-217

ribitol.sp <- mona.spec[which(mona.names=="Ribitol")][[1]]

index <- which(rt.spclist > ribitol_rtmin & rt.spclist < ribitol_rtmax)
matches<-data.frame(ps_nr=numeric())
for (i in 1:length(index)){
    n = index[i]
    if(ribitol_mz %in% spclist.by.sample[[sample_nr]][[n]][,1]){ # in spclist.by.sample[[sample]][[ps]]
        print(n)
        matches<-rbind(matches, n)
        colnames(matches)<-"ps_nr"
    }
}

par(mfrow=c(length(matches$ps_nr),1))
for (match in matches$ps_nr){
    dot_ribitol = compareSpectra(ribitol.sp,
                             spclist.by.sample.spec.norm[[sample_nr]][[match]],
                             fun = "cor")
    plot(ribitol.sp,
        spclist.by.sample.spec.norm[[sample_nr]][[match]],     
        main = paste0("Ribitol: dot product = ", dot_ribitol))
}

# Result: ribitol is probably ps 22
ribitol_ps<-22


# Glucose -----------------------------------------------------------------

#glucose: rt 1100-1200, mz_glucose 321
sample_nr<-16
glucose_rtmin<-1100
glucose_rtmax<-1200
glucose_mz<-321

glucose.sp <- mona.spec[which(mona.names=="D-Glucose")][[1]]

index <- which(rt.spclist > glucose_rtmin & rt.spclist < glucose_rtmax)
matches<-data.frame(ps_nr=numeric())
for (i in 1:length(index)){
    n = index[i]
    if(glucose_mz %in% spclist.by.sample[[sample_nr]][[n]][,1]){ # in spclist.by.sample[[sample]][[ps]]
        print(n)
        matches<-rbind(matches, n)
        colnames(matches)<-"ps_nr"
    }
}

par(mfrow=c(length(matches$ps_nr),1))
for (match in matches$ps_nr){
    dot_glucose = compareSpectra(glucose.sp,
                                 spclist.by.sample.spec.norm[[sample_nr]][[match]],
                                 fun = "cor")
    plot(glucose.sp,
         spclist.by.sample.spec.norm[[sample_nr]][[match]],     
         main = paste0("Glucose: dot product = ", dot_glucose))
}
# Result: glucose is probably ps 12
glucose_ps<-12


# create function to determine ps_nr
# depends on: 
# - rt.spclist
# - spclist.by.sample
# - spclist.by.sample.spec.norm

# required input:
# - reference spectrum
# - rtmin
# - rtmax

ref_spec_glucose <- mona.spec[which(mona.names=="D-Glucose")][[1]]
find_ps_nr<-function(ref_spec, rt, rtmin, rtmax, mz, sample_nr=16, compound_name){
    
    matches<-data.frame(ps_nr=numeric())
    index <- which(rt.spclist > rtmin & rt.spclist < rtmax)
    matches<-data.frame(ps_nr=numeric())
    for (i in 1:length(index)){
        n = index[i]
        if(mz %in% spclist.by.sample[[sample_nr]][[n]][,1]){ 
            matches<-rbind(matches, n)
            colnames(matches)<-"ps_nr"
        }
    }
    par(mfrow=c(length(matches$ps_nr),1), mai=c(0.6,0.8,0.4,0.2))
    plots<-list()
    if (missing(compound_name)){
        print("no compound name provided")
        for (i in 1:length(matches$ps_nr)){
            dot_prod = compareSpectra(ref_spec,
                                      spclist.by.sample.spec.norm[[sample_nr]][[matches$ps_nr[i]]],
                                      fun = "dot")
            plots[[i]]<-plot(glucose.sp,
                             spclist.by.sample.spec.norm[[sample_nr]][[matches$ps_nr[i]]],     
                             main = paste0("dot product = ", round(dot_prod, 3),
                                           " for pseudospectrum ", matches$ps_nr[i]))
        }
    }
    else {
        print(compound_name)
        for (i in 1:length(matches$ps_nr)){
            dot_prod = compareSpectra(ref_spec,
                                      spclist.by.sample.spec.norm[[sample_nr]][[matches$ps_nr[i]]],
                                      fun = "dot")
            plots[[i]]<-plot(glucose.sp,
                             spclist.by.sample.spec.norm[[sample_nr]][[matches$ps_nr[i]]],     
                             main = paste0(compound_name," dot product = ", round(dot_prod, 3),
                                           " for pseudospectrum ", matches$ps_nr[i]))
        }
    }
    return(list(matches, plots))
}

find_ps_nr(ref_spec = ref_spec_glucose, rtmin=1100, rtmax = 1140, sample_nr = 18,
           mz= glucose_mz,compound_name = "glucose")
