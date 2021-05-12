
# Deconvolution --------------------------------------------------------
#start with xset from GCMS-data-analysis-script_ASSEMBLE.r
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

rt.spclist <- sapply(spclist, "[", 1, 77)

#ribitol: rt 1000-1020, mz_ribitol 217-217.4
ribitol.msp <- read.msp("./msp-mona-singlecompounds/ribitol_MoNA002421.msp")

ribitol.sp <- new("Spectrum1", 
                  mz = ribitol.msp[[1]]$pspectrum[,1], 
                  intensity = ribitol.msp[[1]]$pspectrum[,2])

index <- which(rt.spclist > 1000 & rt.spclist < 1020)
for (i in 1:length(index)){
    n = index[i]
    if(217 %in% spclist.by.sample[[13]][[n]][,1]){
        print(n)
    }
}

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

#create spectrum objects for each db pseudospectrum
mona.spec <- lapply(mona.msp,
                    function(x) new("Spectrum1", 
                                    mz = x$pspectrum[,1], 
                                    intensity = x$pspectrum[,2],
                                    centroided = TRUE))
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




