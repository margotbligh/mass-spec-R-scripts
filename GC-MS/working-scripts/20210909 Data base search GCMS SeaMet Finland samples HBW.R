# Compare pseudospectra to MoNA database ---------------
# based on pseudospectra from GCMS-data-analyis-script_ASSEMBLE
# 

# 1. Preparation ----------------------------------------------------------
library(rebus)
library(scatterplot3d)
time_cat_colors<-c(scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                   scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))
time_cat_4colors<-c("grey",scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                    scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))
time_cat_5colors<-c("white", "grey",scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                    scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))
plant_colors<-scico(2, palette="tofino", begin = 0.6, end=0.9, alpha=0.7)


#read in external database 
#mona.msp <- metaMS::read.msp(file = "MoNA-export-GC-MS_Spectra.msp")
#mona.names <- unlist(sapply(mona.msp, "[", 1))
#create spectrum objects for each db pseudospectrum
          ### CAREFUL: takes super long ###
#mona.spec <- lapply(mona.msp,
#                    function(x) new("Spectrum1", 
#                                    mz = x$pspectrum[,1], 
#                                    intensity = x$pspectrum[,2],
#                                    centroided = TRUE))

setwd("C:/Users/admin/ownCloud/Assemble/Data/SeaMet/GCMS data analysis/results/")
load(file = "RData_metaMStesting_20210511.RData")

#### I think you can directly go to 2. after loading the environment ####



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


# 2. define function for pseudospectra matching ---------------------------

rt.spclist <- sapply(spclist, "[", 1, 77) # subset from ps list, row 1, column 77 (contains rt of ps)

# create function to determine ps_nr
# depends on: 
# - rt.spclist
# - spclist.by.sample
# - spclist.by.sample.spec.norm

# required input:
# - reference spectrum
# - rtmin
# - rtmax


find_ps_nr<-function(ref_spec, rt, rtmin, rtmax, mz, 
                     sample_nr=16, compound_name, print= TRUE){
  
  
  if (missing(rtmin)){
    if (missing(rtmax)){
      if (missing(rt)){
        print("Missing retention time info")
        stop()
      }
      else{
        print(paste("Using ", rt, "+/- 30 sec as retention time window"))
        rtmin<-rt-60
        rtmax<-rt+60
      }
    }
    else {
      if (missing(rt)){
        print("Missing retention time info")
        stop()
      }
      else{
        print(paste("Using ", rt, "-30 sec up to ", rtmax, " as retention time window"))
        rtmin<-rt-60
      }
    }
  }
  else{
    if (missing(rtmax)){
      if (missing(rt)){
        print("Missing retention time info")
        stop()
      }
      else{
        print(paste("Using ", rtmin, " up to ", rt ," +30 sec as retention time window"))
        rtmax<-rt+60
      }
    }
  }
  
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
  if (length(matches$ps_nr)==0){
    print(paste("No match found for ", compound_name))
  }
  while(length(matches$ps_nr)>=1){
  
  plots<-list()
  if (missing(compound_name)){
    print("no compound name provided")
    for (i in 1:length(matches$ps_nr)){
      for (j in 1:length(ref_spec)){
      dot_prod = compareSpectra(ref_spec[[j]],
                                spclist.by.sample.spec.norm[[sample_nr]][[matches$ps_nr[i]]],
                                fun = "dot")
      if(print==T){
      plots[[i]]<-plot(ref_spec[[j]],
                       spclist.by.sample.spec.norm[[sample_nr]][[matches$ps_nr[i]]],     
                       main = paste0("dotpr=", round(dot_prod, 3),
                                     " for ps ", matches$ps_nr[i], " ref spec ", j))
      }
      }
    }
  }
  else {
    print(compound_name)
    for (i in 1:length(matches$ps_nr)){
      for (j in 1:length(ref_spec)){
      dot_prod = compareSpectra(ref_spec[[j]],
                                spclist.by.sample.spec.norm[[sample_nr]][[matches$ps_nr[i]]],
                                fun = "dot")
      if(print==T){
      plots[[i]]<-plot(ref_spec[[j]],
                       spclist.by.sample.spec.norm[[sample_nr]][[matches$ps_nr[i]]],     
                       main = paste0(compound_name," dotpr=", round(dot_prod, 3),
                                     " for ps ", matches$ps_nr[i], " ref spec ", j))
      }
      }
    }
  }
  if(print==TRUE){
  par(mfrow=c(3,2), mai=c(0.6,0.8,0.4,0.2))
    return(list(matches, plots))
  }
  matches<-NULL
}
}




# Find pseudospectra matching standards -----------------------------------
standard_ps_table<-data.frame(Compound = character(),
                              Pseudospectrum = numeric(),
                              stringsAsFactors = F)

ref_spec_ribitol<-mona.spec[which(mona.names=="Ribitol")]
ribitol_rt<-1008
ribitol_mz<-217

find_ps_nr(ref_spec = ref_spec_ribitol, rt = ribitol_rt, mz = ribitol_mz, compound_name = "Ribitol")

standard_ps_table<-rbind(standard_ps_table, c("Ribitol", 22))
colnames(standard_ps_table)<-c("Compound", "Pseudospectrum")


ref_spec_cholestane<-mona.spec[which(mona.names=="CHOLESTANE")]
cholestane_rt<-1622
cholestane_mz<-373

find_ps_nr(ref_spec = ref_spec_cholestane, rt = cholestane_rt, mz = cholestane_mz, compound_name = "Cholestane")

standard_ps_table<-rbind(standard_ps_table, c("Cholestane", 17))
colnames(standard_ps_table)<-c("Compound", "Pseudospectrum")



ref_spec_mannitol<-mona.spec[which(mona.names=="mannitol")]
mannitol_rt<-1140
mannitol_mz<-319

find_ps_nr(ref_spec = ref_spec_mannitol, rt = mannitol_rt, mz = mannitol_mz, compound_name = "mannitol")

# difficult to say if ps 12 or 16 corresponds to mannitol
# dotproduct and commons 0.63 and 29 for ps 12 and 0.57 and 59 for ps 16
# rt of ps 12 is 1140, fits ref value best
standard_ps_table<-rbind(standard_ps_table, c("mannitol", 12))
colnames(standard_ps_table)<-c("Compound", "Pseudospectrum")

ref_spec_mannose<-mona.spec[which(mona.names=="D-Mannose")]
mannose_rt<-1126
mannose_mz<-319 # from list should look for 203

find_ps_nr(ref_spec = ref_spec_mannose, rt = mannose_rt, mz = mannose_mz, compound_name = "Mannose")
# good hits for ps 12, 16 (with mz319) and 54 (with mz203)
# ps 16 fits best in terms of rt
standard_ps_table<-rbind(standard_ps_table, c("Mannose", 16))
colnames(standard_ps_table)<-c("Compound", "Pseudospectrum")



ref_spec_glucose <- mona.spec[which(mona.names=="D-Glucose")]
glucose_rt<-1126
glucose_mz<-203

find_ps_nr(ref_spec = ref_spec_glucose, rt=1120, sample_nr = 18,
           mz= glucose_mz,compound_name = "glucose")



# look for sugars not in standard -----------------------------------------

rt.spclist
glucosamine<-mona.spec[which(mona.names=="D-Glucosamine")]

glucosamine_table<-data.frame(sample<-character(),
                              pseudospectrum<-numeric(),
                              retention_time<-numeric(),
                              dot_product<-numeric(),
                              stringsAsFactors = F)
for (sample_nr in 1:75){
  for (i in 1:264){
    if(!is.null(spclist.by.sample.spec.norm[[sample_nr]][[i]])){
    product_temp<-compareSpectra(glucosamine[[1]],
               spclist.by.sample.spec.norm[[sample_nr]][[i]],
               fun = "dot")
    pseudospectrum_temp<-i
    retention_time_temp<-rt.spclist[i]
    sample_temp<-colnames(spclist[[1]])[as.numeric(sample_nr)+1]
    glucosamine_table<-rbind(glucosamine_table, c(sample_temp, pseudospectrum_temp,
                                                retention_time_temp, product_temp))
    colnames(glucosamine_table)<-c("sample", "pseudospectrum", 
                                   "retention_time", "dot_product")
    }
  }

}

glucosamine_hits_across_samples<-glucosamine_table %>%
  dplyr::group_by(pseudospectrum, retention_time) %>%
  dplyr::summarise(mean_dot_product<-mean(as.numeric(dot_product)), 
                   sdt_dev_dot_product<-sd(as.numeric(dot_product)), 
                   n=n())

glucosamine<-mona.spec[which(mona.names=="D-Glucosamine")]
find_ps_nr(ref_spec=glucosamine, rt = 1140, mz = 320, sample_nr = 1)



# 3. STANDARD COMPOUNDS TARGETED ------------------------------------------


# match compounds from table with standards ----------------------------------------------

# import information on standards

standard_compounds<-read.csv("C:/Users/admin/ownCloud/Assemble/Data/SeaMet/GCMS_extendedcompoundlist_forR.csv",
                             sep=",", h=T, stringsAsFactors = F)
standard_compounds<-standard_compounds[complete.cases(standard_compounds[,2:4]),]
# create reference db subset


for (h in 1:length(standard_compounds[,1])){
  compound_name_temp<-standard_compounds[h,1]
  ref_spec_temp <- mona.spec[which(mona.names==compound_name_temp)]
  if (length(ref_spec_temp)==0){
    print(paste("No reference spectrum found for ", compound_name_temp))
    next
  }
  if (is.na(standard_compounds[h,3])){
    rt_temp<-standard_compounds[h,2]
  }
  else{
    rt_temp<-standard_compounds[h,3]
  }
  mz_temp<-standard_compounds[h,4]
  find_ps_nr(ref_spec = ref_spec_temp, rtmin = rt_temp-100, rtmax = rt_temp+100, 
             mz = mz_temp,
             compound_name = compound_name_temp, sample_nr = 16,
             print=T)
}


# Compound and pseudospectra matching manually ----------------------------

standard_compounds$pseudospectrum<-NA
standard_compounds$pseudospectrum[standard_compounds$Compound=="valine"]<-NA
standard_compounds$pseudospectrum[standard_compounds$Compound=="Thymine"]<-NA
standard_compounds$pseudospectrum[standard_compounds$Compound=="threonine"]<-32
standard_compounds$pseudospectrum[standard_compounds$Compound=="sucrose"]<-23 # okay matches for 4, 23, 63, rts off
standard_compounds$pseudospectrum[standard_compounds$Compound=="succinic acid"]<-NA
standard_compounds$pseudospectrum[standard_compounds$Compound=="serine"]<-14
standard_compounds$pseudospectrum[standard_compounds$Compound=="ribose"]<-27
standard_compounds$pseudospectrum[standard_compounds$Compound=="ribitol"]<-22
standard_compounds$pseudospectrum[standard_compounds$Compound=="proline"]<-25
standard_compounds$pseudospectrum[standard_compounds$Compound=="Phenylalanine"]<-NA
standard_compounds$pseudospectrum[standard_compounds$Compound=="Myo-Inositol"]<-21
standard_compounds$pseudospectrum[standard_compounds$Compound=="Mannitol"]<-90
standard_compounds$pseudospectrum[standard_compounds$Compound=="malic acid"]<-NA
standard_compounds$pseudospectrum[standard_compounds$Compound=="lysine"]<-NA #30?
standard_compounds$pseudospectrum[standard_compounds$Compound=="lactic acid"]<-NA
standard_compounds$pseudospectrum[standard_compounds$Compound=="isoleucine"]<-20
standard_compounds$pseudospectrum[standard_compounds$Compound=="glutamic acid"]<-NA
standard_compounds$pseudospectrum[standard_compounds$Compound=="fumaric acid"]<-NA
standard_compounds$pseudospectrum[standard_compounds$Compound=="alanine"]<-NA
standard_compounds$pseudospectrum[standard_compounds$Compound=="D-Galactose"]<-16
standard_compounds$pseudospectrum[standard_compounds$Compound=="trehalose"]<-4 # okay matches for 4, 23 
standard_compounds$pseudospectrum[standard_compounds$Compound=="D-Mannose"]<-54
standard_compounds$pseudospectrum[standard_compounds$Compound=="D-Glucose"]<-12
standard_compounds$pseudospectrum[standard_compounds$Compound=="fructose"]<-35
standard_compounds$pseudospectrum[standard_compounds$Compound=="CHOLESTANE"]<-17
standard_compounds$pseudospectrum[standard_compounds$Compound=="fucose"]<-30


rt.spclist.table<-as.data.frame(rt.spclist)
rt.spclist.table$ps_nr<-seq(1,264,1)

standard_compounds$rt_pseudospectrum<-NA
for (i in 1:length(standard_compounds$pseudospectrum)){
  ps_nr_temp<-standard_compounds$pseudospectrum[i]
  standard_compounds$rt_pseudospectrum[i]<-rt.spclist.table$rt.spclist[rt.spclist.table$ps_nr==ps_nr_temp]
}


rt.spclist[12]
rt.spclist.table<-as.data.frame(rt.spclist)
rt.spclist.table$ps_nr<-seq(1,264,1)

# resolving them sugars ######
#galactose
galactose_rt<-rt.spclist.table[rt.spclist.table$rt.spclist<=1250 & 
                   rt.spclist.table$rt.spclist>=1000,]
ref_spec<-mona.spec[which(mona.names=="D-Galactose")]
for (i in galactose_rt$ps_nr){
 temp_common<-compareSpectra(ref_spec[[2]],
               spclist.by.sample.spec.norm[[16]][[i]],
               fun = "common")
 temp_dot<-compareSpectra(ref_spec[[2]],
                             spclist.by.sample.spec.norm[[16]][[i]],
                             fun = "dot")
 print(paste(i, "common:",temp_common, "dotprod:", temp_dot, "rt:", galactose_rt$rt.spclist[galactose_rt$ps_nr==i]))
}
ps16_galactose<-data.frame(list(spclist.by.sample.spec.norm[[16]][[16]]@mz, 
                                spclist.by.sample.spec.norm[[16]][[16]]@intensity))
colnames(ps16_galactose)<-c("mz", "int")

mannitol_rt<-rt.spclist.table[rt.spclist.table$rt.spclist<=1200 & 
                                 rt.spclist.table$rt.spclist>=1120,]
ref_spec<-mona.spec[which(mona.names=="mannitol")]
for (i in mannitol_rt$ps_nr){
  temp_common<-compareSpectra(ref_spec[[1]],
                              spclist.by.sample.spec.norm[[18]][[i]],
                              fun = "common")
  temp_dot<-compareSpectra(ref_spec[[1]],
                           spclist.by.sample.spec.norm[[18]][[i]],
                           fun = "dot")
  print(paste(i, "common:",temp_common, "dotprod:", temp_dot, "rt:", mannitol_rt$rt.spclist[mannitol_rt$ps_nr==i]))
}

fucose_rt<-rt.spclist.table[rt.spclist.table$rt.spclist<=1060 & 
                                rt.spclist.table$rt.spclist>=980,]
ref_spec<-mona.spec[which(mona.names=="fucose")]
for (i in fucose_rt$ps_nr){
  temp_common<-compareSpectra(ref_spec[[2]],
                              spclist.by.sample.spec.norm[[3]][[i]],
                              fun = "common")
  temp_dot<-compareSpectra(ref_spec[[2]],
                           spclist.by.sample.spec.norm[[3]][[i]],
                           fun = "dot")
  print(paste(i, "common:",temp_common, "dotprod:", temp_dot, "rt:", fucose_rt$rt.spclist[fucose_rt$ps_nr==i]))
}

# is ps 30 fucose?
find_ps_nr(ref_spec = ref_spec, rt=1036, mz=277, sample_nr = 6, compound_name = "fucose")
# could be
# Calculate table with sample data ----------------------------------------

# tryout
plots<-NULL
for (compound in standard_compounds$Compound){
  ref_spec_temp<-mona.spec[which(mona.names==compound)]
  pseudospectrum_temp<-standard_compounds$pseudospectra[standard_compounds$Compound==compound]
  for (j in 1:length(ref_spec_temp)){
    dot_prod = compareSpectra(ref_spec_temp[[j]],
                              spclist.by.sample.spec.norm[[i]][[pseudospectrum_temp]],
                              fun = "dot")
    plots[[i]]<-plot(ref_spec_temp[[j]],
                     spclist.by.sample.spec.norm[[2]][[pseudospectrum_temp]],     
                     main = paste0(compound," dotpr=", round(dot_prod, 3),
                                   " for ps ", pseudospectrum_temp, " ref spec ", j))
  }
}

# all samples

# try to fix by changing quant ions
standard_compounds$mz[standard_compounds$Compound=="D-Glucose"]<-205
standard_compounds$mz[standard_compounds$Compound=="D-Mannose"]<-203
standard_compounds$mz[standard_compounds$Compound=="ribose"]<-205

#
#remove double entries
standard_compounds_cleaned<-standard_compounds[-17,]
standard_compounds_cleaned<-standard_compounds_cleaned[-18,]

sample_names<-colnames(spclist[[1]])[2:76]
sample_compounds<-data.frame(Compound=character(),
                             Pseudospectrum=numeric(),
                             Retention_time=numeric(),
                             Sample_name=character(),
                             Dot_product=numeric(),
                             Quantion_mz=numeric(),
                             Quantion_intensity=numeric(),
                             stringsAsFactors = F)
for (i in 1:length(sample_names)){
  sample_name_temp<-sample_names[i]

  for (compound in standard_compounds_cleaned$Compound){
    if (is.na(standard_compounds_cleaned$pseudospectrum[standard_compounds_cleaned$Compound==compound])==TRUE){
      next
    }
    ref_spec_temp<-mona.spec[which(mona.names==compound)][[1]]
    pseudospectrum_temp<-standard_compounds_cleaned$pseudospectrum[standard_compounds_cleaned$Compound==compound]
    retention_time_temp<-standard_compounds_cleaned$rt_pseudospectrum[standard_compounds_cleaned$Compound==compound]
    dot_prod_temp = compareSpectra(ref_spec_temp,
                            spclist.by.sample.spec.norm[[i]][[pseudospectrum_temp]],
                            fun = "dot")
    quantion_mz_temp<-standard_compounds_cleaned$mz[standard_compounds_cleaned$Compound==compound]
    quantion_intensity_temp<-spclist.by.sample.spec[[i]][[pseudospectrum_temp]]@intensity[spclist.by.sample.spec[[i]][[pseudospectrum_temp]]@mz==quantion_mz_temp]
    sample_compounds<-rbind(sample_compounds, c(compound, pseudospectrum_temp, retention_time_temp,
                                                sample_name_temp, dot_prod_temp, quantion_mz_temp,
                                                quantion_intensity_temp))
    colnames(sample_compounds)<-c("Compound", "Pseudospectrum", "Retention_time",
                                  "Sample_name", "Dot_product", "Quantion_mz",
                                  "Quantion_intensity")
    }
}

#normalize for intensity of ribitol
sample_compounds$ribitol_norm<-NA
for (sample in sample_compounds$Sample_name){
  ribitol_temp<-
    sample_compounds$Quantion_intensity[sample_compounds$Sample_name==sample & sample_compounds$Compound=="ribitol"]
  sample_compounds$ribitol_norm[sample_compounds$Sample_name==sample]<-
    as.numeric(sample_compounds$Quantion_intensity[sample_compounds$Sample_name==sample])/as.numeric(ribitol_temp)
}

#### extract metadata from sample names #### 

meta_temp=data.frame(plant=character(),
                     water_body=character(),
                     timepoint=character(),
                     treatment=character(), 
                     time_treat=character())

water_body_temp=character()
plant_temp=character()
timepoint_temp=character()
treatmeant_temp=character()
time_treat_temp=character()

for (i in 1:length(sample_compounds$Compound)){
  if (grepl(sample_compounds$Sample_name[i], pattern = "^ASW_")){
    if (grepl(sample_compounds$Sample_name[i], pattern = "^ASW_met")){
      plant_temp<-"none"
      water_body_temp<-"Artificial seawater"
      timepoint_temp<-"none"
      treatment_temp<-"metabolite mix"
      time_treat_temp<-"ASW metmix"
    }
    if (grepl(sample_compounds$Sample_name[i], pattern = "^ASW_blank")){
      plant_temp<-"none"
      water_body_temp<-"Artificial seawater"
      timepoint_temp<-"none"
      treatment_temp<-"blank"
      time_treat_temp<-"ASW blank"
    }
  }
  if (grepl(sample_compounds$Sample_name[i], pattern = "^BLK_")){
    plant_temp<-"none"
    water_body_temp<-"Solvent"
    timepoint_temp<-"none"
    treatment_temp<-"blank"
    time_treat_temp<-"solvent blank"
  } 
  # start loop for benthic chamber samples
  if (grepl(sample_compounds$Sample_name[i], pattern = "^[A-Z]_")){
        name_parts_temp<-unlist(str_split(sample_compounds$Sample_name[i], "_"))
    if (name_parts_temp[1]=="A"){
      plant_temp<-"Fucus"  
    }
    if (name_parts_temp[1]=="S"){
      plant_temp<-"Zostera"
    }    
    if (grepl(name_parts_temp[2], pattern = "^WC")){
      water_body_temp<-"water column"
    }
    if (grepl(name_parts_temp[2], pattern = "^PW")){
      water_body_temp<-"pore water"
    }      
    if (grepl(name_parts_temp[2], pattern="I$|l$")){  #there was a typo in the file name
      timepoint_temp<-"initial"
      time_treat_temp<-"initial"
    }
    
    if (grepl(name_parts_temp[2], pattern="E$")){
      timepoint_temp<-"end"
      name_subparts_temp<-unlist(strsplit(name_parts_temp[2], split=""))
      if(name_subparts_temp[3]<=3){
        treatment_temp="light"
        time_treat_temp<-"end_light"
      }
      if(name_subparts_temp[3]>=4){
        treatment_temp="dark"
        time_treat_temp<-"end_dark"
      }
      
    }
  }  # end loop for benthic chamber samples
    
  meta_temp<-rbind(meta_temp, c(as.character(plant_temp), as.character(water_body_temp),
                                as.character(timepoint_temp), 
                                as.character(treatment_temp), as.character(time_treat_temp)))
  colnames(meta_temp)<-c("plant","water_body", "timepoint", "treatment", "time_treat")
}

sample_compounds$plant<-meta_temp$plant
sample_compounds$water_body<-meta_temp$water_body
sample_compounds$timepoint<-meta_temp$timepoint
sample_compounds$treatment<-meta_temp$treatment
sample_compounds$time_treat<-meta_temp$time_treat


# explorative plotting

time_cat_colors<-c(scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.8),
                   scico(1, palette = "davos",begin=0.1,end=0.2, alpha=0.8),
                   scico(1, palette = "lajolla",begin=0.2,end=0.3, alpha=0.8))
sample_compounds_wo_int_std<-sample_compounds[which(sample_compounds$Compound %in% c(
  "fructose", "D-Glucose", "Mannitol", "D-Mannose", "ribose", "D-Galactose","isoleucine",
  "Myo-Inositol", "proline", "serine",
  "sucrose",    "trehalose", "fucose")
),]
sample_compounds_wo_int_std$water_body_plant<-paste(sample_compounds_wo_int_std$water_body,
                                                    sample_compounds_wo_int_std$plant,
                                                    sep = "_")
sample_compounds_for_plot<-sample_compounds_wo_int_std[sample_compounds_wo_int_std$water_body_plant!="pore water_Zostera",]
sample_compounds_for_plot$time_treat<-factor(sample_compounds_for_plot$time_treat, 
                                             levels = c("none" ,"initial", "end_dark", "end_light"))
sample_compounds_for_plot$Compound<-factor(sample_compounds_for_plot$Compound,
                                           levels =   c("fructose", "fucose", "D-Galactose", "D-Glucose","isoleucine", "Mannitol", "D-Mannose", 
                                           "Myo-Inositol", "proline", "ribose", "serine",
                                           "sucrose", "trehalose"),
                                           labels= c("fructose", "fucose", "galactose", "glucose",
                                                     "isoleucine", "mannitol", "mannose", 
                                                     "myo-inositol", "proline", "ribose", "serine",
                                                     "sucrose", "trehalose"))
sample_compounds_for_plot$time_treat<-factor(sample_compounds_for_plot$time_treat,
                                             levels=c("initial", "end_dark", "end_light"),
                                             labels=c("Seawater", "Dark inc.", "Light inc."))
# FIGURE compounds in samples ----------------------------------------



fig_samples<-ggplot(data=sample_compounds_for_plot[sample_compounds_for_plot$timepoint!="none",])+
  geom_bar(aes(x=Compound, y=sqrt(as.numeric(ribitol_norm)), fill=plant),
           position = position_dodge(), stat = "summary", colour="black")+
  geom_point(aes(x=Compound, y=sqrt(as.numeric(ribitol_norm)), fill=plant), 
             position = position_dodge(.9), pch=21, size=2)+
  scale_y_continuous(name=expression(sqrt(Ribitol-normalized~intensity~(a.u.))))+
  scale_fill_manual(name=NULL,  breaks=c("Fucus", "Zostera"),
                    labels=c(expression(italic(Fucus~vesiculosus)), 
                             expression(italic(Zostera~marina))),
                    values = plant_colors)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.y = element_text(size=font_size, family = font_family, color = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.y = element_rect(fill = c(time_cat_colors[1], 
                                                   time_cat_colors[2], time_cat_colors[3])),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "#262626"))+
  facet_wrap(~time_treat, nrow=3, strip.position = "right")
h <- ggplot_gtable(ggplot_build(fig_samples))
stripr <- which(grepl('strip-r', h$layout$name))
fills <- time_cat_colors
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', h$grobs[[i]]$grobs[[1]]$childrenOrder))
  h$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(h)

tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210902_fig_gc-ms.tiff", 
     height=14, width = 28, units = "cm", res=300)
grid.draw(h)
dev.off()
pdf("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210902_fig_gc-ms.pdf", 
    height=5, width = 12)
grid.draw(h)
dev.off()

# plot blanks and standards -----------------------------------------------

fig_blanks_stds<-ggplot(data=sample_compounds[sample_compounds$plant=="Artificial seawater",])+
  geom_bar(aes(x=Compound, y=as.numeric(ribitol_norm), fill=time_treat),
           position = position_dodge(), stat = "summary", colour="black")+
  geom_point(aes(x=Compound, y=as.numeric(ribitol_norm), fill=time_treat), 
             position = position_dodge(.9), pch=21, size=2)+
  scale_y_continuous(name="Ribitol-normalized intensity (a.u.)", limits = c(0,0.01))+
  #  scale_y_continuous(name=expression(paste(Concentration~(microg~L^-1))))+
#  scale_fill_manual(name=NULL,  breaks=c("initial", "end_dark", "end_light"),
#                    labels=c("before incubation",
#                             "incubation in dark",
#                             "incubation in light"),
#                    values = time_cat_colors)+
  facet_grid(rows="treatment")+
  theme_bw()

plot_grid(fig_blanks_stds, fig_samples)


#plot porewater vs. watercolumn

fig_samples_pw_vs_wc<-ggplot(data=sample_compounds[sample_compounds$plant=="Zostera" ,])+
  geom_bar(aes(x=Compound, y=as.numeric(ribitol_norm), fill=time_treat),
           position = position_dodge(), stat = "summary", colour="black")+
  geom_point(aes(x=Compound, y=as.numeric(ribitol_norm), fill=time_treat), 
             position = position_dodge(.9), pch=21, size=2)+
  scale_y_continuous(name="Ribitol-normalized intensity (a.u.)", limits = c(0,0.005))+
  #  scale_y_continuous(name=expression(paste(Concentration~(microg~L^-1))))+
  scale_fill_manual(name=NULL,  breaks=c("initial", "end_dark", "end_light"),
                    labels=c("before incubation",
                             "incubation in dark",
                             "incubation in light"),
                    values = time_cat_colors)+
  facet_grid(rows="water_body")+
  theme_bw()


# 4. ACCUMULATION UNTARGETED------------------------------------------------------------



# preselection of peaks based on treatment --------------------------------

# selection pattern:
# before < dark < light

# a. calculate means first, select using means

# abusing meta_temp to get sample metadata:

meta_abuse<-meta_temp[seq(1,675,9),]


# melting xsannotate peaklist to combine with sample info

samplenames_dotdotdot<-gsub(" - ", "...", samplenames)
samplenames_dotdotdot<-gsub(".mzML", "", samplenames_dotdotdot)


peaklist<-getPeaklist(an, intval="into")

allpks_melted<-melt(peaklist, id.vars = c("mz", "rt", "pcgroup"),
                    measure.vars = samplenames_dotdotdot)

meta_rep<-data.frame(plant=character(),
                     water_body=character(),
                     timepoint=character(),
                     treatment=character(),
                     time_treat=character(),
                     stringsAsFactors = F)
nr_peaks<-8960

for (i in 1:length(meta_abuse[,1])){
  start_temp<-nr_peaks*(i-1)+1
  end_temp<-nr_peaks*i
  meta_rep[start_temp:end_temp,]<-meta_abuse[i,]
}

allpks_winfo<-data.frame(allpks_melted, meta_rep)
allpks_winfo$FT_nr<-rep(seq(1,8960,1),75)

meanpks_winfo<-allpks_winfo %>%
  dplyr::group_by(FT_nr, mz, rt, pcgroup, plant, time_treat) %>%
  dplyr::summarise(mean_int=mean(value))


# let's try filtering with 2xASW blank < Fucus initial < Fucus end_light

candidates1<-data.frame(FT_nr=numeric(),
                        mz=numeric(),
                        rt=numeric(),
                        pseudospec=numeric(),
                        mean_asw_blanks=numeric(),
                        mean_Fv_initial=numeric(),
                        mean_Fv_end_light=numeric(),
                        mean_Fv_end_dark=numeric(),
                        stringsAsFactors = F)

for (ft in 1:max(meanpks_winfo$FT_nr)){
  if(   is.numeric(meanpks_winfo$mean_int[meanpks_winfo$FT_nr==ft & 
                                            meanpks_winfo$time_treat=="initial" &
                                            meanpks_winfo$plant=="Fucus"]) &&
        is.numeric(meanpks_winfo$mean_int[meanpks_winfo$FT_nr==ft & 
                                                                  meanpks_winfo$time_treat=="end_light" &
                                                                  meanpks_winfo$plant=="Fucus"]) &&
        is.numeric(meanpks_winfo$mean_int[meanpks_winfo$FT_nr==ft & 
                                                                 meanpks_winfo$time_treat=="end_dark" &
                                                                 meanpks_winfo$plant=="Fucus"])
        ){
  mz_temp<-meanpks_winfo$mz[meanpks_winfo$FT_nr==ft][1]
  rt_temp<-meanpks_winfo$rt[meanpks_winfo$FT_nr==ft][1]
  pseudospec_temp<-meanpks_winfo$pcgroup[meanpks_winfo$FT_nr==ft][1]
  mean_asw_blanks_temp<-meanpks_winfo$mean_int[meanpks_winfo$FT_nr==ft & 
                                                 meanpks_winfo$time_treat=="ASW blank"]
  if(is.na(mean_asw_blanks_temp)){
    mean_asw_blanks_temp<-0
  }
  mean_Fv_initial_temp<-meanpks_winfo$mean_int[meanpks_winfo$FT_nr==ft & 
                                                 meanpks_winfo$time_treat=="initial" &
                                                 meanpks_winfo$plant=="Fucus"]
  if(is.na(mean_Fv_initial_temp)){
    mean_Fv_initial_temp<-0
  }
  mean_Fv_end_light_temp<-meanpks_winfo$mean_int[meanpks_winfo$FT_nr==ft & 
                                                   meanpks_winfo$time_treat=="end_light" &
                                                   meanpks_winfo$plant=="Fucus"]
  mean_Fv_end_dark_temp<-meanpks_winfo$mean_int[meanpks_winfo$FT_nr==ft & 
                                                 meanpks_winfo$time_treat=="end_dark" &
                                                 meanpks_winfo$plant=="Fucus"]
  if(is.na(mean_Fv_end_dark_temp)){
    mean_Fv_end_dark_temp<-0
  }
  if(!is.na(mean_Fv_end_light_temp) &&
     mean_asw_blanks_temp*2 < mean_Fv_end_light_temp){
    if(mean_Fv_initial_temp < mean_Fv_end_light_temp){
      candidates1<-rbind(candidates1, c(ft, mz_temp, rt_temp,pseudospec_temp, mean_asw_blanks_temp,
                                        mean_Fv_initial_temp, mean_Fv_end_light_temp,
                                        mean_Fv_end_dark_temp))
      colnames(candidates1)<-c("FT_nr", "mz", "rt", "pseudospec_nr", "mean_asw_blanks", "mean_Fv_initial",
                               "mean_Fv_end_light", "mean_Fv_end_dark")
    }
  }
  }
}

pseudospectra_count<-dplyr::group_by(candidates1, pseudospec_nr)%>%
 dplyr::summarise(n=n())


candidates1.1<-candidates1[which(candidates1$pseudospec_nr %in% 
                                   pseudospectra_count$pseudospec_nr[pseudospectra_count$n>=10]),]

# b. do not calculate means first -----------------------------------------

# candidates2 contains features for which the intensity in the lowest end_light replicate
# exceeds the intensity in the highest asw blank replicate
# and the lowest end_light replicate exceeds the mean of initials

# this resulted in 522 features from 154
# pseudospectra_count2 contains the number of obtained features per pseudospectrum

candidates2<-data.frame(FT_nr=numeric(),
                        mz=numeric(),
                        rt=numeric(),
                        pseudospec=numeric(),
                        mean_asw_blanks=numeric(),
                        sd_asw_blanks=numeric(),
                        mean_Fv_initial=numeric(),
                        sd_Fv_initial=numeric(),
                        mean_Fv_end_light=numeric(),
                        sd_Fv_end_light=numeric(),
                        mean_Fv_end_dark=numeric(),
                        sd_Fvend_dark=numeric(),
                        stringsAsFactors = F)
allpks_winfo_noNA<-allpks_winfo
allpks_winfo_noNA[is.na(allpks_winfo_noNA)]<-0

for (i in 1:max(allpks_winfo_noNA$FT_nr)){
  mz_temp<-allpks_winfo_noNA$mz[allpks_winfo_noNA$FT_nr==i][1]
  rt_temp<-allpks_winfo_noNA$rt[allpks_winfo_noNA$FT_nr==i][1]
  pseudospec_temp<-allpks_winfo_noNA$pcgroup[allpks_winfo_noNA$FT_nr==i][1]
  
  min_light<-min(allpks_winfo_noNA$value[allpks_winfo_noNA$FT_nr==i & 
                                           allpks_winfo_noNA$plant =="Fucus" &
                                           allpks_winfo_noNA$time_treat=="end_light"])
  max_blank<-max(allpks_winfo_noNA$value[allpks_winfo_noNA$FT_nr==i & 
                                           allpks_winfo_noNA$plant =="none" &
                                           allpks_winfo_noNA$time_treat=="ASW blank"])
  if (min_light>max_blank){
    mean_initial<-mean(allpks_winfo_noNA$value[allpks_winfo_noNA$FT_nr==i & 
                                             allpks_winfo_noNA$plant =="Fucus" &
                                             allpks_winfo_noNA$time_treat=="initial"])
    if (min_light>mean_initial){
      mean_asw_blanks<-mean(allpks_winfo_noNA$value[allpks_winfo_noNA$FT_nr==i & 
                                                          allpks_winfo_noNA$plant =="none" &
                                                          allpks_winfo_noNA$time_treat=="ASW blank"])
      sd_asw_blanks<-sd(allpks_winfo_noNA$value[allpks_winfo_noNA$FT_nr==i & 
                                                           allpks_winfo_noNA$plant =="none" &
                                                           allpks_winfo_noNA$time_treat=="ASW blank"])
      mean_Fv_initial<-mean(allpks_winfo_noNA$value[allpks_winfo_noNA$FT_nr==i & 
                                                     allpks_winfo_noNA$plant =="Fucus" &
                                                     allpks_winfo_noNA$time_treat=="initial"])
      sd_Fv_initial<-sd(allpks_winfo_noNA$value[allpks_winfo_noNA$FT_nr==i & 
                                                      allpks_winfo_noNA$plant =="Fucus" &
                                                      allpks_winfo_noNA$time_treat=="initial"])
      mean_Fv_end_light<-mean(allpks_winfo_noNA$value[allpks_winfo_noNA$FT_nr==i & 
                                               allpks_winfo_noNA$plant =="Fucus" &
                                               allpks_winfo_noNA$time_treat=="end_light"])
      sd_Fv_end_light<-sd(allpks_winfo_noNA$value[allpks_winfo_noNA$FT_nr==i & 
                                                        allpks_winfo_noNA$plant =="Fucus" &
                                                        allpks_winfo_noNA$time_treat=="end_light"])
      mean_Fv_end_dark<-mean(allpks_winfo_noNA$value[allpks_winfo_noNA$FT_nr==i & 
                                                        allpks_winfo_noNA$plant =="Fucus" &
                                                        allpks_winfo_noNA$time_treat=="end_dark"])
      sd_Fv_end_dark<-sd(allpks_winfo_noNA$value[allpks_winfo_noNA$FT_nr==i & 
                                                    allpks_winfo_noNA$plant =="Fucus" &
                                                    allpks_winfo_noNA$time_treat=="end_dark"])
      
      candidates2<-rbind(candidates2, c(i, mz_temp, rt_temp, pseudospec_temp, mean_asw_blanks, sd_asw_blanks,
                                        mean_Fv_initial, sd_Fv_initial, mean_Fv_end_light, sd_Fv_end_light,
                                        mean_Fv_end_dark, sd_Fv_end_dark))
      colnames(candidates2)<-c("FT_nr", "mz", "rt", "pseudospec", "mean_asw_blanks", "sd_asw_blanks",
                               "mean_Fv_initial", "sd_Fv_initial", "mean_Fv_end_light", "sd_Fv_end_light",
                               "mean_Fv_end_dark", "sd_Fv_end_dark")
    }
  }
}

pseudospectra_count2<-dplyr::group_by(candidates2, pseudospec)%>%
  dplyr::summarise(n=n())


# SAVEPOINT: workspace image saved ----------------------------------------

save.image(file = "Workspace_GCMS_DBsearch.RData")
load("Workspace_GCMS_DBsearch.RData")


# attempt annotation by comparison to mona.spec -> fail ------------------

hits_ps32<-data.frame(ref_compound=character(),
                      dot_product=numeric(),
                      common=numeric(),
                      stringsAsFactors = F)
for (j in 1:length(mona.spec)){
  dot_prod = compareSpectra(mona.spec[[j]],
                            spclist.by.sample.spec.norm[[5]][[32]],
                            fun = "dotproduct")
  common = compareSpectra(mona.spec[[j]],
                            spclist.by.sample.spec.norm[[5]][[32]],
                            fun = "common", tolerance=0.1)
  name_temp<-mona.names[[j]]
  dp_temp<-dot_prod
  hits_ps32<-rbind(hits_ps32, c(name_temp, as.numeric(dp_temp), as.numeric(common)))
  colnames(hits_ps32)<-c("ref_compound", "dot_product", "common")
}

hits_ps32$dot_product<-as.numeric(hits_ps32$dot_product)
hits_ps32$common<-as.numeric(hits_ps32$common)


# list of m/z values of interest: 203, 205, 217, 307, 319, 321, 361, 377
ions_of_interest<-c(203, 205, 217, 307, 319, 321, 361, 377)
for (i in ions_of_interest){
  hits<-unlist(grep(pattern=i, candidates2$mz))
  if (length(hits)!=0){
    for (j in 1:length(hits)){
      print(candidates2$pseudospec[hits[j]])
    }
  }
}


# Tentative annotations ---------------------------------------------------


# seeking annotation for ps nr 13, 26, 74, 78, 105

rt.spclist
mannitol<-mona.spec[which(mona.names=="Mannitol")]

mannitol_table<-data.frame(sample<-character(),
                              pseudospectrum_candidatelist<-numeric(),
                           pseudospectrum_spclist<-numeric(),
                              retention_time<-numeric(),
                              dot_product<-numeric(),
                           common<-numeric(),
                              stringsAsFactors = F)
for (sample_nr in 1:75){
  for (ps in c(13, 26, 74, 78, 105)){
    rt_query<-candidates2$rt[candidates2$pseudospec==ps][1]
    rt_query_range<-number_range(as.numeric(rt_query)-0.5, as.numeric(rt_query)+0.5)
    ps_query<-grep(pattern=rt_query_range, rt.spclist)
      
        for(hits in ps_query){
          if(!is.null(spclist.by.sample.spec.norm[[sample_nr]][[hits]])){
      product_temp<-compareSpectra(mannitol[[1]],
                                   spclist.by.sample.spec.norm[[sample_nr]][[hits]],
                                   fun = "dot")
      common_temp<-compareSpectra(mannitol[[1]],
                                   spclist.by.sample.spec.norm[[sample_nr]][[hits]],
                                   fun = "common")
      pseudospectrum_candidatelist<-ps
      pseudospectrum_spclist<-hits
      retention_time_temp<-rt.spclist[hits]
      sample_temp<-colnames(spclist[[1]])[as.numeric(sample_nr)+1]
      mannitol_table<-rbind(mannitol_table, c(sample_temp, pseudospectrum_candidatelist,
                                              pseudospectrum_spclist,
                                                    retention_time_temp, product_temp, common_temp))
      colnames(mannitol_table)<-c("sample", "pseudospectrum_candidatelist", "pseudospectrum_spclist", 
                                "retention_time", "dot_product", "common")
        }
    }
  
  }
}


# results (summary_mannitol) suggest that candidate ps 13 is mannitol
# plot mannitol
spcs_mannitol_plot<-allpks_melted[allpks_melted$pcgroup=="13",]

spcs_mannitol_plotA<-spcs_mannitol_plot[which(spcs_mannitol_plot$variable %in% sample_names[1:12]),]
spcs_mannitol_plotA$variable<-factor(spcs_mannitol_plotA$variable, 
                                        levels=c("A_WC4E_20201120_022...Copy", "A_WC5E_20201118_023...Copy", "A_WC6E_20201118_028...Copy",
                                                 "A_WC1E_20201119_060...Copy", "A_WC2E_20201119_058...Copy" , "A_WC3E_20201120_013...Copy",
                                                
                                                "A_WC1I_20201120_006...Copy", "A_WC2I_20201119_069...Copy", "A_WC3I_20201118_034...Copy",
                                                "A_WC4I_20201120_019...Copy", "A_WC5l_20201118_013...Copy", "A_WC6I_20201119_052...Copy"))
spcs_mannitol_plotA$sample_group<-spcs_mannitol_plotA$variable
levels(spcs_mannitol_plotA$sample_group)<-c(rep("A_WC_Edark", 3),rep("A_WC_Elight",3), rep("A_WC_I", 6))
spcs_mannitol_plotA$positionX<-spcs_mannitol_plotA$variable
levels(spcs_mannitol_plotA$positionX)<-c(0.8, 1, 1.2, 1.8, 2, 2.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4)
spcs_mannitol_plotA$positionXnum<-NA
for (i in 1:length(spcs_mannitol_plotA$positionX)){
  xnum<-as.numeric(as.character(spcs_mannitol_plotA$positionX[i]))
  spcs_mannitol_plotA$positionXnum[i]<-as.numeric(xnum)
}


scatterplot3d(x = spcs_mannitol_plotA$positionXnum, y = as.numeric(spcs_mannitol_plotA$mz), z = spcs_mannitol_plotA$value, type = "h",
              pch = "", grid = F, angle=50)

ggplot(candidates2[candidates2$pseudospec==13,])+
  geom_bar(aes(x=))


fucose<-mona.spec[which(mona.names=="fucose")][1]

fucose_table<-data.frame(sample<-character(),
                           pseudospectrum_candidatelist<-numeric(),
                           pseudospectrum_spclist<-numeric(),
                           retention_time<-numeric(),
                           dot_product<-numeric(),
                         common<-numeric(),
                           stringsAsFactors = F)
for (sample_nr in 1:75){
  for (ps in c(13, 26, 74, 78, 105)){
    rt_query<-candidates2$rt[candidates2$pseudospec==ps][1]
    rt_query_range<-number_range(as.numeric(rt_query)-0.5, as.numeric(rt_query)+0.5)
    ps_query<-grep(pattern=rt_query_range, rt.spclist)
    
    for(hits in ps_query){
      if(!is.null(spclist.by.sample.spec.norm[[sample_nr]][[hits]])){
        product_temp<-compareSpectra(fucose[[1]],
                                     spclist.by.sample.spec.norm[[sample_nr]][[hits]],
                                     fun = "dot")
        common_temp<-compareSpectra(fucose[[1]],
                                     spclist.by.sample.spec.norm[[sample_nr]][[hits]],
                                     fun = "common")
        pseudospectrum_candidatelist<-ps
        pseudospectrum_spclist<-hits
        retention_time_temp<-rt.spclist[hits]
        sample_temp<-colnames(spclist[[1]])[as.numeric(sample_nr)+1]
        fucose_table<-rbind(fucose_table, c(sample_temp, pseudospectrum_candidatelist, 
                                            pseudospectrum_temp,
                                            retention_time_temp, product_temp, common_temp))
        colnames(fucose_table)<-c("sample", "pseudospectrum_candidatelist", "pseudospectrum_spclist", 
                                    "retention_time", "dot_product", "common")
      }
    }
    
  }
}
# results (summary_fucose) suggest that candidate ps 78/spclist ps 90, rt 1069 could be fucose


# looking for a match for candidate pseudospectrum 26 (spclist ps 24)

db_candidate_spectra_ps26<-list()
db_candidate_compounds_ps26<-vector(mode="character")
hits=0
mz=c(205, 217, 334)
for (entry in 1:length(mona.spec)){
  if (mz[1] %in% mona.spec[[entry]]@mz){
    if (mz[2] %in% mona.spec[[entry]]@mz){
      if(mz[3] %in% mona.spec[[entry]]@mz){
        hits=hits+1
         db_candidate_spectra_ps26[[hits]]<-mona.spec[entry]
         db_candidate_compounds_ps26[hits]<-mona.names[entry]
         print(hits)
         
      }
   
    }
  }
}
candidateps26_hits<-data.frame(sample=character(),
                     compound=character(),
                     dot_product=numeric(),
                     common=numeric(),
                     stringsAsFactors = F)
for (sample_nr in 1:75){
  for (i in 1:length(db_candidate_spectra_ps26)){
    sample_temp<-colnames(spclist[[1]])[as.numeric(sample_nr)+1]
    compound_temp<-db_candidate_compounds_ps26[i]
  product_temp<-compareSpectra(db_candidate_spectra_ps26[[i]][[1]],
                               spclist.by.sample.spec.norm[[sample_nr]][[24]],
                               fun = "dot")
  common_temp<-compareSpectra(db_candidate_spectra_ps26[[i]][[1]],
                              spclist.by.sample.spec.norm[[sample_nr]][[24]],
                              fun = "common")
  
 candidateps26_hits<-rbind(candidateps26_hits, c(sample_temp, compound_temp,
                                           product_temp, common_temp))
  colnames(candidateps26_hits)<-c("sample", "compound", "dot_product", "common")
  }
}
candidateps26_hits$dot_product<-as.numeric(candidateps26_hits$dot_product)
candidateps26_hits$common<-as.numeric(candidateps26_hits$common)
candidateps26_hits_select<-candidateps26_hits[candidateps26_hits$common>=10,]
candidateps26_hits_select<-candidateps26_hits_select[candidateps26_hits_select$dot_product>=0.5,]
unique(candidateps26_hits_select$compound)
# candidate ps 26 (spclist ps 24), rt 1104 sec / 18.4 min
# the results suggest that candidate ps 26 could be maltose (but the retention time doesn't really fit)
# or L-Iditol, N-Acetyl-glucosamine, Trehalose, D-Xylulose, D-Galacturonic acid
# or sinapic acid


for (ps in c(13, 26, 74, 78, 105)){
  rt_query<-candidates2$rt[candidates2$pseudospec==ps][1]
  rt_query_range<-number_range(as.numeric(rt_query)-0.5, as.numeric(rt_query)+0.5)
  ps_query<-grep(pattern=rt_query_range, rt.spclist)
  
  for(hits in ps_query){
    if(!is.null(spclist.by.sample.spec.norm[[1]][[hits]])){
      product_temp<-compareSpectra(fucose[[1]],
                                   spclist.by.sample.spec.norm[[1]][[hits]],
                                   fun = "dot")
      common_temp<-compareSpectra(fucose[[1]],
                                  spclist.by.sample.spec.norm[[1]][[hits]],
                                  fun = "common")
      pseudospectrum_candidatelist<-ps
      pseudospectrum_spclist<-hits
      retention_time_temp<-rt.spclist[hits]
    }
  }
}
# look at ion intensities in ps 105
par(mfrow=c(3,2))
plots<-list()
for(i in 1:6){
  plots[[i]]<-plot(spclist.by.sample.spec.norm[[i]][[24]])
  return(list(plots))
}
plots


db_candidate_spectra_ps105<-list()
db_candidate_compounds_ps105<-vector(mode="character")
hits=0
mz=c(129, 217, 243)
for (entry in 1:length(mona.spec)){
  if (mz[1] %in% mona.spec[[entry]]@mz){
    if (mz[2] %in% mona.spec[[entry]]@mz){
      if(mz[3] %in% mona.spec[[entry]]@mz){
        hits=hits+1
        db_candidate_spectra_ps105[[hits]]<-mona.spec[entry]
        db_candidate_compounds_ps105[hits]<-mona.names[entry]
        print(hits)
        
      }
      
    }
  }
}
candidateps105_hits<-data.frame(sample=character(),
                               compound=character(),
                               dot_product=numeric(),
                               common=numeric(),
                               stringsAsFactors = F)
# selection of with 70-95% highest intensity in mean of 
for (sample_nr in 1:75){
  for (i in 1:length(db_candidate_spectra_ps105)){
    sample_temp<-colnames(spclist[[1]])[as.numeric(sample_nr)+1]
    compound_temp<-db_candidate_compounds_ps105[i]
    product_temp<-compareSpectra(db_candidate_spectra_ps105[[i]][[1]],
                                 spclist.by.sample.spec.norm[[sample_nr]][[24]],
                                 fun = "dot")
    common_temp<-compareSpectra(db_candidate_spectra_ps105[[i]][[1]],
                                spclist.by.sample.spec.norm[[sample_nr]][[24]],
                                fun = "common")
    
    candidateps105_hits<-rbind(candidateps105_hits, c(sample_temp, compound_temp,
                                                    product_temp, common_temp))
    colnames(candidateps105_hits)<-c("sample", "compound", "dot_product", "common")
  }
}
candidateps105_hits$dot_product<-as.numeric(candidateps105_hits$dot_product)
candidateps105_hits$common<-as.numeric(candidateps105_hits$common)
candidateps105_hits_select<-candidateps105_hits[candidateps105_hits$common>=10,]
candidateps105_hits_select<-candidateps105_hits_select[candidateps105_hits_select$dot_product>=0.5,]
unique(candidateps105_hits_select$compound)



# 5. looking for phlorotannins -----------------------------------------------


standard_ps_table<-data.frame(Compound = character(),
                              Pseudospectrum = numeric(),
                              stringsAsFactors = F)

ref_spec_phloroglucinol<-mona.spec[which(mona.names=="phloroglucinol")]

plot(ref_spec_phloroglucinol[[1]])
phloroglucinol_rt<-1008
phloroglucinol_mz<-342

find_ps_nr(ref_spec = ref_spec_phloroglucinol, rtmin = 800, rtmax = 1400, mz = phloroglucinol_mz, 
           sample_nr = 1, compound_name = "phloroglucinol")

standard_ps_table<-rbind(standard_ps_table, c("phloroglucinol", 22))
colnames(standard_ps_table)<-c("Compound", "Pseudospectrum")
