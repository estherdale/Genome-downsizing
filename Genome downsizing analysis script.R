##### Genome downsizing analysis #####
# Written by Esther Dale, last run 16 April 2024
#section 6 adapted from a script by Matt Larcombe

### Sections ###
# 1. checking botanical names
# 2. preparing occurrence data
# 3. preparing downsizing data
# 4. plotting figure 1
# 5. phylogenetic path analysis
# 6. density plots
# 7. phylo PCA analysis

###########################################
####### 1. checking botanical names #######
###########################################
# standardises the botanical names in genome downsizing dataset and ALLOTB phylogeny
# prune ALLOTB phylogeny to focal species
library(ape)
library(dplyr)
library(phytools)
library(taxize) #for gnr resolve

#read in genome downsizing data
gd.df <- read.csv("Genome size April 24.csv")
gd.df$Species <- gsub(" $", "", gd.df$Species) #remove spaces present at ends of some of the names
#exclude Pratia (because it has multiploid species)
gd.df <- gd.df[!gd.df$Clade=="Pratia/Lobelia",]
#remove an exotic species
gd.df <- gd.df[!gd.df$Scientific.name.original=="Leptinella scariosa",] #remove an exotic species

gd.names.checked <- gnr_resolve(sci=paste(gd.df$Genus, gd.df$Species, sep=" "), canonical = T, best_match_only = T)

# read in ALLOTB phylogeny from Smith & Brown 2018 DOI https://doi.org/10.1002/ajb2.1019
allotb.tr <- read.tree("ALLOTB.tre")
allotb.df <- as.data.frame(allotb.tr$tip.label, stringsAsFactors = F)
colnames(allotb.df) <- "tip_lab"
allotb.df$genus <- sapply(strsplit(allotb.df$tip_lab, split = "_"), `[`, 1) #make genus column

#subset to genera in gd data
allotb.df <- allotb.df[allotb.df$genus %in% unique(c(gd.df$Genus, gd.df$Clade, "Pratia","Lobelia","Anthosachne", "Elymus")),]
allotb.df$tip_lab <- gsub("_", " ", allotb.df$tip_lab) #remove underscore

# run through gnrs
allotb.checked <- gnr_resolve(sci=allotb.df$tip_lab, canonical = T, best_match_only = T) #this takes hours
saveRDS(allotb.checked, file="GNR output ALLOTB names genomic downsizing genera.Rdata")
#allotb.checked <- readRDS(file="GNR output ALLOTB names genomic downsizing genera.Rdata")

#prune phylogeny to genera in downsizing data
gd.genera.allotb.tr <-keep.tip(allotb.tr, tip=gsub(" ","_", allotb.checked$user_supplied_name))  
pdf("GD genera allotb new.pdf", height=80, width=80)
plot(gd.genera.allotb.tr, type="fan")
tiplabels()
dev.off()

# change names in phylo to checked names
allotb.checked <- allotb.checked[!duplicated(allotb.checked$matched_name2),] #remove duplicates i.e. matched to genus
allotb.tr <- keep.tip(allotb.tr, which(allotb.tr$tip.label %in% gsub(" ","_", allotb.checked$user_supplied_name))) #drop duplicates from phylo
allotb.tr$tip.label <- gsub("_", " ", allotb.tr$tip.label)
allotb.tr$tip.label <- allotb.checked$matched_name2[match(unlist(allotb.tr$tip.label), allotb.checked$user_supplied_name)]

#rename some tips to match or act as proxies for species in the genome downsizing data
allotb.tr$tip.label <- gsub("Azorella elegans", "Azorella hookeri", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Azorella coekaynei", "Azorella cockaynei", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Azorella dichopelala", "Azorella nitens", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Azorella microdonta", "Azorella colensoi", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Anthosachne multiflora", "Anthosachne kingiana multiflora", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Brachyscome walshii", "Brachyscome radicata", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Coprosma perpusilla", "Coprosma perpusilla perpusilla", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Deyeuxia pooides", "Deyeuxia quadriseta", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Deyeuxia crinita", "Deyeuxia aff quadriseta", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Deyeuxia angustifolia", "Deyeuxia avenoides", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Hydrocotyle zongoana", "Hydrocotyle heteromeria", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Hydrocotyle novae-zealandiae", "Hydrocotyle novae-zeelandiae", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Hydrocotyle novae-zeelandiae$", "Hydrocotyle novae-zeelandiae novae-zeelandiae", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Hierochloe pluriflora", "Hierochloe redolens", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Lachnagrostis pilosa", "Lachnagrostis pilosa pilosa", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Leptinella dioica subsp. monoica", "Leptinella dioica dioica 2", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Leptinella pectinata subsp. pectinata", "Leptinella pectinata pectinata 2", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Plantago turficola", "Plantago picta", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Poa himalayana", "Poa aff cita a", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Poa borneensis", "Poa aff cita b", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Pratia nummularia", "Pratia old man", allotb.tr$tip.label)
allotb.tr$tip.label <- gsub("Veronica pauciramosa", "Veronica odora", allotb.tr$tip.label)

#change names in genome downsizing data to checked names
gd.df$binomial<- gd.names.checked$matched_name2[match(unlist(paste(gd.df$Genus, gd.df$Species, sep=" ")), gd.names.checked$user_supplied_name)]

#edit names that aren't formal species
gd.df$binomial[gd.df$Species=="aff. quadriseta" &gd.df$Genus=="Deyeuxia"] <- "Deyeuxia aff quadriseta"
gd.df$binomial[gd.df$Species=="aff cita (a)" &gd.df$Genus=="Poa"] <- "Poa aff cita a"
gd.df$binomial[gd.df$Species=="aff cita (b)" &gd.df$Genus=="Poa"] <- "Poa aff cita b"
gd.df$binomial[gd.df$Species=="novae-zealandiae var novae-zealandiae" &gd.df$Genus=="Hydrocotyle"] <- "Hydrocotyle novae-zeelandiae novae-zeelandiae"

#identify species that need multiple tips and add number name
polyploids <- gd.df[gd.df$Base.species==0,]
multiploids <- polyploids[duplicated(polyploids$binomial),]
multiploids$binomial <- gsub("$", " 2", multiploids$binomial) #add 2 to end of binomial to differentiate from matching tip
multiploids$binomial[grep("Leptinella squalida mediana 2", multiploids$binomial)[2]] <- gsub("2", "3", multiploids$binomial[grep("Leptinella squalida mediana 2", multiploids$binomial)[2]]) 
multiploids$binomial[grep("Leptinella squalida mediana 2", multiploids$binomial)[2]] <- gsub("2", "4", multiploids$binomial[grep("Leptinella squalida mediana 2", multiploids$binomial)[2]]) 
polyploids <- rbind(polyploids[!duplicated(polyploids$binomial),], multiploids)

#identify which tips are in the phylogeny
polyploids$allotb <- rep(0, nrow(polyploids))
polyploids$allotb[polyploids$binomial %in% allotb.tr$tip.label] <-1 

# prune phylogeny to species that match with genome downsizing species
allotb.tr.pruned <- keep.tip(allotb.tr, tip = polyploids$binomial[polyploids$allotb==1])
allotb.tr.pruned$node.label <- NULL

#add new tips for multiploid species
multiploid.list <- polyploids$binomial[grep("[[:digit:]]$", polyploids$binomial)]
allotb.tr.enhanced <- allotb.tr.pruned
#cycle through adding tips (except #3 and 4 Leptinella squalida mediana will do manually)
for(i in c(1:4,7:9)){
  allotb.tr.enhanced <- bind.tip(allotb.tr.enhanced, tip.label = multiploid.list[i], where=grep(gsub(" [[:digit:]]$", "",multiploid.list[i]), allotb.tr.enhanced$tip.label), edge.length = min(allotb.tr.pruned$edge.length))
}
#manually add Leptinella squalida mediana 3 and 4
allotb.tr.enhanced <- bind.tip(allotb.tr.enhanced, tip.label = multiploid.list[5], where=getMRCA(allotb.tr.enhanced, grep(gsub(" [[:digit:]]$", "",multiploid.list[7]), allotb.tr.enhanced$tip.label)),edge.length = min(allotb.tr.pruned$edge.length) )
allotb.tr.enhanced <- bind.tip(allotb.tr.enhanced, tip.label = multiploid.list[6], where=getMRCA(allotb.tr.enhanced, grep(gsub(" [[:digit:]]$", "",multiploid.list[8]), allotb.tr.enhanced$tip.label)),edge.length = min(allotb.tr.pruned$edge.length) )
#manually add a tip for Leptinella squalida
LEPsqua.node <- getMRCA(allotb.tr.enhanced, grep("Leptinella squalida",allotb.tr.enhanced$tip.label))
#tip length is difference between one of the tips and the node where tip will be added
LEPsqua.length <- nodeheight(allotb.tr.enhanced, grep("Leptinella squalida squalida", allotb.tr.enhanced$tip.label))-nodeheight(allotb.tr.enhanced, LEPsqua.node)
allotb.tr.enhanced <- bind.tip(allotb.tr.enhanced, tip.label = "Leptinella squalida", where=LEPsqua.node,edge.length = LEPsqua.length )

#combine fixed names for polyploids with unchanged base species
gd.df <- rbind(gd.df[gd.df$Base.species==1,], polyploids[,1:11])

#fix clade spelling mistakes
gd.df$Clade <- gsub("Elymus/Anthoschne", "Anthosachne/Connorochloa", gd.df$Clade) #rename anthosachne clade
gd.df$Clade <- gsub("Lachnogrostis", "Lachnagrostis", gd.df$Clade) #fix spelling

#save
write.tree(allotb.tr.pruned, file="ALLOTB_gd_species_checked.tre")
write.tree(allotb.tr.enhanced, file="ALLOTB_gd_species_checked_with_multiploids.tre")
saveRDS(gd.df, file="Genome downsizing data names checked.Rdata")

## make table for supplementary materials ##
table.df <- gd.df[,c(1,9,11)]
table.df<-  table.df %>% group_by(Clade, Ploidy) %>% summarise(species=n())
habit.table <- unique(gd.df[,c(1,5)])
table.df <- left_join(table.df, unique(gd.df[,c(1,5)]), by="Clade") 
families <- tax_name(get="family", table.df$Clade)
table.df$Family <- families$family
table.df$Family[table.df$Clade %in% c("Anthosachne/Connorochloa", "Deyeuxia", "Hierochloe")] <- "Poaceae"
table.df$Family[table.df$Clade == "Azorella"] <- "Apiaceae"
table.df$Family[table.df$Clade == "Libertia"] <- "Iridaceae"
table.df$Family[table.df$Clade %in% c("Brachyscome", "Leptinella")] <- "Asteraceae"
table.df <- table.df[order(table.df$Family),] #sort by family
#add chromosome numbers
table.df$chromosome.number <- rep(0,nrow(table.df))
for(i in 1:nrow(table.df)){
  chromosomes.i <- unique(gd.df$X2n[gd.df$Clade==table.df$Clade[i]&gd.df$Ploidy==table.df$Ploidy[i]])
  if(length(chromosomes.i)>1){chromosomes.i <- paste(min(chromosomes.i), max(chromosomes.i), sep="-")}
  table.df$chromosome.number[i] <- chromosomes.i
}
table.df <- table.df[,c(5,1,4,2,6,3)] 
write.csv(table.df, file="Table of number of species by ploidy level and clade.csv")

###################################################
########## 2. preparing occurrence data ###########
###################################################
### prepping AVH data and GBIF for TTR modelling for polyploid series
# AVH data for records in NZ and GBIF data from anything in NZ, subset to the species of interest
# checks species names, filter out points in the sea/missing/rounded etc.

library(plyr)
library(dplyr)
library(CoordinateCleaner)
library(taxize)
library(rgbif)
library(stringr)

#read in focal species list
species.list <- readRDS("data/Genomic downsizing data names checked.Rdata")
species.list <- species.list[!grepl("[[:digit:]]", species.list$binomial),] #exclude multiploids (have numbers in binomials)
species.list <- species.list[!grepl(" aff ", species.list$binomial),] #exclude affiliated taxa
species.list$words <- sapply(strsplit(species.list$binomial, " "), length)
species.list <- species.list[species.list$words %in% 2:3,] #include genus-only taxa
species.list$subspecies <- rep("FALSE", nrow(species.list))
species.list$subspecies[species.list$words>2] <- "TRUE" #identifier for subspecies and vars
gd.species <- unique(species.list$binomial)

#read in AVH data and filter to focal genera
#AVH data from DOI https://doi.org/10.26197/5f101cb43c935 and DOI https://doi.ala.org.au/doi/c2c95558-3f35-444e-af3e-fd75b06d2deb
avh.data <- read.csv("records-2021-09-07.csv", stringsAsFactors = F) #read in csv
focal.genera <- unique(c(species.list$Genus, "Elymus", "Anthosachne", "Hymenanthera", "Lobelia", "Hebe", "Schizeilema"))
avh.data <- avh.data[avh.data$genus %in% focal.genera,] #subset to genera in data
avh.data <- avh.data[,which(colnames(avh.data) %in% c("scientificName", "decimalLongitude", "decimalLatitude"))] #remove NAs
avh.data <- avh.data[complete.cases(avh.data),]
avh.data$scientificName <- gsub("^Coprosma baueri$", "Coprosma lucida", avh.data$scientificName) #change Coprosma baueri to lucida
avh.data$scientificName <- gsub("^Libertia pulchella$", "Libertia micrantha", avh.data$scientificName)
avh.data$scientificName <- gsub("^Lobelia perpusilla$", "Pratia perpusilla", avh.data$scientificName)
avh.data$scientificName <- gsub("^Lobelia angulata$", "Pratia angulata", avh.data$scientificName)
avh.data$scientificName <- gsub("^Lobelia arenaria$", "Pratia arenaria", avh.data$scientificName)

#check AVH botanical names and subset to focal species
avh.checked <- gnr_resolve(sci=unique(avh.data$scientificName), canonical = T, best_match_only = T)
avh.data$binomial <-avh.checked$matched_name2[match(avh.data$scientificName, avh.checked$user_supplied_name)]
avh.data <- avh.data[avh.data$binomial %in% gd.species,]
saveRDS(avh.checked, file="AVH names checked genomic downsizing genera.Rdata")

# download GBIF data for all plant samples from NZ and filter to focal genera 
gbif.data <- occ_download_get(key = "0011617-210819072339941", overwrite = TRUE) %>%  #get key from dataset downloaded on gbif
  occ_download_import(gbif.data, quote="", na.strings = c("", NA))
gbif.data <- gbif.data[gbif.data$genus %in% focal.genera,] #subset to focal genera
gbif.data <- gbif.data[,which(colnames(gbif.data) %in% c("scientificName", "decimalLongitude", "decimalLatitude"))]
gbif.data <- gbif.data[complete.cases(gbif.data),]

#check GBIF botanical names and subset to focal species
species.list.gbif.raw<- unique(gbif.data$scientificName)
##check for issues with gnr request, last name printed is throwing an error
for(i in 249:length(species.list.gbif.raw)){
  print(species.list.gbif.raw[i])
  test <- gnr_resolve(sci=species.list.gbif.raw[i],canonical = T, best_match_only = T)
}
species.list.gbif.raw <- species.list.gbif.raw[-248] #remove name that is causing an error
gbif.checked <- gnr_resolve(sci=species.list.gbif.raw, canonical = T, best_match_only = T)
gbif.data$binomial <-gbif.checked$matched_name2[match(gbif.data$scientificName, gbif.checked$user_supplied_name)]
gbif.data <- gbif.data[gbif.data$binomial %in% gd.species,]

#combine datasets#
avh.data <- avh.data[,c(3,1,2,4)] #rearrange columns to same order as gbif
combined.data <- rbind(avh.data, gbif.data)

### data filtering ###
# check for coords in capital city, near institutions, in the sea etc.
combined.data <- clean_coordinates(combined.data, lon="decimalLongitude", lat="decimalLatitude", species = "binomial", tests=c("capitals", "centroids", "equal", "institutions", "zeros"), value = "clean")
combined.data <- cc_sea(combined.data, lon="decimalLongitude", lat="decimalLatitude", value="clean") #remove occurrence records in the sea, did it separately cos sometimes the sea doesn't work

#### subset columns and rearrange columns ###
combined.data <- combined.data[,4:2] #need name then long then lat
colnames(combined.data) <- c("sp", "lon", "lat")
combined.data <- distinct(combined.data) #exclude duplicates

# check for subspecies
combined.data$words <- sapply(strsplit(combined.data$sp, " "), length)
combined.data <- combined.data[combined.data$words%in%c(2,3),] #subset to records with two (species) or three (subspecies/varieties) words

# check how many points per taxa and exclude with too few
pr.counts <- as.data.frame(table(combined.data$sp))
species.list<- full_join(species.list, pr.counts, by=c("binomial"= "Var1"))
colnames(species.list)[14] <- "pr.points"
combined.data <- combined.data[which(combined.data$sp %in% species.list$binomial[species.list$pr.points>10]),]

combined.data <- combined.data[,1:3]#remove words column

### Split into different files by species ###
combined.data$sp <- gsub(" ", "_", combined.data$sp) #change spaces to underscores
a<-levels(as.factor(combined.data$sp)) #list of species to split into

for(i in 1:length(a)) 
{
  d <- combined.data[combined.data$sp==a[i],2:3] #subset to ith species and lon and lat columns
  saveRDS(d, file = paste("PR.", as.character(a[i]),".Rdata",sep=""))
  print(a[i])
  print(nrow(d))
}

# now the TTR.sdm model needs to be run on these PR files
# in the repo we've included the TTR.sdm output data frame "sdm.stats.table.jun22.Rdata"

###################################################
########## 3. preparing downsizing data ###########
###################################################
# calculates genome downsizing metrics for all non-base ploidy species (i.e. polyploids)
# also combines genome downsizing data and TTR data, and calculates difference between base and focal species for TTR traits
#includes clade ages and speciation rates for clades

library(tidyr)
library(dplyr)
library(geiger)

#load genome downsizing data
gd.df <- readRDS("Genome downsizing data names checked.Rdata")

# add in clade age
clade.df <- read.csv(file="Genus size in NZ genomic downsizing.csv")
ages.df <- read.csv("NZ flora colonisation ages Heenan.csv")
#this age csv comes from Heenan and McGlone 2019 DOI https://doi.org/10.1080/0028825X.2019.1632356

#add ages to clade size data and make new row that combined Anthosachne and Connorchloa as a single clade (they have the same age values already)
ages.df$genus <- sapply(strsplit(ages.df$Taxon, " "), `[`, 1)
ages.df <- ages.df[,c(13,4,7)]
ages.df <- ages.df[-c(57,232,261,358,359),] #remove duplicates we don't want of Lachnagrostis, Plantago, Lobelia, Azorella
clade.df <- clade.df[c(1:18,1),] #add another anthosachne row to act for the clade
clade.df <- left_join(clade.df, ages.df, by=c("Genus"="genus"))
clade.df$Genus[clade.df$Genus=="Anthosachne"][2] <- "Anthosachne/Connorochloa" #rename second Anthosachne entry to be for the clade
clade.df$Clade.size[clade.df$Genus=="Anthosachne/Connorochloa"] <- sum(clade.df$Clade.size[clade.df$Genus %in% c("Anthosachne", "Connorochloa")]) #sum genus counts together for anthosachne/connorochloa clade count

#calculate diversification rate and add to gd.df
clade.df$diversification <- bd.ms(time=clade.df$Crown.age, n=clade.df$Clade.size, missing=0, crown=T)
gd.df <- left_join(gd.df, clade.df, by=c("Clade"="Genus"))

## read in TTR data
load("sdm.stats.table.jun22.Rdata")
ttr.df <- sdm.stats.table
ttr.df$sp <- gsub("_", " ", ttr.df$sp)
#exclude species with very poor fits
good.fit <- ttr.df$d.fneg<=ttr.df$d.tpos
confusion.matrix <- cbind(ttr.df[,1:5], good.fit)
confusion.matrix <- left_join(gd.df[1:200,11:10], confusion.matrix, by=c("binomial"="sp"))
confusion.matrix$ttr.not.fitted <- !complete.cases(confusion.matrix$d.tpos)  
confusion.matrix$good.fit <- gsub("FALSE", "0", confusion.matrix$good.fit)
confusion.matrix$good.fit <- gsub("TRUE", "1", confusion.matrix$good.fit)
confusion.matrix$ttr.not.fitted <- gsub("FALSE", "0", confusion.matrix$ttr.not.fitted)
confusion.matrix$ttr.not.fitted <- gsub("TRUE", "1", confusion.matrix$ttr.not.fitted)
write.csv(confusion.matrix, file="TTR confusion matrix.csv")

ttr.df <- ttr.df[good.fit,]
gd.df <- left_join(gd.df, ttr.df[,c(1:8,33:56)], by=c("binomial"="sp"))

#separate base species and calculate mean DNA
base.df <- gd.df[gd.df$Base.species==1,] %>% group_by(Clade) %>%   summarise(
  base.total.DNA = mean (as.numeric(as.character(Total.DNA)), na.rm=T),
  base.2n = median(as.numeric(as.character(X2n)), na.rm=T),
  base.chromosome.size = (mean (as.numeric(as.character(Total.DNA)), na.rm=T)/median(as.numeric(as.character(X2n)), na.rm=T)),
  base.ploidy = median(as.numeric(as.character(Ploidy)), na.rm=T),
  base.species.number = n(),
  base.NZ.count = mean(NZ.count, na.rm=T),
  base.resample.count = mean (resample.count, na.rm=T),
  base.tmax1.back = mean(tmax1.back, na.rm=T),
  base.tmax2.back = mean(tmax2.back, na.rm=T),
  base.tmax3.back = mean(tmax3.back, na.rm=T),
  base.tmax4.back = mean(tmax4.back, na.rm=T),
  base.q1.back = mean(q1.back, na.rm=T),
  base.q2.back = mean(q2.back, na.rm=T),
  base.w11.back = mean(w11.back, na.rm=T),
  base.w12.back = mean(w12.back, na.rm=T),
  base.ns1.back = mean(ns1.back, na.rm=T),
  base.ns2.back = mean(ns2.back, na.rm=T),
  base.tmean1.back = mean(tmean1.back, na.rm=T),
  base.tmean2.back = mean(tmean2.back, na.rm=T),
  base.w21.back = mean(w21.back, na.rm=T),
  base.w22.back = mean(w22.back, na.rm=T),
  base.w23.back = mean(w23.back, na.rm=T),
  base.w24.back = mean(w24.back, na.rm=T),
  base.nsoil1.back = mean(nsoil1.back, na.rm=T),
  base.nsoil2.back = mean(nsoil2.back, na.rm=T),
  base.tmin1.back = mean(tmin1.back, na.rm=T),
  base.tmin2.back = mean(tmin2.back, na.rm=T),
  base.tmin3.back = mean(tmin3.back, na.rm=T),
  base.tmin4.back = mean(tmin4.back, na.rm=T),
  base.tmean21.back = mean(tmean21.back, na.rm=T),
  base.tmean22.back = mean(tmean22.back, na.rm=T)
  
)

# add base values to genome downsizing df
downsizing.df <- gd.df[gd.df$Base.species==0,c(1:9,11:46)]
downsizing.df <- left_join(downsizing.df, base.df, by="Clade")

#calculate downsizing metrics
downsizing.df$base.per.copy <- downsizing.df$base.total.DNA/downsizing.df$base.ploidy
downsizing.df$per.copy <- downsizing.df$Total.DNA/downsizing.df$Ploidy
downsizing.df$no.downsizing <- downsizing.df$Ploidy*downsizing.df$base.per.copy
downsizing.df$downsizing <- (downsizing.df$base.per.copy*downsizing.df$Ploidy)-downsizing.df$Total.DNA
downsizing.df$downsizing.percent <- 100*downsizing.df$downsizing/downsizing.df$no.downsizing

#save
saveRDS(downsizing.df, file="Downsizing metrics with TTR results.Rdata")
saveRDS(gd.df, file="Genome downsizing data with TTR results.Rdata")
saveRDS(ttr.df, file="TTR data poor fits removed.Rdata")

##########################################
########## 4. plotting figure 1 ##########
##########################################
# Plot Fig 1 percent genome downsizing ~ ploidy and habit
# fits pgls testing genome downsizing~ploidy*habit effect
library(caper)

# read in data
downsizing.df <- readRDS("Downsizing metrics with TTR results.Rdata")

# set colours
life.form.cols <- c("#D55E00", "#F0E442", "#56B4E9")
border.cols <- c("#660000", "#CC9900", "#003366")

#change ploidy to a factor and set levels so boxplot is to scale
ploidy.factor <- c(downsizing.df$Ploidy, seq(from=4, to=38, by=2)) #add in all the numbers we want in the axis
ploidy.factor <- as.factor(as.numeric(ploidy.factor)) #set as factor
ploidy.factor <- ploidy.factor[1:nrow(downsizing.df)] #subset back to actual ploidy data

box.width <- seq(from=0.01, to=0.5, length=18)
box.width2 <- c(rep(0.02,18))
box.width3 <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)

# pgls to test downsizing, ploidy and habit relationship
gd.cd <- comparative.data(allotb.tr.pruned, downsizing.df[,c(5,9,10,81)], names.col = "binomial")
pgls1.model <- pgls(downsizing.percent~Ploidy*Life.form, data=gd.cd, lambda = "ML")
pgls2.model <- pgls(downsizing.percent~Life.form, data=gd.cd, lambda = "ML")
pgls3.model <- pgls(downsizing.percent~Ploidy+Life.form, data=gd.cd, lambda = "ML")
pgls4.model <- pgls(downsizing.percent~Ploidy, data=gd.cd, lambda = "ML")

pgls1 <- summary(pgls1.model)
pgls2 <- summary(pgls2.model)

AIC(pgls2.model, pgls3.model, pgls1.model, pgls4.model)

## plotting ###
pdf(height=6, width = 6, file="Genome downsizing by ploidy and lifeform boxplots.pdf")
par(mfrow=c(3,1), mar=c(0,5,4.1,3))
boxplot(downsizing.df$downsizing.percent[downsizing.df$Life.form=="forb"]~ploidy.factor[downsizing.df$Life.form=="forb"],col=life.form.cols[3], border=border.cols[3], xaxt="n", xlab="", ylab="Genome downsizing (%)",ylim=c(90,-20), boxwex=0.3)
segments(x0=1, y0=pgls2$coefficients[1,1],x1=12,y1=pgls2$coefficients[1,1], col=unique(border.cols)[3], lwd=1)
abline(h=0, lty=2, col="grey60")
text(0.2,-15, "(a)")
text(17,84, "forbs", col=border.cols[3], pos=4)

par(mar=c(2.05,5,2.05,3))
boxplot(downsizing.df$downsizing.percent[downsizing.df$Life.form=="grass"]~ploidy.factor[downsizing.df$Life.form=="grass"],col=life.form.cols[2], border=border.cols[2], xaxt="n", xlab="",  ylab="Genome downsizing (%)",ylim=c(90,-20), boxwex=0.3)
segments(x0=1,y0=(pgls2$coefficients[1,1]+pgls2$coefficients[2,1]),x1=18, y1=(pgls2$coefficients[1,1]+pgls2$coefficients[2,1]),col=unique(border.cols)[2], lwd=1)
abline(h=0, lty=2, col="grey60")
text(0.2,-15, "(b)")
text(17,84, "grasses", col=border.cols[2], pos=4)

par(mar=c(4.1,5,0,3))
boxplot(downsizing.df$downsizing.percent[downsizing.df$Life.form=="woody"]~ploidy.factor[downsizing.df$Life.form=="woody"],col=life.form.cols[1], border=border.cols[1], xlab="Ploidy", ylab="Genome downsizing (%)",ylim=c(90,-20), boxwex=0.3)
segments(x0=1,y0=(pgls2$coefficients[1,1]+pgls2$coefficients[3,1]),x1=15, y1=(pgls2$coefficients[1,1]+pgls2$coefficients[3,1]),col=unique(border.cols)[1], lwd=1)
abline(h=0, lty=2, col="grey60")
text(0.2,-15, "(c)")
text(17,84, "woody", col=border.cols[1], pos=4)

dev.off()

downsizing.means <- downsizing.df %>% group_by(Life.form) %>% summarise(mean.downsizing=mean(downsizing), mean.downsizing.percent=mean(downsizing.percent), std.dev.downsizing=sd(downsizing), std.dev.downsizing.percent=sd(downsizing.percent), n=length(downsizing), median.downsizing=median(downsizing), median.percent.downsizing=median(downsizing.percent)) 
downsizing.means$std.err.downsizing <- downsizing.means$std.dev.downsizing/sqrt(downsizing.means$n)
downsizing.means$std.err.downsizing.percent <- downsizing.means$std.dev.downsizing.percent/sqrt(downsizing.means$n)

########################################
#### 5. phylogenetic path analysis #####
########################################
# identifies traits that correlate with genome downsizing and then compares various
# models with a phylogenetic path analysis
## correlation plot for genome downsizing vars

library(ape)
library(dplyr)
library(ellipse)
library(ggcorrplot)
library(phylopath)
library(phytools)
library(Rcpp)

#load data and phylogeny
downsizing.df <- readRDS(file="Downsizing metrics with TTR results.Rdata")
size.df <- read.csv("Genomic downsizing species heights.csv")
tr<- allotb.tr.pruned

#add median leaf size to downsizing data
size.df$leaf.length.median <- rowMeans(size.df[,6:7], na.rm = T) #calculate median
#use similar species for values of affinities
size.df$leaf.length.median[size.df$binomial=="Deyeuxia aff quadriseta"] <- size.df$leaf.length.median[size.df$binomial=="Deyeuxia quadriseta"]
size.df$leaf.length.median[grepl("Poa aff cita", size.df$binomial)] <- size.df$leaf.length.median[size.df$binomial=="Poa cita"]
downsizing.df <- left_join(downsizing.df, size.df[,c(3,5,10)], by="binomial")
#make binary version of life form (phylopath can only handle binary categorical vars)
downsizing.df$Life.form.full <- downsizing.df$Life.form
downsizing.df$Life.form <- gsub("forb","herbaceous",downsizing.df$Life.form)
downsizing.df$Life.form <- gsub("grass","herbaceous",downsizing.df$Life.form)

#subset downsizing.df to exclude multiploid taxa
downsizing.df <- downsizing.df[!grepl("2", downsizing.df$binomial),]
downsizing.df <- downsizing.df[!grepl("3", downsizing.df$binomial),]
downsizing.df <- downsizing.df[!grepl("4", downsizing.df$binomial),]
downsizing.df <- downsizing.df[-which(downsizing.df$binomial=="Leptinella squalida"),]
tr <- keep.tip(tr, tip=downsizing.df$binomial)
rownames(downsizing.df) <- downsizing.df$binomial

#correlation matrix
traits.cor <- cor(downsizing.df[,c(81, 7:9,46:49,12,14, 82:83)], use="pairwise.complete.obs")
colnames(traits.cor) <- c("genome downsizing","genome size", "chromosome number", "ploidy", "base genome size","base chromosome number","base chromosome size",
                            "base ploidy","clade time since colonisation","clade diversification rate", "max height","leaf length")
rownames(traits.cor) <- c("genome downsizing","genome size", "chromosome number", "ploidy", "base genome size","base chromosome number","base chromosome size",
                          "base ploidy","clade time since colonisation","clade diversification rate", "max height","leaf length")

p.mat <- cor_pmat(downsizing.df[,c(81, 7:9,46:49,12,14, 82:83)])
p.mat.inverse <- 1-p.mat #inverse of p values so can indicate significant values not insignificant
colnames(p.mat.inverse) <- c("genome downsizing","genome size", "chromosome number", "ploidy", "base genome size","base chromosome number","base chromosome size",
                          "base ploidy","clade time since colonisation","clade diversification rate", "max height","leaf length")
rownames(p.mat.inverse) <- c("genome downsizing","genome size", "chromosome number", "ploidy", "base genome size","base chromosome number","base chromosome size",
                          "base ploidy","clade time since colonisation","clade diversification rate", "max height","leaf length")

#pdf("Trait correlation heatmap1.pdf", height=10, width=10)
cor.plot <-  ggcorrplot(traits.cor, hc.order = F, outline.color = "white", type="lower", p.mat = p.mat.inverse, sig.level = 0.95,
                         pch.cex = 7, pch.col = "black", colors = c("#046C9A", "white","#F2300F"), pch="*")
ggsave(
  "Trait correlation heatmap.pdf",
  plot = cor.plot,
  width = 10,
  height = 10 )


# without multiploids
models<- define_model_set(
  # all
  all.1=c(downsizing.percent~base.ploidy+base.2n+Life.form+Max_height_m,
           Life.form~base.2n,
           Max_height_m~base.ploidy+base.2n+Life.form,
          base.2n~base.ploidy
  ),
  
  
  all.2=c(downsizing.percent~base.ploidy+base.2n+Life.form+Max_height_m,
          base.ploidy~Life.form+Max_height_m,
          base.2n~Life.form+Max_height_m+base.ploidy),
  
  all.3=c(downsizing.percent~base.ploidy+base.2n+Life.form+Max_height_m,
          base.ploidy~Life.form,
          base.2n~Life.form+base.ploidy),
  
  all.4=c(downsizing.percent~base.ploidy+base.2n+Life.form+Max_height_m,
          Life.form~base.ploidy+base.2n,
          Max_height_m~Life.form
  ),
  #  genetic direct
  genetic.direct.1=c(downsizing.percent~base.ploidy+base.2n,
          base.ploidy~Life.form+Max_height_m,
          base.2n~Life.form+Max_height_m+base.ploidy,
          Max_height_m~Life.form
  ),
  
  genetic.direct.2=c(downsizing.percent~base.ploidy+base.2n,
                     base.ploidy~Life.form,
                     base.2n~Life.form+base.ploidy
  ),
  

#ecological direct
ecological.direct.1=c(downsizing.percent~Life.form+Max_height_m,
         Life.form~base.ploidy+base.2n,
         Max_height_m~base.ploidy+base.2n+Life.form,
         base.2n~base.ploidy
),

ecological.direct.2=c(downsizing.percent~Life.form+Max_height_m,
                      Life.form~base.2n,
                      Max_height_m~base.2n+Life.form
),


  # just ecol
  ecological.only.1=c(downsizing.percent~Life.form+Max_height_m,
            Max_height_m~Life.form
  ),
ecological.only.2=c(downsizing.percent~Life.form
  ),

  #just genetic
  genetic.only.1=c(downsizing.percent~base.ploidy+base.2n,
          base.2n~base.ploidy),
  genetic.only.2=c(downsizing.percent~base.ploidy+base.2n),
  #common across all models
  .common=NULL
)
results<- phylo_path(models, data=downsizing.df, tree=tr, model='lambda', lower.bound=0, upper.bound=1, btol=12.4, method="logistic_MPLE")

pdf("Path analysis model weights.pdf")
plot(summary(results))
dev.off()

#set up variable labels and positions
positions <- data.frame(
  name = c('base.2n', 'base.ploidy', 'Life.form', 'Max_height_m', 'downsizing.percent'),
  x = c(3,4,3,4,5),
  y = c(3,3,1,1,2)
)

labels <-c("base\nploidy", "base\n2n", "life\nform", "max\nheight","genome\ndownsizing")
names(labels) <- c("base.ploidy", "base.2n", "Life.form", "Max_height_m", "downsizing.percent")


pdf("Phylogenetic path models genomic downsizing.pdf", height=15, width=17)
plot_model_set(models, manual_layout =  positions, labels = labels, box_x = 16, box_y = 16, text_size = 3, nrow=3)
dev.off()

#supported models all.1 and all.3, but all.3 is simpler so will select

#get estimates for the selected model
estimates1 <- choice(results, choice="all.3", boot=1000)
#estimates2 <- choice(results, choice="all.1", boot=1000)

plot1 <- coef_plot(estimates1, reverse_order = TRUE) + 
  ggplot2::coord_flip() + 
  ggplot2::theme_bw()+ggtitle("model all.3") 

saveRDS(estimates1, file="Estimated paths model all.3.Rdata")

pdf("Genome downsizing intrinsic phylogenetic path analysis.pdf", height=7, width=9)
plot(estimates1,box_x = 16, box_y = 16, manual_layout = positions, text_size = 3, labels=labels, curvature =0.05)
dev.off()


######################################
########## 6. density plots ##########
######################################
#estimates niche size for base species and polyploids with high and low genome downsizing
# uses TTR predicted distributions

library(arm)
library(dplyr)
#----------Make niche area density plots using bayesglm ----------------------------

#load genome downsizing and TTR data
gd.df <- readRDS("Genome downsizing data with TTR results.Rdata")
downsizing.df <- readRDS("Downsizing metrics with TTR results.Rdata")

ttr2 <- gd.df

#subset and combine ploidy and ttr data, rename base species to labels
ttr2$Base.species <- gsub("^1$", "base", ttr2$Base.species)
ttr2$Base.species <- gsub("^0$", "polyploid", ttr2$Base.species)
ttr2 <- ttr2[complete.cases(ttr2),] #exclude species without TTR data

# add in % downsizing for polyploids
ttr2 <- left_join(ttr2, downsizing.df[,c(10, 81)], by="binomial")
ttr2$downsizing.percent[ttr2$Base.species=="base"] <- NA #make sure any species with both base and ployploid don't get downsizing % for base row

### repeat for diff life forms ####
life.forms <- c("woody", "grass", "forb", "all")

bayesglm.res.df <- data.frame("life.form"=character(), "niche.type"=character(), "downsizing.level"=character(), "credible.interval.lower"=numeric(), "credible.interval.upper"=numeric())

for(i in 1:4){
 if(i==4){gd.i <- ttr2}else{gd.i <- ttr2[ttr2$Life.form==life.forms[i],]}
   
  gd.i$downsizing.levels <- as.character(cut(gd.i$downsizing.percent, breaks =c(-20,median(gd.i$downsizing.percent, na.rm = T), 80), labels=c("low","high") ))
  gd.i$downsizing.levels[gd.i$Base.species=="base"] <-"base" 
  
  f.wp.i <- bayesglm (NZ.count ~ factor(downsizing.levels),
                      prior.scale=Inf, prior.df=Inf,
                      data=gd.i)
  f.op.i<- bayesglm (resample.count ~  factor(downsizing.levels),
                     prior.scale=Inf, prior.df=Inf,
                     data=gd.i)
   
  fwp.sim.c.i <- coef(sim(f.wp.i,n.sims=50000))
  fop.sim.c.i <- coef(sim(f.op.i,n.sims=50000))
  
  myadjust<-1 #this scales the y axis so you can see the plot, may need to fiddle with this
  #constant plus polyploid/life form dummy variables
  base.wp.i<-density(fwp.sim.c.i[,1],adjust=myadjust)                 
  polyploid.low.wp.i<-density(fwp.sim.c.i[,1] + fwp.sim.c.i[,3],adjust=myadjust) 
  polyploid.high.wp.i<-density(fwp.sim.c.i[,1]+ fwp.sim.c.i[,2],adjust=myadjust) 
  
  base.op.i<-density(fop.sim.c.i[,1],adjust=myadjust)                 
  polyploid.low.op.i<-density(fop.sim.c.i[,1] + fop.sim.c.i[,3],adjust=myadjust) 
  polyploid.high.op.i<-density(fop.sim.c.i[,1] + fop.sim.c.i[,2],adjust=myadjust) 
  
  nz.95ci.low <- quantile(fwp.sim.c.i[,1] + fwp.sim.c.i[,3], c(0.025,0.975))
  nz.95ci.high <- quantile(fwp.sim.c.i[,1]+ fwp.sim.c.i[,2], c(0.025,0.975))
  
  resampled.95ci.low <- quantile(fop.sim.c.i[,1] + fop.sim.c.i[,3], c(0.025,0.975))
  resampled.95ci.high <- quantile(fop.sim.c.i[,1]+ fop.sim.c.i[,2], c(0.025,0.975))
  
  new.row.resample <- rbind(c(life.forms[i], "resample", "high", resampled.95ci.high), c(life.forms[i], "resample", "low", resampled.95ci.low))
  new.row.nz <- rbind(c(life.forms[i], "nz", "high", nz.95ci.high), c(life.forms[i], "nz", "low", nz.95ci.low))
  
  bayesglm.res.df[c(nrow(bayesglm.res.df)+1,nrow(bayesglm.res.df)+2),] <- new.row.resample
  bayesglm.res.df[c(nrow(bayesglm.res.df)+1,nrow(bayesglm.res.df)+2),] <- new.row.nz
  
  if(life.forms[i]=="woody"){
    base.wp.woody <- base.wp.i
    polyploid.low.wp.woody <- polyploid.low.wp.i
    polyploid.high.wp.woody <- polyploid.high.wp.i
    
    base.op.woody <- base.op.i
    polyploid.low.op.woody <- polyploid.low.op.i
    polyploid.high.op.woody <- polyploid.high.op.i
  }
  if(life.forms[i]=="grass"){
    base.wp.grass <- base.wp.i
    polyploid.low.wp.grass <- polyploid.low.wp.i
    polyploid.high.wp.grass <- polyploid.high.wp.i
    
    base.op.grass <- base.op.i
    polyploid.low.op.grass <- polyploid.low.op.i
    polyploid.high.op.grass <- polyploid.high.op.i
  }
  if(life.forms[i]=="forb"){
    base.wp.forb <- base.wp.i
    polyploid.low.wp.forb <- polyploid.low.wp.i
    polyploid.high.wp.forb <- polyploid.high.wp.i
    
    base.op.forb <- base.op.i
    polyploid.low.op.forb <- polyploid.low.op.i
    polyploid.high.op.forb <- polyploid.high.op.i
    
  }
  if(life.forms[i]=="all"){
    base.wp.all <- base.wp.i
    polyploid.low.wp.all <- polyploid.low.wp.i
    polyploid.high.wp.all <- polyploid.high.wp.i
    
    base.op.all <- base.op.i
    polyploid.low.op.all <- polyploid.low.op.i
    polyploid.high.op.all <- polyploid.high.op.i
    
  }
}

life.form.names <- c("Woody", "Grasses", "Forbs", "All")
downsizing.names <- c("Low","Medium", "High")
base.cols <- c("#D55E00", "#F0E442", "#56B4E9", "grey")
polyploid.cols <- c("#660000", "#CC9900", "#003366", "black")
#polyploid.level.cols <- c("#9e9ac8","#756bb1","#54278f")
nz.range.base <- list(base.wp.woody, base.wp.grass, base.wp.forb, base.wp.all)
nz.range.polyploid.high <- list(polyploid.high.wp.woody, polyploid.high.wp.grass, polyploid.high.wp.forb, polyploid.high.wp.all)
nz.range.polyploid.low <- list(polyploid.low.wp.woody, polyploid.low.wp.grass, polyploid.low.wp.forb, polyploid.low.wp.all)
resampled.range.base <- list(base.op.woody, base.op.grass, base.op.forb, base.op.all)
resampled.range.polyploid.high <- list(polyploid.high.op.woody, polyploid.high.op.grass, polyploid.high.op.forb, polyploid.high.op.all)
resampled.range.polyploid.low <- list(polyploid.low.op.woody, polyploid.low.op.grass, polyploid.low.op.forb, polyploid.low.op.all)

#-------------Projected   

pdf(file="Density TTR base and life form low and high downsizing.pdf",width=8,height=6)
par(mfrow=c(2,2))

### resampled niche for each  life form ###
for(i in length(life.forms):1){
  if(i==4){par(mar=c(2, 5.1, 3.1, 0))}
  if(i==3){par(mar=c(2, 2, 3.1, 3.1))}
  if(i==2){par(mar=c(5.1, 5.1, 0, 0))}
  if(i==1){par(mar=c(5.1,2, 0,3.1))}
  
  
  plot(1,1,xlim=c(1, 150000),ylim=c(-0.000005,0.00007),type="n",axes=F, xlab="",ylab="")
  axis(2, at=c(0, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5))
  axis(1)
 
  #credible intervals
  arrows(y0=-0.000002, x0=as.numeric(bayesglm.res.df$credible.interval.lower[bayesglm.res.df$life.form==life.forms[i]&bayesglm.res.df$niche.type=="resample"&bayesglm.res.df$downsizing.level=="low"]),
         y1=-0.000002, x1=as.numeric(bayesglm.res.df$credible.interval.upper[bayesglm.res.df$life.form==life.forms[i]&bayesglm.res.df$niche.type=="resample"&bayesglm.res.df$downsizing.level=="low"]),
         code=3, angle=90, length=0.05, col=base.cols[i], lwd = 1.5)
  
  arrows(y0=-0.000004, x0=as.numeric(bayesglm.res.df$credible.interval.lower[bayesglm.res.df$life.form==life.forms[i]&bayesglm.res.df$niche.type=="resample"&bayesglm.res.df$downsizing.level=="high"]),
         y1=-0.000004, x1=as.numeric(bayesglm.res.df$credible.interval.upper[bayesglm.res.df$life.form==life.forms[i]&bayesglm.res.df$niche.type=="resample"&bayesglm.res.df$downsizing.level=="high"]),
         code=3, angle=90, length=0.05, col=polyploid.cols[i], lwd=1.5)
  
 lines(resampled.range.polyploid.low[[i]],col=base.cols[i],lwd=4)
  lines(resampled.range.polyploid.high[[i]],col=polyploid.cols[i],lwd=4)
  text(1000, 6e-5,life.form.names[i], font = 2, pos=4)
  
  text(0, 6.5e-5,paste(letters[4:1][i], ")", sep=""))
  
  legend(1000, 6e-05,fill=c(base.cols[i], polyploid.cols[i]),
         legend=c("Low", "High"),
         bty="n",hor=F) 
  
  if(i %in% c(2:1)){
    mtext("Projected resampled physiological range size\n(no. of 1km grid cells)",side=1, line=3)
    
  }
  if(i %in% c(4,2)){
    mtext("Probability density", side=2, line=3)
    
  }
}


dev.off()

##################################
###### 7. phylo PCA analysis######
##################################
# phylogenetic PCA of niche traits (TTR) by degree of genome downsizing and habit
#subset phylo and data to species with TTR data
downsizing.ttr <- downsizing.df[complete.cases(downsizing.df$resample.count),]
downsizing.ttr <- downsizing.ttr[!downsizing.ttr$Species=="squalida",] #not tip for Leptinella squalida
rownames(downsizing.ttr) <- downsizing.ttr$binomial
tr.ttr <- keep.tip(tr, rownames(downsizing.ttr))
downsizing.ttr <- downsizing.ttr[downsizing.ttr$binomial%in%tr.ttr$tip.label,]

# add in TTR trait differences
#calculate trait differences
# differences are base-polyploid
ttr.vars <- colnames(downsizing.ttr)[20:45]

for(i in 1:length(ttr.vars)){
  cols.i <- grep(ttr.vars[i], colnames(downsizing.ttr))
  new.col.no <- ncol(downsizing.ttr)+1
  downsizing.ttr[,new.col.no] <- downsizing.ttr[,cols.i[2]]-downsizing.ttr[,cols.i[1]]
  colnames(downsizing.ttr)[new.col.no] <- paste(ttr.vars[i], ".difference", sep="")
}
# subset TTR variables
ttr.pars <- downsizing.ttr[,84:107]
colnames(ttr.pars) <- gsub(".back.difference", "", colnames(ttr.pars))
pca.ttr <- phyl.pca(tr.ttr, ttr.pars)

for(i in 1:ncol(ttr.pars)){
  var.explained.i <- pca.ttr$Eval[i,i]/sum(pca.ttr$Eval) #variance explained by each PC
  if(i==1){var.explained <- var.explained.i }else{var.explained <- c(var.explained, var.explained.i)}
  
}
#scree plot
pdf("TTR trait differences pcs scree plot.pdf")
plot(var.explained, ylab="Explained variance", xlab="Principal component", pch=19)
lines(1:24, var.explained)
dev.off()
#go for 8 principal components


##### plotting PCA ############
#define colours
library(grDevices)
border.cols <- c("#660000", "#CC9900", "#003366")
habit.cols <- downsizing.ttr$Life.form
habit.cols <- gsub("woody", "#D55E00", habit.cols)
habit.cols <- gsub("grass", "#F0E442", habit.cols)
habit.cols <- gsub("forb", "#56B4E9", habit.cols)
bg.cols <- habit.cols
bg.cols <- gsub("#D55E00","#660000",bg.cols)
bg.cols <- gsub("#F0E442","#CC9900",bg.cols)
bg.cols <- gsub("#56B4E9","#003366",bg.cols)
habit.shape <- downsizing.ttr$Life.form
habit.shape <- gsub("woody", "24", habit.shape)
habit.shape <- gsub("grass", "21", habit.shape)
habit.shape <- gsub("forb", "22", habit.shape)

color.ramp<- colorRamp(colors = c("#046C9A", "white", "#F2300F"), space="rgb", bias=1)
gd.cols <-rgb(color.ramp((downsizing.ttr$downsizing.percent/200)+0.5), maxColorValue = 255) 
gd.cols2 <-rgb(color.ramp((downsizing.ttr$downsizing.percent/140)+0.5), maxColorValue = 255) 

gd.levels.cols <- rgb(color.ramp((seq(from=-10, to=60, length=5)/200)+0.5), maxColorValue = 255) 
gd.levels.cols <- rgb(color.ramp((seq(from=-10, to=60, length=5)/140)+0.5), maxColorValue = 255) 

color.ramp2<- colorRamp(colors = c("white", "#F2300F"), space="rgb")
niche.size.cols <-rgb(color.ramp2(downsizing.ttr$resample.count/max(downsizing.ttr$resample.count)), maxColorValue = 255) 
niche.levels.cols <- rgb(color.ramp2(seq(from=25000, to=175000, length=4)/max(downsizing.ttr$resample.count)), maxColorValue = 255) 


###ellipses for GD
# calculate correlations betw pc1 and pc2 for each downsizing range
tab <- matrix(c(pca.ttr$S[,1], pca.ttr$S[,2]), ncol=2)
cneg10 <- cor(tab[downsizing.ttr$downsizing.percent<0,])
c10 <- cor(tab[downsizing.ttr$downsizing.percent>0&downsizing.ttr$downsizing.percent<20,])
c30 <- cor(tab[downsizing.ttr$downsizing.percent>20&downsizing.ttr$downsizing.percent<40,])
c50 <- cor(tab[downsizing.ttr$downsizing.percent>40,])

ellipse.neg10 <- ellipse(cneg10, centre=colMeans(tab[downsizing.ttr$downsizing.percent<0,]), 
                         level=0.95, scale=c(sd(tab[downsizing.ttr$downsizing.percent<0,1]),
                                             sd(tab[downsizing.ttr$downsizing.percent<0,2])))
ellipse.10 <- ellipse(c10, centre=colMeans(tab[downsizing.ttr$downsizing.percent>0&downsizing.ttr$downsizing.percent<20,]), 
                      level=0.95, scale=c(sd(tab[downsizing.ttr$downsizing.percent>0&downsizing.ttr$downsizing.percent<20,1]),
                                          sd(tab[downsizing.ttr$downsizing.percent>0&downsizing.ttr$downsizing.percent<20,2])))
ellipse.30 <- ellipse(c30, centre=colMeans(tab[downsizing.ttr$downsizing.percent>20&downsizing.ttr$downsizing.percent<40,]), 
                      level=0.95, scale=c(sd(tab[downsizing.ttr$downsizing.percent>20&downsizing.ttr$downsizing.percent<40,1]),
                                          sd(tab[downsizing.ttr$downsizing.percent>20&downsizing.ttr$downsizing.percent<40,2])))
ellipse.50 <- ellipse(c50, centre=colMeans(tab[downsizing.ttr$downsizing.percent>40,]), 
                      level=0.95, scale=c(sd(tab[downsizing.ttr$downsizing.percent>40,1]),
                                          sd(tab[downsizing.ttr$downsizing.percent>40,2])))

pdf(file="phyl.pca TTR by downsizing level.pdf")
#downsizing
plot(1, xlab="PC1", ylab="PC2", xlim=c(-320, 300), ylim=c(-250,320), type="n", asp=1)
polygon(ellipse.neg10, col=NA, border=gd.levels.cols[1], lwd=3)
polygon(ellipse.10, col=NA, border=gd.levels.cols[2], lwd=3)
polygon(ellipse.30, col=NA, border=gd.levels.cols[3], lwd=3)
polygon(ellipse.50, col=NA, border=gd.levels.cols[4], lwd=3)

#arrows(x0=rep(0,nrow(pca.ttr$L)), y0=rep(0,nrow(pca.ttr$L)),x1=300*(pca.ttr$L[,1]), y1=300*(pca.ttr$L[,2]), code=2, length = 0.1)
points(pca.ttr$S[,1], pca.ttr$S[,2], pch=as.numeric(habit.shape), bg=gd.cols, cex=1.25 )
#text(300*(pca.ttr$L[,1]), 300*(pca.ttr$L[,2]), rownames(pca.ttr$L),cex=0.7, pos=c(2,2,1,1,1,1,1,2,3,1,1,2,4,4,1,4,4,1,1,2,2,2,4,3))

#legend
legend(x = -340, y = 310,
       legend = rep(NA, 25),
       fill = colorRampPalette(colors = c("#F2300F", "white", "#046C9A"))(25),
       border = NA,
       bty="n",
       y.intersp = 0.1,
       cex = 1, text.font = 2)
points(rep(-315,3), c(170,150,130), pch=c(22,21,24))
text(c(-340, -340), c(320,190), c("Genome downsizing (%)", "Life form"), pos=4, cex=0.8)
text(c(-315, -315, -315), c(300, 267.5,235), c("70", "0", "-70"), pos=4, cex=0.75)
text(c(-315, -315, -315), c(170, 150,130), c("forb", "grass", "woody"), pos=4, cex=0.75)

dev.off()

pdf(file="phyl.pca biplot TTR.pdf")

plot(pca.ttr$Evec, type="n", ylim=c(-0.55,0.55), xlim=c(-0.55,0.55))
arrows(x0=rep(0,nrow(pca.ttr$Evec)), y0=rep(0,nrow(pca.ttr$Evec)), x1=pca.ttr$Evec[,1], y1=pca.ttr$Evec[,2], length = 0.1, code=2)
text(pca.ttr$Evec[,1], pca.ttr$Evec[,2], rownames(pca.ttr$Evec), pos=c(3,4,1,1,1,4,2,3,4,1,4,1,2,2,3,4,2,2,2,2,2,2,2,4))

dev.off()

## plot GD by components ##
pca.df <- cbind(downsizing.ttr[,c(5,10,81)], pca.ttr$S)
ttr.cd <- comparative.data(tr, pca.df, names.col = "binomial")

pca.gd.results.df <- data.frame("Principal component"=character(), "slope"=numeric(),"F"=numeric(), "df1"=numeric(), "df2"=numeric(),"p.value"=numeric(), "R.sq"=numeric())

## without habit ##
pdf("Principal components and genomic downsizing plots.pdf", height=15, width=7)
par(mfrow=c(4,2))
for(i in 1:8){
  ttr.cd.i <- comparative.data(tr, pca.df[,c(1:3,(i+3))], names.col = "binomial")
  colnames(ttr.cd.i$data)[3] <- "PC"
  pgls.i <- pgls(downsizing.percent~PC, data=ttr.cd.i, lambda = "ML")
  summary.i <- summary(pgls.i)
  
  #if(i %in% 1:8){
    plot(pca.ttr$S[,i], downsizing.ttr$downsizing.percent, pch=19, xlab=paste("PC", i, sep=""), ylab="Genome downsizing (%)", xlim=c(-300, 300), ylim=c(-50, 100))
   if(summary.i$coefficients[2,4]<0.05){abline(pgls.i) } 
    text(-300,-20, paste("R2 = ", round(summary.i$r.squared,3)), pos=4)
    text(-300,-30, paste("p value = ", round(summary.i$coefficients[2,4],3)), pos=4)
    
  #   }
  
  new.row <- c(i,summary.i$coefficients[2,1],summary.i$fstatistic, summary.i$coefficients[2,4],summary.i$r.squared)
  pca.gd.results.df[i,] <- new.row
  
  rm(new.row)
}
dev.off()

pca.gd.results.df[,2:7] <- round(pca.gd.results.df[,2:7], 3)
write.csv(pca.gd.results.df, file="PCA pgls results no habit.csv")


## with habit ##
pca.gd.habit.results.df <- data.frame("Principal component"=character(), "F"=numeric(), "df1"=numeric(), "df2"=numeric(),"R2"=numeric(),"p.value.PC"=numeric(), "p.value.grass"=numeric(), "p.value.woody"=numeric(),"p.value.interaction.grass"=numeric(), "p.value.interaction.woody"=numeric(),
                                      "effect.size.PC"=numeric(), "effect.size.grass"=numeric(), "effect.size.woody"=numeric(),"effect.size.interaction.grass"=numeric(), "effect.size.interaction.woody"=numeric())

for(i in 1:8){
    ttr.cd.i <- comparative.data(tr, pca.df[,c(1:3,(i+3))], names.col = "binomial")
  colnames(ttr.cd.i$data)[3] <- "PC"
  pgls.i <- pgls(downsizing.percent~PC*Life.form, data=ttr.cd.i, lambda = "ML")
  summary.i <- summary(pgls.i)
  new.row <- c(i,summary.i$fstatistic, summary.i$r.squared,summary.i$coefficients[2:6,4],summary.i$coefficients[2:6,2])
  pca.gd.habit.results.df[i,] <- new.row
  
  rm(new.row)
  }

#adjust p-values for multiple comparisons
pca.gd.habit.results.df$adjusted.p.value.PC <- rep(NA, nrow(pca.gd.habit.results.df))
pca.gd.habit.results.df$adjusted.p.value.PC[1:8]<- p.adjust(pca.gd.habit.results.df$p.value.PC[1:8], method ="holm" )

pca.gd.habit.results.df$adjusted.p.value.grass <- rep(NA, nrow(pca.gd.habit.results.df))
pca.gd.habit.results.df$adjusted.p.value.grass[1:8]<- p.adjust(pca.gd.habit.results.df$p.value.grass[1:8], method ="holm" )

pca.gd.habit.results.df$adjusted.p.value.woody <- rep(NA, nrow(pca.gd.habit.results.df))
pca.gd.habit.results.df$adjusted.p.value.woody[1:8]<- p.adjust(pca.gd.habit.results.df$p.value.woody[1:8], method ="holm" )

pca.gd.habit.results.df$adjusted.p.value.interaction.grass <- rep(NA, nrow(pca.gd.habit.results.df))
pca.gd.habit.results.df$adjusted.p.value.interaction.grass[1:8]<- p.adjust(pca.gd.habit.results.df$p.value.interaction.grass[1:8], method ="holm" )
pca.gd.habit.results.df$adjusted.p.value.interaction.woody <- rep(NA, nrow(pca.gd.habit.results.df))
pca.gd.habit.results.df$adjusted.p.value.interaction.woody[1:8] <- p.adjust(pca.gd.habit.results.df$p.value.interaction.woody[1:8], method ="holm" )

pca.gd.habit.results.df[,2:15] <- round(pca.gd.habit.results.df[,2:15], 3)

write.csv(t(pca.gd.habit.results.df), file="PCA pgls results.csv")
