
###Smaller version of the earlier "AssembleMasterDNAtable" script; omitting DUP names and anything not passing quality control

library(data.table)
library(taxize)

#set working dir
setwd("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/TideEffectHoodCanal_Gallego/Analysis/blasting")


#read in the taxonomic annotation (modified output from MEGAN/BLAST)
taxonomy1=read.csv("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/TideEffectHoodCanal_Gallego/Analysis/blasting/BLASTed_COI_tides_OTUs_20170930_min224_Genus-ex.txt", header=F)

#assemble master data frame
master= taxonomy1
	names(master)=c("OTU_name_swarm", "OTU_abundance","OTU_taxon_MEGAN")


# # ######if desired, lookup taxonomic lineages for unique family names (or names of whatever rank)#######

familyNames=unique(master$OTU_taxon_MEGAN)
	familyNames = familyNames[-which(familyNames%in%c("No hits", "Not assigned", NA))]
	familyNames <- gsub("environmental samples <","", familyNames)
	familyNames <- gsub(">","", familyNames)

# #note that if there are multiple valid ncbi taxa matching a given name, this loop will require user input to choose from among them; don't walk away.
# #note also that if the whole lineage is an NA (i.e., if the name searched isn't in NCBI database), the loop after this one will fail -- you won't be able to get the taxonomy out easily.  Suggest changing the annotations in the taxonomy1 file to something sensible
lineages=list(NA)
for (i in 1:length(familyNames)){
		tempx=classification(familyNames[i], db="ncbi")
	lineages[i]=ifelse(isTRUE(tempx), NA, tempx)
}
 saveRDS(lineages, "Familylineages_20171012.rds")	
# #NOTE: the ITIS database uses capital letters ("Phylum"), while NCBI uses lowercase ("phylum")

#lineages=readRDS("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Data/blast_20151116_1714/Familylineages.rds")

 outtable =data.frame(rep(NA, length(familyNames)), rep(NA, length(familyNames)),rep(NA, length(familyNames)),rep(NA, length(familyNames)),rep(NA, length(familyNames)),rep(NA, length(familyNames)),rep(NA, length(familyNames)))

 for (i in which(!is.na(lineages))){
	 outtable[i,] <- lineages[[i]][[1]][match(c("kingdom", "phylum", "class", "order", "family","genus", "species"), lineages[[i]][[2]])]
 }
 names(outtable)=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#augment taxonomy 
	outtable[which(outtable $Class=="Dinophyceae"),2]<-"Dinoflagellata"
	outtable[which(outtable $Phylum=="Bacillariophyta"),1]<-"Diatoms"
	outtable[which(outtable $Phylum=="Dinoflagellata"),1]<-"Dinoflagellates"
	outtable[which(outtable $Phylum=="Proteobacteria"),1]<-"Bacteria"
	outtable[which(outtable $Phylum=="Phaeophyceae"),1]<-"Heterokonta"
	outtable[which(outtable $Class=="Oomycetes"),1]<-"Heterokonta"	
	outtable[which(outtable $Class=="Florideophyceae"),2]<-"Rhodophyta"	
	outtable[which(outtable $Class=="Florideophyceae"),1]<-"Rhodophyta"
	outtable[which(outtable $Class=="Spirotrichea"),1]<-"Protista"
	outtable[which(outtable $Class=="Spirotrichea"),2]<-"Cliophora"

# #if desired, manually take a look at those lineages that failed...
# which(is.na(lineages))
# n=which(is.na(lineages))[3]
	# tempx=classification(familyNames[n], db="ncbi")
	# lineages[n]=ifelse(isTRUE(tempx), NA, tempx)
# #then re-run outtable assembly above

# #write small lookup table
 write.csv(outtable, "LineageLookupTable.csv")

#outtable=outtable[order(outtable[,1], outtable[,2], outtable[,3], outtable[,4]),] #sort alphabetically
#	outtable=outtable[-which(is.na(outtable[,1])),]

#outtable=read.csv("/Users/rpk/GoogleDrive/NWFSC_iMac_blast/blast_20151204/LineageLookupTable.csv", row.names=1)


stripTaxonomy<-function(x, rank="Genus") {outtable[match(master[x,3], unique(familyNames)),rank]} 


master$OTU_taxon_species <- stripTaxonomy(1:nrow(master), "Species")
master$OTU_taxon_genus <- stripTaxonomy(1:nrow(master), "Genus")
master$OTU_taxon_family <- stripTaxonomy(1:nrow(master), "Family")
master$OTU_taxon_order <- stripTaxonomy(1:nrow(master), "Order")
master$OTU_taxon_class <- stripTaxonomy(1:nrow(master), "Class")
master$OTU_taxon_phylum <- stripTaxonomy(1:nrow(master), "Phylum")
master$OTU_taxon_kingdom <- stripTaxonomy(1:nrow(master), "Kingdom")

##add lineage information to the master table  (this is inefficient, in that it writes many duplicate names to the master file, but might be handy if looking at an individual OTU)

# #infer lineages from existing Genus-level information
# master=data.frame(master, MeganAssignmentSpecies , MeganAssignmentGenus) ; names(master)[4:5]=c("OTU_taxon_species","OTU_taxon_genus")
	# masterfamilies=outtable[,5][ifelse(is.na(master$OTU_taxon_genus), NA, match(master$OTU_taxon_genus,outtable$Genus))]
	# masterphyla=outtable[,2][ifelse(is.na(master$OTU_taxon_genus), NA, match(master$OTU_taxon_genus,outtable$Genus))] 		
	# masterclasses=outtable[,3][ifelse(is.na(master$OTU_taxon_genus), NA, match(master$OTU_taxon_genus,outtable$Genus))] 
	# masterorders=outtable[,4][ifelse(is.na(master$OTU_taxon_genus), NA, match(master$OTU_taxon_genus,outtable$Genus))] 
	# masterkingdoms=outtable[,1][ifelse(is.na(master$OTU_taxon_genus), NA, match(master$OTU_taxon_genus,outtable$Genus))] 

# #backfill info from higher taxa, if those were the best Megan could do
# for (i in 1:nrow(master)){
	# if(!is.na(MeganAssignmentFamily[i])) masterfamilies[i] <- outtable[,5][match(MeganAssignmentFamily[i], outtable[,5])]
	# if(!is.na(MeganAssignmentOrder[i])) masterorders[i] <- outtable[,4][match(MeganAssignmentOrder[i], outtable[,4])]
	# if(!is.na(MeganAssignmentClass[i])) masterclasses[i] <- outtable[,3][match(MeganAssignmentClass[i], outtable[,3])]
	# if(!is.na(MeganAssignmentPhylum[i])) masterphyla[i] <- outtable[,2][match(MeganAssignmentPhylum[i], outtable[,2])]
	# if(!is.na(MeganAssignmentKingdom[i])) masterphyla[i] <- outtable[,1][match(MeganAssignmentKingdom[i], outtable[,1])]
# }


		# master=data.frame(master, masterfamilies, masterorders, masterclasses, masterphyla, masterkingdoms)
			# names(master)[5:9]<-c("OTU_taxon_family","OTU_taxon_order", "OTU_taxon_class", "OTU_taxon_phylum", "OTU_taxon_kingdom")

##further backfilling ; 	#where you have Order info, fill in Class and Phylum

for (i in which(is.na(master$OTU_taxon_family)&!is.na(master$OTU_taxon_order))){
	master$OTU_taxon_class[i]<- outtable[,3][match(master$OTU_taxon_order[i], outtable[,4])]
	master$OTU_taxon_phylum[i]<- outtable[,2][match(master$OTU_taxon_order[i], outtable[,4])]
	master$OTU_taxon_kingdom[i]<- outtable[,1][match(master$OTU_taxon_order[i], outtable[,4])]
}

	#where you have only Class info, fill in Phylum
for (i in which(is.na(master$OTU_taxon_family)&is.na(master$OTU_taxon_order)&!is.na(master$OTU_taxon_class))){
	master$OTU_taxon_phylum[i]<- outtable[,2][match(master$OTU_taxon_class[i], outtable[,2])]
	master$OTU_taxon_kingdom[i]<- outtable[,1][match(master$OTU_taxon_class[i], outtable[,2])]
}

#for safety
master[is.na(master$OTU_taxon_MEGAN),4:8]<-NA

#add in number of dups that go into each OTU, useful for later calculations
# master=data.frame(as.numeric(nDUPs_in_OTU[match(row.names(master), names(nDUPs_in_OTU))]), master)
	# names(master)[1]="N_DUPs"


#manual fix for taxa that are outside of the designations of class, order, etc
oddNames=unique(master[which(!is.na(master$OTU_taxon_MEGAN)&is.na(master$OTU_taxon_phylum)),]$OTU_taxon_MEGAN)
	oddNames_class=NA
	oddNames_phylum=NA
	oddNames_kingdom=NA
for (i in 1:length(oddNames)){
oddNames_class[i]<- ifelse(length(lineages[[min(grep(oddNames[i],lineages))]][,1][lineages[[min(grep(oddNames[i],lineages))]][,2]=="class"])==0, NA, lineages[[min(grep(oddNames[i],lineages))]][,1][lineages[[min(grep(oddNames[i],lineages))]][,2]=="class"])
oddNames_phylum[i]<- ifelse(length(lineages[[min(grep(oddNames[i],lineages))]][,1][lineages[[min(grep(oddNames[i],lineages))]][,2]=="phylum"])==0, NA, lineages[[min(grep(oddNames[i],lineages))]][,1][lineages[[min(grep(oddNames[i],lineages))]][,2]=="phylum"])
oddNames_kingdom[i]<- ifelse(length(lineages[[min(grep(oddNames[i],lineages))]][,1][lineages[[min(grep(oddNames[i],lineages))]][,2]=="kingdom"])==0, NA, lineages[[min(grep(oddNames[i],lineages))]][,1][lineages[[min(grep(oddNames[i],lineages))]][,2]=="kingdom"])
}

for (i in 1:length(oddNames)){
master[master$OTU_taxon_MEGAN%in%oddNames[i],"OTU_taxon_class"]<-as.character(oddNames_class[i])
master[master$OTU_taxon_MEGAN%in%oddNames[i],"OTU_taxon_phylum"]<-as.character(oddNames_phylum[i])
master[master$OTU_taxon_MEGAN%in%oddNames[i],"OTU_taxon_kingdom"]<-as.character(oddNames_kingdom[i])
}

master[which(master$OTU_taxon_phylum=="Arthropoda"),"OTU_taxon_kingdom"]<-"Metazoa"
#master[which(master$OTU_taxon_MEGAN =="Alphaproteobacteria"),8:9]<-c("Proteobacteria","Bacteria")
#master[which(master$OTU_taxon_MEGAN =="Flavobacteriaceae"),9]<-"Bacteria"

#higher taxa becomes problematic for some of these lineages
#head(master[!is.na(master$OTU_taxon_family)&is.na(master$OTU_taxon_phylum),],20)


#write file out
write.csv(master, "masterDNAtable_20171012.csv", row.names=F, quote=F)


taxonCounts=as.data.frame(matrix(NA, nrow=length(unique(master$OTU_taxon_phylum)), ncol=7))
	names(taxonCounts)=c("Kingdom", "Phylum","Classes", "Orders", "Families", "Genera", "OtherRank")
for (i in 1:length(unique(master$OTU_taxon_phylum))){
	foc=unique(master$OTU_taxon_phylum)[i]
	taxonCounts[i,1]=outtable[!is.na(outtable$Phylum)&outtable$Phylum==foc,1][1]
	taxonCounts[i,2]=as.character(foc)
		taxonCounts[i,3]= length(unique(master[master$OTU_taxon_phylum == foc,7]))-sum(is.na(unique(master[master$OTU_taxon_phylum == foc,7])))
		taxonCounts[i,4]= length(unique(master[master$OTU_taxon_phylum == foc,6]))-sum(is.na(unique(master[master$OTU_taxon_phylum == foc,6])))
		taxonCounts[i,5]= length(unique(master[master$OTU_taxon_phylum == foc,5]))-sum(is.na(unique(master[master$OTU_taxon_phylum == foc,5])))
		taxonCounts[i,6]= length(unique(master[master$OTU_taxon_phylum == foc,4]))-sum(is.na(unique(master[master$OTU_taxon_phylum == foc,4])))
				taxNames=unique(master[master$OTU_taxon_phylum == foc,3]) ; taxNames=taxNames[!is.na(taxNames)]
				classified=unique(unlist(master[master$OTU_taxon_phylum == foc,4:7])) ; classified<-classified[!is.na(classified)]
		taxonCounts[i,7] = length(taxNames[is.na(match(taxNames, classified))])
	}
taxonCounts[taxonCounts[,6]==1,c(4:5)]<-1
#taxonCounts[taxonCounts[,5]==0,5]<-NA
taxonCounts<-taxonCounts[!is.na(taxonCounts[,1]),]
taxonCounts<-taxonCounts[order(taxonCounts$Kingdom, taxonCounts$Phylum),]

write.csv(taxonCounts, "TaxonCountsSummaryTable_20171012.csv")



