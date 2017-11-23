#species occupancy modeling
#what's the probability of seeing some number of counts of eDNA when the species is not, in fact, present?  i.e., what's the false-positive rate?
#see Lahoz-Monfort 2015 Mol Ecol Res

setwd("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/TideEffectHoodCanal_Gallego/Analysis/OccupancyModeling")

			### Required packages
			library(unmarked)  # Maximum likelihood estimation for SODMs 
			library(vegan)
			library(jagsUI)    # Bayesian analysis with JAGS
			library(data.table)
			library(digest)
			library(tidyr)
			library(dplyr)
			library(googlesheets)
			library(Biostrings) 
			# This function is later used for printing simulation results on the screen
			printsummres<-function(thetahat,thename="estimated parameter"){
			  cat(paste0("\n",thename,": mean = ",round(mean(thetahat),2),
			              ", 2.5% = ",round(quantile(thetahat,prob=0.025),2),
			              ", 97.5% = ",round(quantile(thetahat,prob=0.975),2)))
			}
		
			#function for geometric means, from https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
			gm_mean = function(x, na.rm=TRUE){
		  	exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
			}

			getSample<-function(x){strsplit(x, "_")[[1]][1]}
			
			ceilingMean<-function(x) ceiling(mean(x))
			
			Col2RN<-function(df, x){row.names(df)<-df[,x]; df<-df[-x]; return(df)}

			
##Load data, etc			


				#MiSeq full run metadata
				  	meta=read.csv("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/TideEffectHoodCanal_Gallego/Data/COI/FullRun/metadata_wTides.csv")
				  	sampleInfo<-as.data.frame(gs_read(gs_title("Step1_Samples"))) #field data for covariates
				  		names(sampleInfo)[1]<-"original_sample"
				  		meta<-merge(meta, sampleInfo, by="original_sample", all=T)
				  		meta<-meta[1:99,]
				  	lib_tags= paste0("ID1=",meta$pri_index_name, ";ID2A=",meta$tag_sequence, ";ID2B=",as.character(reverseComplement(DNAStringSet(meta$tag_sequence))))  
				  	  meta<-meta[-4,]
				  	lib_tags<-lib_tags[-4]  #one missing, for some reason
				      meta$lib_tags<-lib_tags #add lib tags to metadata
				      meta$run<-"20170727_MiSeq"
					
				
				###get metadata from nano run, and merge
				    metaNano=read.csv("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/TideEffectHoodCanal_Gallego/Data/COI/Nano_banzai_out_20170719_1420/metadata_correctedRamon_20170926.csv")
				    sampleInfo<-as.data.frame(gs_read(gs_title("Step1_Samples"))) #field data for covariates
				    names(sampleInfo)[1]<-"original_sample"
				    metaNano <-merge(metaNano, sampleInfo, by="original_sample", all=T)
				    metaNano<-metaNano[!is.na(metaNano$sample_name),]
				    levels(metaNano$tag_sequence)<-c(levels(metaNano$tag_sequence), "noTag")
				    metaNano$tag_sequence[metaNano$tag_number=="No_Tag"]<-"noTag"
				      lib_tags_nano= paste0("ID1=", metaNano $pri_index_seq, ";", metaNano $tag_sequence)  
				    metaNano$lib_tags<-lib_tags_nano
				    metaNano$run<-"20170710_Nano"
				
				#create merged metadatafile
				allMeta<-merge(meta, metaNano, all=T)
				    allMeta$sample_name<-as.character(allMeta$sample_name)
				    #need to disambiguate sample names across runs and homogenize naming conventions
				    allMeta$sample_name[99:106]<-paste0("TW20170311", gsub("TW11", "",allMeta$sample_name[99:106]))
				    allMeta$sample_name[107:110]<-paste0("LL20170312", gsub("LL12", "",allMeta$sample_name[107:110]))
				    allMeta$sample_name[111]<-"Ostrich4"
				    
				    #interpolate/harmonize salinity measurements [we measured with two different instruments -- a Hannah multi-probe (here, called "YSI", because that was shorter), and a hand-held analog refractometer.]
				    mod<-lm(Salinity_YSI~Salinity_refractometer, data=allMeta)
				    mod2<-lm(Salinity_refractometer ~Salinity_YSI, data=allMeta)
				    	missingYSI<-which(is.na(allMeta$Salinity_YSI) &  !is.na(allMeta$Salinity_refractometer))
				  	missingRefract<-which(is.na(allMeta$Salinity_refractometer) &  !is.na(allMeta$Salinity_YSI))  
				    allMeta$Salinity_YSI[missingYSI]<-round((allMeta$Salinity_refractometer[missingYSI]*mod[[1]][2])+mod[[1]][1],1)
				    allMeta$Salinity_refractometer[missingRefract]<-round((allMeta$Salinity_YSI[missingRefract]*mod2[[1]][2])+mod2[[1]][1],1)
				
				#format date/time in metadata
				 allMeta$DateTime<- as.POSIXct(strptime(paste(allMeta$Date, allMeta$Time), format="%Y%m%d %H:%M:%S"))
				    
				    
				    
				    
				#######################################################
				
				##########LOAD OTU Data  ##############################
				#OTU data
				allOTUs<-
					read.table("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/TideEffectHoodCanal_Gallego/Data/COI/CombinedNanoFullRun/OTUs_swarm/OTU.map") %>%
					filter(V2%in%allMeta$lib_tags) 
				allOTUs$V2<-allMeta $sample_name[match(allOTUs$V2, allMeta$lib_tags)]
				allOTUs<-spread(allOTUs, key=V2, value=V3, fill=0)
					allOTUs <- Col2RN(allOTUs, 1)

	#to harmonize w earlier code			
	OTUs=allOTUs  ; rm(allOTUs)
	meta=allMeta
	
	#filter singletons, if desired
	OTUs<-OTUs[rowSums(OTUs)>1,]

		bottles= substr(meta$sample_name,1,11) #make human-readable and uniform within bottle sample
		insuff<-names(which(table(bottles)<3))
		tooMany<-names(which(table(bottles)>3)) #names of samples with more than three replicates; use to subsample so that there are exactly 3

	REPLICATEDmeta<-meta[-which(bottles%in%insuff),] #NOTE: assumes meta and OTU table are in the same order; this removes ALL instances of insuff in data and metadata
	REPLICATEDmeta<-REPLICATEDmeta[-match(tooMany, substr(REPLICATEDmeta$sample_name, 1, 11)),] #just removes first instance of tooMany (knocking samples w 4 replicates down to three... TODO: generalize this)
				
	REPLICATED <-OTUs[, REPLICATEDmeta$sample_name] 
	
	Nreplicates<-3  #hard-coded; number of replicates per sample. TODO -- figure out a way to infer this from the data. could use table(colnames(OTUs))

# #rarefy and quick QC for read-depth
	# set.seed(49)
	# REPLICATEDmeta<-REPLICATEDmeta[-which(colSums(REPLICATED)<20000),]
	# REPLICATED<-REPLICATED[,-which(colSums(REPLICATED)<20000)]  #gets rid of very-low abundance CANON reads
		# REPLICATED<-t(rrarefy(t(REPLICATED), min(colSums(REPLICATED))))		#scale(REPLICATED, center=F, scale=colSums(REPLICATED))
		# REPLICATED<-REPLICATED[-which(rowSums(REPLICATED)<2),] #remove singletons

#create convenience variables
	Nsites=ncol(REPLICATED)/Nreplicates
#	Annotated<-intersect(row.names(REPLICATED), row.names(tax))
#	sampleDates<-as.Date(REPLICATEDmeta$DATE)
	sampleNames<-substr(colnames(REPLICATED),1,11) #create vector of sample names
	siteNames<-unique(sampleNames)  #substr(colnames(REPLICATED),1,2)
#	yearvector<-format(sampleDates, "%Y")
	dayvector<-substr(colnames(REPLICATED),9,10)

REPLICATED<-t(as.data.table(REPLICATED))
	colnames(REPLICATED)<-row.names(OTUs)
	ReplicatedMeans<-aggregate.data.frame(REPLICATED, by=list(sampleNames), ceilingMean)  #create df; mean number of reads per OTU per sample
		 ReplicatedMeans<-Col2RN(ReplicatedMeans, 1) ; ReplicatedMeans<-t(ReplicatedMeans) ; ReplicatedMeans<-ReplicatedMeans[, unique(sampleNames)]

# #rather than OTU-specific occupancy params, it might be good to have taxon-specific params. Here is a df that helps to do that.
# taxReadsTable<-REPLICATED[Annotated,]
	# tax<-tax[row.names(tax)%in%Annotated,]
	# taxReadsTable<-aggregate(taxReadsTable, by=list(tax[,2]), sum)
	# taxReadsMeans<-aggregate.data.frame(t(taxReadsTable[,2:ncol(taxReadsTable)]), by=list(rep(siteNames, each=3)), ceilingMean)  #create df; mean number of reads per taxon per site
		# row.names(taxReadsMeans)<-taxReadsMeans[,1]; taxReadsMeans <-taxReadsMeans[,-1] ; taxReadsMeans <-t(taxReadsMeans); taxReadsMeans <-taxReadsMeans[,siteNames]
		# row.names(taxReadsMeans)<-taxReadsTable[,1]



##############################################
##############################################
##############################################
##############################################

#needs A binary detection matrix where rows=sites, & columns = replicates.
myRealData = REPLICATED ; myRealData[myRealData>0]<-1 ; myRealData=t(myRealData)
	#order rows
	#myRealData= myRealData[order(row.names(myRealData)),]
	#remove OTUs that don't occur at least twice in this dataset
	myRealData= myRealData[,colSums(myRealData)>1]
		#filter ReplicatedMeans accordingly, to match myRealData
		#ReplicatedMeans<-ReplicatedMeans[row.names(ReplicatedMeans)%in%colnames(myRealData),]

#simultaneously create a data-output long-df capturing OTU name and model info, etc.
#this will have one row for each OTU at each of Nsites
	dataOUT<-data.frame(
        OTUname=rep(row.names(myRealData), each=Nsites), 
	      bottle=rep(siteNames, times=nrow(myRealData)),
	#date=rep(unique(sampleDates), times=ncol(myRealData)),
	Nrep=rep(Nreplicates, each=Nsites*nrow(myRealData)),
	MeanReads=as.vector(t(ReplicatedMeans))
	)


# #for taxon-specific occupancy [note: should make into a long/tidy df rather than this awk format]
# taxData = taxReadsTable[,-1] ; taxData[taxData>0]<-1 ; taxData =t(taxData); colnames(taxData)<-taxReadsTable[,1]
	# taxData = taxData[,colSums(taxData)>1]
		# #filter ReplicatedMeans accordingly, to match myRealData
		# taxReadsMeans <-taxReadsMeans[row.names(taxReadsMeans)%in%colnames(taxData),]

# #create analogous data for taxon-specific detection rates
# taxDataOut <-data.frame(Taxon=rep(colnames(taxData), each=Nsites), 
	# site=rep(siteNames, times=ncol(taxData)),
	# date=rep(unique(sampleDates), times=ncol(taxData)),
	# Nrep=rep(Nreplicates, each=Nsites*ncol(taxData)),
	# MeanReads=as.vector(t(taxReadsMeans))
	# )

#treating each transect as its own site ; reshape data for model fitting
OTUlist=list()
for (i in 1:nrow(myRealData)){
	OTU=c(myRealData[i,])
		dim(OTU)<-c(Nreplicates, ncol(myRealData)/Nreplicates)  #reshape ; note these params will depend upon the input data in question
	OTUlist[[i]]=t(OTU)
}
Ndetections<-unlist(lapply(OTUlist, rowSums)) #collapses matrix, to create vector of N detections per bottle of water
dataOUT<-data.frame(dataOUT, Ndetections)


# Taxlist=list()
# for (i in 1:ncol(taxData)){
	# TAX=c(taxData[,i])
		# dim(TAX)<-c(Nreplicates, nrow(taxData)/Nreplicates)  #reshape ; note these params will depend upon the input data in question
	# Taxlist[[i]]=t(TAX)
# }
# NtaxDetections<-unlist(lapply(Taxlist, rowSums))
# taxDataOut<-data.frame(taxDataOut, NtaxDetections)


										##############################################
										##############################################
										#The Actual Model
										##############################################
										##############################################
										##############################################
										
										####Bayesian version
										### Run first this part once to create the file with the JAGS model
										sink("RoyleLink_prior.txt")
										cat("model {
										    # Priors
										    psi ~ dunif(0,1)
										    p11 ~ dunif(0.01,1)
										    p10 ~ dunif(0.001, p10_max)
										    
										    # Likelihood 
										    for (i in 1:S){
										    z[i] ~ dbern(psi)
										    p[i] <- z[i]*p11 + (1-z[i])*p10
										    for (j in 1:K){
										    Y[i,j] ~ dbern(p[i])
										    }
										    }
										    } ",fill=TRUE)
										sink()
										
										
										model_Bayesian <- function(datalist,nOTUs=length(datalist), S=Nsites, K=Nreplicates, doprint=FALSE,p10_max=0.10,
										                                     ni=3000,nt=2,nc=1,nb=1000,myparallel=TRUE) {   
										  psihat<-p11hat<-p10hat<-rep(nOTUs)
										  modelSummaries<-list()
										  for(ii in 1: nOTUs){
										    #cat("\r", ii, "of", nOTUs,"   ")
										    #if (doprint) cat("\r", ii, "of", nOTUs,"   ")    
										    hh<-datalist[[ii]]
										    # fit the model    
										    jags.inits <-function()(list(psi=runif(1,0.05,0.95),p11=runif(1, 0.01,1),p10=runif(1,0.001,p10_max)))
										    jags.data  <-list(Y=hh,S=S,K=K,p10_max=p10_max)
										    jags.params<-c("psi","p11","p10")
										    model<-jags(data = jags.data, inits = jags.inits, parameters.to.save= jags.params, 
										                model.file= "RoyleLink_prior.txt", n.thin= nt, n.chains= nc, 
										                n.iter= ni, n.burnin = nb, parallel=myparallel)  #, working.directory= getwd()
											    # jpeg(paste0(format(Sys.time(), "%H_%M_%S"),"_ModelParamsPlot_",names(focalOTUlist)[ii],".jpg"))
											     		# plot(model)
											    	# dev.off()
										    # extract results (medians of the marginal posteriors)
										    psihat[ii] <- model$summary["psi","50%"]
										    p11hat[ii] <- model$summary["p11","50%"]
										    p10hat[ii] <- model$summary["p10","50%"]    
										    modelSummaries[[ii]]<-model$summary
										  }
										  if (doprint){
										    printsummres(psihat,thename="estimated psi")
										    printsummres(p11hat,thename="estimated p11")
										    printsummres(p10hat,thename="estimated p10")
										  }
										  saveRDS(modelSummaries, paste0(format(Sys.time(), "%H_%M_%S"),"_ModelSummaries.rds"))
											BayesResults<-list(psihat=psihat,p11hat=p11hat,p10hat=p10hat,modelSummaries=modelSummaries)
										  return(BayesResults)
										}
										#############################
										#############################
										#############################

# focalTaxlist= Taxlist
	# #TaxonModelSummaries<-model_Bayesian(focalTaxlist)
# TaxonModelSummaries<-readRDS("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/MBON/m2w/OccupancyModeling/COI/12_31_20_ModelSummaries.rds")
	# names(TaxonModelSummaries)<-colnames(taxData)


focalOTUlist= OTUlist #[sample(1:length(OTUlist), 10)]  #select dataset from OTUs for analysis
modelSummaries<-readRDS("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/TideEffectHoodCanal_Gallego/Analysis/OccupancyModeling/23_01_36_ModelSummaries.rds")
#modelSummaries = model_Bayesian(focalOTUlist)  #  run model and store results
	names(modelSummaries)<-row.names(myRealData)

##having estimated the parameters, plug into species occupancy model and see what prob of occurrence is given X detections:
#define probability of occupancy, given parameters, according to binomial distribution and equation in Lahoz-Monfort
ProbOcc=function(x, psi, p11, p10, K){
	(psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
	}

#calculate this probability for each OTU; NOTE -- structure of modelSummaries list changes depending upon whether you generate it as above (as a list) or read it in using readRDS(), as below. 

psimeans=vector(); psi2.5=vector(); psi97.5<-vector()
p11means=vector(); p112.5=vector(); p1197.5<-vector()
p10means=vector(); p102.5=vector(); p1097.5<-vector()

#create an occupancy prob for each OTU
ProbOccOut=NA
for (i in 1:length(modelSummaries)){
	K = 3 #replicates per site 
	psi = modelSummaries[[i]][1] #prob of site occupancy
	p11 = modelSummaries[[i]][2] #prob of detection, given occupancy (i.e., true positive rate)
	p10 = modelSummaries[[i]][3] #prob of detection, given lack of occupancy (i.e., false positive rate)
	psilower=modelSummaries[[i]][9] #2.5% CI
	p11lower=modelSummaries[[i]][10]
	p10lower=modelSummaries[[i]][11]
	psiupper=modelSummaries[[i]][25]
	p11upper=modelSummaries[[i]][26]
	p10upper=modelSummaries[[i]][27]	
	nObs= max(rowSums(focalOTUlist[[i]], na.rm=T)[rowSums(focalOTUlist[[i]], na.rm=T)>0]) #sum(OTUlist[[i]], na.rm=T) #treating all occurrences together #min(rowSums(OTUlist[[i]], na.rm=T)[rowSums(OTUlist[[i]], na.rm=T)>0]) #here, using the min number of detections of that OTU across all transects, all sites... essentially asking if there is any transect for which this OTU is likely a false positive.  There are other ways of doing this.
ProbOccOut[i]=ProbOcc(nObs, psi, p11, p10, K)  

psimeans=c(psimeans, rep(psi, Nsites)) #; psi2.5=vector(); psi97.5<-vector()
p11means=c(p11means, rep(p11, Nsites)) #; p112.5=vector(); p1197.5<-vector()
p10means=c(p10means, rep(p10, Nsites)) #; p102.5=vector(); p1097.5<-vector()
psi2.5=c(psi2.5, rep(psilower, Nsites))
p112.5=c(p112.5, rep(p11lower, Nsites))
p102.5=c(p102.5, rep(p10lower, Nsites))
psi97.5<-c(psi97.5, rep(psiupper, Nsites))
p1197.5<-c(p1197.5, rep(p11upper, Nsites))
p1097.5<-c(p1097.5, rep(p10upper, Nsites))
}

#ProbOccOut #at OTU level; that is, the overall probability of the OTU being real...regardless of individual site of detection.

# # #get the same, but for each taxon, rather than each OTU.  Note this is a way more efficient manner of computing these values than I did above, which was dumb.
# psi_tax<-unlist(TaxonModelSummaries)[seq(1, length(unlist(TaxonModelSummaries)), 36)]; psi_tax<-rep(psi_tax, each=Nsites)
# p11_tax<-unlist(TaxonModelSummaries)[seq(2, length(unlist(TaxonModelSummaries)), 36)]; p11_tax <-rep(p11_tax, each=Nsites)
# p10_tax<-unlist(TaxonModelSummaries)[seq(3, length(unlist(TaxonModelSummaries)), 36)]; p10_tax <-rep(p10_tax, each=Nsites)

# #create df for taxon-specific probabilities
# taxDataOut <-data.frame(taxDataOut, 
	# psi_tax,
	# p11_tax,
	# p10_tax,
	# ProbOcc=ProbOcc(NtaxDetections, psi_tax, p11_tax, p10_tax, rep(3, length(psi_tax)))
	# )

	# saveRDS(taxDataOut, file=paste0("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/MBON/m2w/OccupancyModeling/COI/ProbOcc1COI_", format(Sys.time(), "%H_%M_%S"), "taxon.RDS"))


#below, calc site-specific occupancy probabilities and add to df
dataOUT<-data.frame(dataOUT, psimeans, p11means, p10means, psi2.5, p112.5, p102.5, psi97.5, p1197.5, p1097.5, probOcc<-ProbOcc(dataOUT$Ndetections, psimeans, p11means, p10means, dataOUT$Nrep))
	names(dataOUT)[ncol(dataOUT)]<-"probOcc"

#add taxonomy to dataOUT
# dataOUT$taxon<-tax[match(dataOUT$OTUname, tax[,1]),2]
# dataOUT$TaxHash<-tax[match(dataOUT$OTUname, tax[,1]),3]

dataOUT$locus="COI"

#write data out; much more efficient as Rdata than as csv
#write.csv(dataOUT, "/Users/rpk/GoogleDrive/Kelly_Lab/Projects/MBON/m2w/OccupancyModeling/COI/ProbOccCOI_dataOUT_rarefied.csv")
save(dataOUT, file=paste0("ProbOccCOI_dataOUT", format(Sys.time(), "%H_%M_%S"),".Rdata"))



# # #focalOTUlist[[which.min(ProbOccOut)]]
# ProbOccOut = data.frame(colnames(myRealData), ProbOccOut) ; names(ProbOccOut)=c("OTUname", "ProbOccupancy")
# write.csv(ProbOccOut, paste0(format(Sys.time(), "%H_%M_%S"),"_OTU_probs_wholeDataset.csv"), row.names=F)


# b=ProbOccOut[ProbOccOut[,2]>0.8,1]  #here, applying an 80% likelihood cutoff, which is a bit unsatisfying
# d=row.names(OTUs)%in%b
# write.csv(OTUs[d,], "OTUs_BayesianVetted_wholeDataset.csv", row.names=T, quote=F)

