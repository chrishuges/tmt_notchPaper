# TODO: Add comment
# 
# Author: cshughes
###############################################################################
######################
####required libraries
######################
library(mzR)
library(XML)
library(parallel)
library(RColorBrewer)
######################
####set the output directory for write files
######################
setwd(dir="/Users/cshughes/Documents/projects/tmtPlexing/notchPaper_concentrationTests/Routput/")
######################
####required functions
######################
quant_ms3 = function(pks, massError, ncpu=1) {
	## pks, given the peak information,i.e matrix with mz in first col and intensity in second column
	## return quantification for all the channels as a matrix
	
	###############################################################################
	## quantification with ion scan
	################################################################################    
	## some settings
	reporterTags = c(126.1277261,127.124761,128.1344357,129.1314706,130.1411453,131.1381802)
	names(reporterTags) = c("126","127","128","129","130","131")
	## function to extract the intensity
	extractFun = max    
	toleranceMass = massError ###20*1/1e6  ## 20 ppm 
	
	myrange = cbind(reporterTags* (1-toleranceMass), reporterTags* (1+toleranceMass))
	## make sure the range don't overlap  
	stopifnot(all(diff(as.vector(sapply(1:nrow(myrange), function(i) myrange[i,]))) >0))
	myrangeList = lapply(1:nrow(myrange), function(i) myrange[i,])
	
	pksMaxIntMS3 = do.call(rbind, mclapply(
					pks, function(thisPk){
						## for each MS3 spec extract the maximum in the given range
						sapply(myrangeList, function(x) {                    
									sel = thisPk[,1] >= x[1] & thisPk[,1] <= x[2]
									ifelse(all(!sel), NA, extractFun(thisPk[sel,2]))
								})        
					}, mc.cores=ncpu
			))
	colnames(pksMaxIntMS3) = names(reporterTags)
	pksMaxIntMS3
}
get_inj <- function(xml, name) {
	xpathSApply(xml, sprintf("//x:cvParam[@accession='%s']", name),
			xmlGetAttr, "value",
			namespaces = c(x = "http://psi.hupo.org/ms/mzml"))
}
######################
####examine the raw data files first
######################
#get the list of the input files
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/notchPaper_concentrationTests", pattern=".mzML", full.names=TRUE)
#create an empty list for the output
ms3set = list()
#set the error for TMT calculations
setError = 50*1/1e6
#loop over the infiles
#this function gets a bit hairy if you try to do all the files at once
for (i in 1:6){
	#read the ms file
	ms = openMSfile(infiles[i])
	#get the ms file as an xml file
	xml <- xmlParse(infiles[i])
	#get the header
	hd = cbind(header(ms), ionInjectionTime = as.double(get_inj(xml, "MS:1000927")))
	#find the location of MS3 spectra
	hdMS3 = hd[hd$msLevel ==3,]
	#print number of MS3 spectra
	message(paste('Number of MS3 spectra = ',nrow(hdMS3),sep=''))
	#quant the data
	eSet = as.data.frame(quant_ms3(peaks(ms, hdMS3$seqNum), setError, ncpu=4))
	hdMS3[,23:28] = eSet
	#calculate the raw number of ions
	hdMS3[,23:28] = apply(hdMS3[,23:28],2,function(x) x*hdMS3$ionInjectionTime/1000)
	#remove rows with all expression values as NA
	hdMS3 = subset(hdMS3, rowSums(is.na(hdMS3[,23:28]))<4)
	colnames(hdMS3)[23:28] = c('x126','x127','x128','x129','x130','x131')
	#calculate number of NA values
	hdMS3$missing = rowSums(is.na(hdMS3[,23:28]))
	#output the data
	ms3set[[i]] = hdMS3
	#ms3set[[i]] = hd
	#name the entry
	names(ms3set)[i] = infiles[i]
	#add an output counter
	message(paste('finished ',i,' files',sep=''))
}


#save the data objet
#saveRDS(lsetout, 'ch_OvC_CellLines_CIT-OT_processedIntensities-mzML.rds')

######################
####plot the notch
######################
#lsetout = readRDS('ch_Notch_ConcentrationPanel_Histograms.rds')
#grab some colors
lcols = brewer.pal(6,'Blues')
#make the plot for ion counts
pdf('ch_Notch_ConcentrationPanel_ionCounts_histograms.pdf')
for (i in 1:6){
	hist(log10(ms3set[[i]][,23]),breaks=200, main = names(ms3set)[i], xlab = 'log10(TMT126 Intensity)', lwd=0.05, xlim = c(1,5),
		border = 'gray40',
		col = 'gray95')
	abline(v = mean(log10(ms3set[[i]][,23]),na.rm=TRUE),lty=2,col=lcols[6],lwd=4)
}
dev.off()

lcols = brewer.pal(6,'Blues')
#make the plot for intensity
pdf('ch_Notch_ConcentrationPanel_intensities_histograms.pdf')
for (i in 1:6){
	hist(log10(ms3set[[i]][,23]),breaks=200, main = names(ms3set)[i], xlab = 'log10(TMT126 Intensity)', lwd=0.05, xlim = c(2,7),
			border = 'gray40',
			col = 'gray95')
	abline(v = mean(log10(ms3set[[i]][,23]),na.rm=TRUE),lty=2,col=lcols[6],lwd=4)
}
dev.off()

######################
####get the proportion of scans above the notch
######################
nth126 = hist(log10(ms3set[[1]][,23]),breaks=200,plot=FALSE)
#calculate the proportion above the notch
p126 = sum(nth126$counts[which(nth126$counts == min(nth126$counts[20:40]))[1]:length(nth126$counts)])/sum(nth126$counts)









