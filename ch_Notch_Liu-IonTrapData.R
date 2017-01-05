# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(mzR)
library(XML)
library(parallel)
library(RColorBrewer)
#set working directory
setwd(dir="/Users/cshughes/Documents/projects/tmtPlexing/publishedData/Routput/")
####required functions ###need to change the quant_ms3 for Liu data because they actually use the 10plex reagents
quant_ms3 = function(pks, massError, ncpu=1) {
	## pks, given the peak information,i.e matrix with mz in first col and intensity in second column
	## return quantification for all the channels as a matrix
	
	###############################################################################
	## quantification with ion scan
	################################################################################    
	## some settings
	reporterTags = c(126.1277261,127.131079,128.1344357,129.137787,130.1411453,131.1381802)
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
###process the files to get the reporter ion data out
#get the list of the input files
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/publishedData/Liu", pattern="\\.mzML", full.names=TRUE)
#create an empty list for the output
ms3set = list()
#set the error for TMT calculations, (50*1/1e6 for OT, 500*1/1e6 for IT)
setError = 500*1/1e6
#loop over the infiles
#this function gets a bit hairy if you try to do all the files at once
for (i in 2:6){
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
	#hdMS3[,23:28] = apply(hdMS3[,23:28],2,function(x) x*hdMS3$ionInjectionTime/1000)
	#remove rows with all expression values as NA
	#hdMS3 = subset(hdMS3, rowSums(is.na(hdMS3[,25:30]))<4)
	colnames(hdMS3)[23:28] = c('x126','x127','x128','x129','x130','x131')
	#calculate number of NA values
	hdMS3$missing = rowSums(is.na(hdMS3[,23:28]))
	#output the data
	ms3set[[i]] = hdMS3
	#name the entry
	names(ms3set)[i] = infiles[i]
	#print out some interesting values
	message('Total number of non-NA MS3 spectra = ',length(which(!is.na(ms3set[[i]][,23]))))
	message('Total number of non-NA TMT126 spectra = ',length(which(!is.na(ms3set[[i]][,23]))))
	message('Total number of non-NA TMT127 spectra = ',length(which(!is.na(ms3set[[i]][,24]))))
	message('Total number of non-NA TMT128 spectra = ',length(which(!is.na(ms3set[[i]][,25]))))
	message('Total number of non-NA TMT129 spectra = ',length(which(!is.na(ms3set[[i]][,26]))))
	message('Total number of non-NA TMT130 spectra = ',length(which(!is.na(ms3set[[i]][,27]))))
	message('Total number of non-NA TMT131 spectra = ',length(which(!is.na(ms3set[[i]][,28]))))
	#add an output counter
	message(paste('finished ',i,' files',sep=''))
	message('')
}



#bind all of the data together into a single frame
lsetout = data.frame()
for (x in 2:6){
	lset = as.data.frame(ms3set[[x]])
	lsetout = rbind(lsetout,lset)
}

###plot the reporter ion data to see the notch
#grab some colors
lcols = brewer.pal(6,'Blues')
#make the plot
pdf('ch_Notch_Liu_IonTrap_Intensity_histogram.pdf')
hist(log10(lsetout[,25]),breaks=200, xlab = 'log10(TMT128 Intensity)', lwd=0.05,
			border = 'gray40',
			col = 'gray55')
	abline(v = median(log10(lsetout[,25]),na.rm=TRUE),lty=2,col=lcols[6],lwd=4)
dev.off()

