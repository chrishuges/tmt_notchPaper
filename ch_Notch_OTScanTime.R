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
setwd(dir="/Users/cshughes/Documents/projects/tmtPlexing/notchPaper_OTscanTime/Routput/")
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
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/notchPaper_OTscanTime", pattern=".mzML", full.names=TRUE)
#create an empty list for the output
ms3set = list()
#set the error for TMT calculations
setError = 50*1/1e6
#loop over the infiles
#this function gets a bit hairy if you try to do all the files at once
for (i in 1:8){
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


###plot the reporter ion data to see the notch
#grab some colors
lcols = brewer.pal(6,'Blues')
#make the plot
pdf('ch_TMTplexing_OTScanTimeTests_MaxFillTimes_IonCount.pdf')
for (i in 1:length(ms3set)){
	x = as.data.frame(ms3set[[i]])
	y = subset(x, ionInjectionTime != max(x$ionInjectionTime))
	z = subset(x, ionInjectionTime == max(x$ionInjectionTime))
	#plot the data
	hist(log10(y[,23]),breaks=200,main=paste(names(ms3set)[i],' non-max injection scans'), xlab = 'log10(TMT126 Intensity)', lwd=0.05, xlim = c(1,5),
			border = 'gray40',
			col = 'gray95')
	abline(v = mean(log10(y[,23]),na.rm=TRUE),lty=2,col=lcols[6],lwd=4)
	text(3,10,nrow(y))
	hist(log10(z[,23]),breaks=200,main=paste(names(ms3set)[i],' max injection scans'), xlab = 'log10(TMT126 Intensity)', lwd=0.05, xlim = c(1,5),
			border = 'gray40',
			col = 'gray75')
	abline(v = mean(log10(z[,23]),na.rm=TRUE),lty=2,col=lcols[6],lwd=4)
	text(3,10,nrow(z))
}
dev.off()


#calculate the proportion above the notch with max injection consideration
for (i in 1:8){
	x = as.data.frame(ms3set[[i]])
	y = subset(x, ionInjectionTime != max(x$ionInjectionTime))
	z = subset(x, ionInjectionTime == max(x$ionInjectionTime))
	message('Number of non-Max scans = ',nrow(y))
	message('Number of max scans = ',nrow(z))
	#make histograms to determine the notch start
	histNoMax = hist(log10(y[,23]),breaks=200,plot=FALSE)
	histMax = hist(log10(z[,23]),breaks=200,plot=FALSE)
	#get the starting position of the notch
	MaxLoc = which(histMax$counts == min(histMax$counts[20:65]))
	NoMaxLoc = which(histNoMax$counts == min(histNoMax$counts[20:65]))
	#calculate the proportion above the notch
	NoMaxProp = sum(histNoMax$counts[NoMaxLoc[NoMaxLoc>20 & NoMaxLoc<65][1]:length(histNoMax$counts)])/sum(histNoMax$counts)
	MaxProp = sum(histMax$counts[MaxLoc[MaxLoc>20 & MaxLoc<65][1]:length(histMax$counts)])/sum(histMax$counts)
	#print out the values
	message(paste('proportion above notch in 126 = ',round(NoMaxProp,2)*100,'%',sep=''))
	message(paste('proportion above notch in 126, max inj = ',round(MaxProp,2)*100,'%',sep=''))
}




#calculate the proportion above the notch across the total data file
ms3scores = data.frame('t126'=integer())
for (i in 1:8){
	x = as.data.frame(ms3set[[i]])
	histMax = hist(log10(x[,23]),breaks=200,plot=FALSE)
	#get the starting position of the notch
	MaxLoc = which(histMax$counts == min(histMax$counts[20:65]))
	#calculate the proportion above the notch
	MaxProp = sum(histMax$counts[MaxLoc[MaxLoc>20 & MaxLoc<65][1]:length(histMax$counts)])/sum(histMax$counts)
	#print out the values
	message(paste('proportion above notch in 126 = ',round(MaxProp,2)*100,'%',sep=''))
	ms3scores[i,1] = round(MaxProp,2)
}

###plot the proportion values
#make the colors
colot1 = brewer.pal(6,"RdBu")[1]
colot2 = brewer.pal(6,"RdBu")[6]
cols1 = c(rep(colot1,4),rep(colot2,4))
cols2 = col2rgb(cols1)
#make the plot
pdf('ch_TMTplexing_OTScanTimeTests_Target-FillTimes_ScoresPlot.pdf')
plot(NA,xlim=c(0.5,4.5),ylim = c(0.80,1))
lines(ms3scores[c(4,1:3),1],col=colot1,lwd=3,type='b',pch=1,cex=1.5)
lines(ms3scores[c(8,5:7),1],col=colot2,lwd=3,type='b',pch=17,cex=1.5,lty=2)
dev.off()



