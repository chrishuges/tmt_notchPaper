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
setwd(dir="/Users/cshughes/Documents/projects/tmtPlexing/notchPaper_dilutionSeries/Routput/")
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
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/notchPaper_dilutionSeries", pattern=".mzML", full.names=TRUE)
#create an empty list for the output
ms3set = list()
#set the error for TMT calculations
setError = 50*1/1e6
#loop over the infiles
#this function gets a bit hairy if you try to do all the files at once
for (i in 1:9){
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





######need to get the scan locations of the spiked peptides
#set working directory
setwd(dir="/Users/cshughes/Documents/projects/tmtPlexing/notchPaper_dilutionSeries/Routput/")
#grab the protein files
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/notchPaper_dilutionSeries", pattern="Proteins\\.txt", full.names=TRUE)
#read in the files into a list
proSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(proSet) = c(sub(".*?1\\_(.*?)(\\_Proteins\\.txt.*|$)", "\\1", infiles))
#do the same for the PSM files
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/notchPaper_dilutionSeries", pattern="PSMs\\.txt", full.names=TRUE)
#read in the files into a list
psmSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(psmSet) = c(sub(".*?1\\_(.*?)(\\_PSMs\\.txt.*|$)", "\\1", infiles))
#do the same for the MS2 files
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/notchPaper_dilutionSeries", pattern="Info\\.txt", full.names=TRUE)
#read in the files into a list
msmsSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(msmsSet) = c(sub(".*?1\\_(.*?)(\\_MSMSSpectrumInfo\\.txt.*|$)", "\\1", infiles))
#and finally for the quan spectra
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/notchPaper_dilutionSeries", pattern="QuanSpectra\\.txt", full.names=TRUE)
#read in the files into a list
quanSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(quanSet) = c(sub(".*?1\\_(.*?)(\\_QuanSpectra\\.txt.*|$)", "\\1", infiles))


######################
####plot the spiked peptides as they approach the notch region
######################
#process the data files to get info of interest out
#make a data holder
tmtData = data.frame()
lsetout = data.frame()
#loop over the files
for (i in 4:6){
	#decide which columns you want
	pepCols = c('Annotated.Sequence','First.Scan','Master.Protein.Accessions','X126','X127','X128','X129','X130','X131')
	pep = as.data.frame(psmSet[[i]])[,pepCols]
	#subset out only spike peptides
	pepSpikes = pep[grepl('SPIKE',pep$Master.Protein.Accessions),]
	#output into tmtData
	tmtData = pepSpikes
	colnames(tmtData)[2] = 'precursorScanNum'
	#need to change these numbers based on what data you want to look at
	x = as.data.frame(ms3set[[i]])
	#keep only scans that correspond with spike peptides
	lsetmerge = merge(tmtData,x,by='precursorScanNum')
	#filter out any that dont give a value in the first channel
	lsetfilt = subset(lsetmerge, !is.na(lsetmerge$x126))
	lsetout = rbind(lsetfilt,lsetout)
}
#make some colors
xCol = col2rgb(brewer.pal(9,'Blues')[7])
#make the plot
pdf('ch_Notch_SpikeDilution_SpikeSignalDistribution_1-30.pdf')
plot(lsetout$acquisitionNum,
		log10(lsetout$x126),
		xlab = 'log10(TMT126 Intensity)',
		col = rgb(xCol[1,],xCol[2,],xCol[3,],75,maxColorValue=255),
		pch = 20,
		cex = 1.25,
		ylim = c(1,5),
		xlim = c(0,80000))
abline(h=2.3,lty=2,lwd=2)
abline(h=3.0,lty=2,lwd=2)
abline(h=1.6,lty=2,lwd=2)
dev.off()




######################
####plot the total set of peptides to visualize the notch region
######################
lsetout = as.data.frame(ms3set[[9]])
lsetfilt = subset(lsetout, !is.na(lsetout$x126))
lsetsamp = lsetfilt[sample(nrow(lsetfilt), 5000), ]
#make some colors
xCol = col2rgb(brewer.pal(9,'Purples')[7])
#make the plot
pdf('ch_Notch_SpikeDilution_TotalSignalDistribution_1-30.pdf')
plot(lsetsamp$acquisitionNum,
		log10(lsetsamp$x126),
		xlab = 'log10(TMT126 Intensity)',
		col = rgb(xCol[1,],xCol[2,],xCol[3,],75,maxColorValue=255),
		pch = 20,
		cex = 1.25,
		ylim = c(1,5),
		xlim = c(0,80000))
abline(h=2.3,lty=2,lwd=2)
dev.off()



######################
####plot the ratios in the 1/30 sample above and below the notch region
######################
#process the data files to get info of interest out
#make a data holder
tmtData = data.frame()
lsetout = data.frame()
#loop over the files
for (i in 1:6){
	#decide which columns you want
	pepCols = c('Annotated.Sequence','First.Scan','Master.Protein.Accessions','X126','X127','X128','X129','X130','X131')
	pep = as.data.frame(psmSet[[i]])[,pepCols]
	#subset out only spike peptides
	pepSpikes = pep[grepl('SPIKE',pep$Master.Protein.Accessions),]
	#output into tmtData
	tmtData = pepSpikes
	colnames(tmtData)[2] = 'precursorScanNum'
	#need to change these numbers based on what data you want to look at
	x = as.data.frame(ms3set[[i]])
	#keep only scans that correspond with spike peptides
	lsetmerge = merge(tmtData,x,by='precursorScanNum')
	#filter out any that are missing any values
	lsetfilt = subset(lsetmerge, rowSums(is.na(lsetmerge[,31:36]))<1)
	lsetout = rbind(lsetfilt,lsetout)
}

#calculate a total intensity
lsetout$totInt = rowSums(lsetout[,31:36],na.rm=TRUE)
lsetout$totSN = rowSums(lsetout[,4:9],na.rm=TRUE)
#get the rows that are below the notch at log10 2.3
up = subset(lsetout, log10(lsetout$x126)>2.3 & log10(lsetout$x126)<3.0)
dn = subset(lsetout, log10(lsetout$x126)<2.3 & log10(lsetout$x126)>1.6)



dn[,31:36] = apply(dn[,31:36],2,function(x) (x/dn$totInt)*24)
up[,31:36] = apply(up[,31:36],2,function(x) (x/up$totInt)*24)
dn[,4:9] = apply(dn[,4:9],2,function(x) (x/dn$totSN)*24)
up[,4:9] = apply(up[,4:9],2,function(x) (x/up$totSN)*24)

dn = dn[order(dn$totInt),]
up = up[order(up$totInt),]







