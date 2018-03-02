
# LOADING LIBRARIES
library(S4Vectors)
library(GenomicRanges)
library(rtracklayer)

# OBJECTS
setClass("gff", representation(name="character", chrom="character", id="character", coords="matrix", strand="character", parent="character"),
	validity=function(object) {
		if (length(object@coords[1,]) != 2) {
			stop("Matrix coords columns number is different to 2")
		}
		else {
			TRUE
		}
	})

setClass("gffList", representation(elementList="list"),
	validity=function(object) {
		 if (length(object@elementList) == 0) {
		 	stop("List is empty!")
		 }
		 else {
			TRUE
		 }
	})


#METHODS

#OBJECTS CONSTRUCTORS

#CREATE GFF-LIST I HAVE A LOT OF THINGS TO FIX REGARDING TO ID, NAME AND PARENT INFORMATION PARSING

setGeneric("create.gffList", function(filename) standardGeneric("create.gffList"))
setMethod("create.gffList", signature("character"),
	function(filename) {
		cat("\nImporting GFF3 file...\n")
		data <- read.delim(filename, header=F, sep="\t", comment.char="#")
		data <- apply(data, 2, as.character)
		data <- apply(data, 2, as.vector)
		dataList <- list()
		ids <- c()
		min <- 0
		max <- length(data[,1])
		pb <- txtProgressBar(min=min, max=max, style=3)
		for (i in 1:length(data[,1])) {
			# EXTRACT chrom, strand, start AND end DATA FROM gff_file
			chrom <- data[i,1]
			strand <- data[i,7]
			start <- data[i,4]
			end <- data[i,5]
			start <- as.numeric(as.character(start))
			end <- as.numeric(as.character(end))
			type <- data[i,3]
			
			# COMMENTS-STRING SPLIT
			# SELECT id AND name
			name <- c()
			id <- c()
			parent <- c()
			identifier <- c()
			for (item in strsplit(strsplit(as.character(data[i,9]), split=";")[[1]], split="=")) {
				if ("Name" %in% item) {
					name <- item[2]
				}
				else if ("Parent" %in% item) {
					parent <- item[2]
				}
				else if ("ID" %in% item) {
					id <- item[2]
				}
			}

			# IF NO "Name" in gff file
			identifier <- ifelse(length(parent) > 0, )
			name <- ifelse(length(name) > 0, name, paste("unknownName_", id, sep=""))
			parent <- ifelse(length(parent) > 0, parent, paste("unknown_", name, sep=""))
	
			# CREATE gff OBJECT OR UPDATE IT
			if (id %in% ids) {
				# UPDATE coords
				if (type %in% types) {
					dataList[[id]][[type]]@coords <- rbind(dataList[[id]][[type]]@coords, c(start, end))
				}
				else {
					types <- c(types, type)
					gff <- new("gff", chrom=chrom, strand=strand, id=id, name=name, parent=parent, coords=t(as.matrix(c(start, end))))
					dataList[[id]][[type]] <- gff
				}
			}
			else {
				# ADD id to ids AND INITIALIZE coords AND gff OBJECT
				types <- c()
				ids <- c(ids, id)
				types <- c(types, type)
				gff <- new("gff", chrom=chrom, strand=strand, id=id, name=name, parent=parent, coords=t(as.matrix(c(start, end))))
				dataList[[id]] <- list()
				dataList[[id]][[type]] <- gff
			}
			min <- min + 1
			setTxtProgressBar(pb, min)
		}
		cat("\n\n")
		gffList <- new("gffList", elementList=dataList)
		gffList
})

##########################################################

setGeneric("import.augustus", function(file) standardGeneric("import.augustus"))
setMethod("import.augustus", signature(file="character"),
	function (file) {
		exonSizes <- c()
		data <- read.delim(file, header=T, sep="\t", row.names=NULL)
		data <- DataFrame(data)
		data@listData$exonStarts <- as.character(data@listData$exonStarts)
		data@listData$exonEnds <- as.character(data@listData$exonEnds)
		for (i in 1:length(data@listData$exonStarts)){
			exonStarts <- as.numeric(strsplit(as.character(data@listData$exonStarts[i]), split=",")[[1]])
	#		data@listData$exonStarts[i] <- paste(c("0", as.character(exonStarts[2:length(exonStarts)])), collapse=",")
			exonEnds <- as.numeric(strsplit(as.character(data@listData$exonEnds[i]), split=",")[[1]])
			diff <- exonEnds - exonStarts
			diff <- as.character(diff)
			diff <- paste(diff, collapse=",")
			exonSizes <- c(exonSizes, diff)
		}
		data@listData$exonSizes <- exonSizes
		data
})

setGeneric("transform.augustus2bed", function(data) standardGeneric("transform.augustus2bed"))
setMethod("transform.augustus2bed", signature("DataFrame"),
	function(data){
		result <- c()
		result <- cbind(as.character(data@listData$chrom),
						as.character(data@listData$txStart),
						as.character(data@listData$txEnd),
						gsub(" ", "_", as.character(data@listData$name2)), # Changed geneName to name2 because xenoRefFlat table has different field name that xenoRefGene: xenoRefFlat -> geneName; xenoRefGene -> name2
						"0", # Score
						as.character(data@listData$strand),
						as.character(data@listData$cdsStart),
						as.character(data@listData$cdsEnd),
						"0", # RGB
						as.character(data@listData$exonCount),
						gsub("$", ",", as.character(data@listData$exonSizes)),
						gsub("$", "", as.character(data@listData$exonStarts)))
		colnames(result) <- c("#chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
	#	colnames(result) <- names(data)
		write.table(result, "aug2bed.out", sep="\t", quote=F, row.names=F, col.names=F)
		result
})

setGeneric("extract.3utr", function(bed) standardGeneric("extract.3utr"))
setMethod("extract.3utr", signature("GenomicRanges"),
	function(bed) {
		gff_result <- c()
		min <- 0
		max <- length(bed)
		pb <- txtProgressBar(min=min, max=max, style=3)
		for(i in 1:length(bed)) {
			chrom <- as(seqnames(bed), "character")[i]
			geneName <- bed[i]@elementMetadata$name
			cdsStart <- start(bed[i]@elementMetadata$thick)+1
			cdsEnd <- end(bed[i]@elementMetadata$thick)
			strand <- as(bed[i]@strand, "character")
			transM <- cbind(start(bed[i]@elementMetadata$blocks)[[1]], end(bed[i]@elementMetadata$blocks)[[1]])
			if (strand == "+") {
				transM <- cbind(ifelse(transM[,1] - cdsEnd <= 0, 0, transM[,1] - cdsEnd), ifelse(transM[,2] - cdsEnd <= 0, 0, transM[,2] - cdsEnd))
				keep <- rowSums(transM) > 0
				transM <- transM[keep,]
				if (class(transM) == "numeric") {
					gff <- c(chrom, ".", "exon", as.character(cdsEnd + transM[1]), as.character(cdsEnd + transM[2]), ".", strand, ".", paste("Parent=", geneName, ";Name=",geneName, "exon_", i, sep=""))
					gff_result <- rbind(gff_result, gff)
				}
				if (class(transM) == "matrix") {
					gff <- c()
					for (j in 1:length(transM[,1])) {
						gff <- rbind(gff, c(chrom, ".", "exon", as.character(cdsEnd + transM[j,1]), as.character(cdsEnd + transM[j,2]), ".", strand, ".", paste("Parent=", geneName, ";Name=exon_", j, sep="")))
					}
					gff_result <- rbind(gff_result, gff)
				}
			}
			else {
				# CODING
				
				transM <- cbind(ifelse(transM[,1] - cdsStart <= 0, transM[,1] - cdsStart, 0), ifelse(transM[,2] - cdsStart <= 0, transM[,2] - cdsStart, 0))
				keep <- rowSums(transM) < 0
				transM <- transM[keep,]
				if (class(transM) == "numeric"){
					gff <- c(chrom, ".", "exon", as.character(cdsStart + transM[1]), as.character(cdsStart + transM[2]), ".", strand, ".", paste("Parent=", geneName, ";Name=exon_", i, sep=""))
					gff_result <- rbind(gff_result, gff)
				}
				if (class(transM) == "matrix") {
					gff <- c()
					for (j in 1:length(transM[,1])) {
						gff <- rbind(gff, c(chrom, ".", "exon", as.character(cdsStart + transM[j,1]), as.character(cdsStart + transM[j,2]), ".", strand, ".", paste("Parent=", geneName, ";Name=exon_", j, sep="")))
					}
					gff_result <- rbind(gff_result, gff)
				}
			}
			min <- min + 1
			setTxtProgressBar(pb, min)
		}
		colnames(gff_result) <- c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
		rownames(gff_result) <- NULL
		gff_result
})
