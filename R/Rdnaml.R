# LIBRARIES

library(S4Vectors)
library(Biostrings)
library(GenomicRanges)
library(msa)

#	LOADING AND PREPARING ORTHOGROUPS SEQUENCES
#	CREATING ORTHOGROUPS DATA STRUCTURE
setGeneric("loadFasta", function(path) standardGeneric("loadFasta"))
setMethod("loadFasta", signature("character"),
	function(path) {
		DNAs <- SimpleList()
		lf <- list.files(path=path, pattern=".ancest.fa")
		
		cat("\nLoading Sequences...\n")
		min <- 0
		max <- length(lf)
		pb <- txtProgressBar(min=min, max=max, style=3)
		
		for (file in lf) {
			dnaSet <- readDNAStringSet(paste(path, file, sep="/"), format="fasta")
			dnaSet <- dnaSet[order(names(dnaSet))]
			DNAs@listData[[file]] <- dnaSet
			min <- min + 1
			setTxtProgressBar(pb, min)
		}

		len <- length(DNAs@listData[[1]])
		for (i in 2:length(DNAs)) {
			if (len != length(DNAs@listData[[i]])) {
				stop("Number of genes between sequences dataset is different! It must be of the same length!")
			}
		}
		
		cat("\nPreparing orthogroups...\n")
		min <- 0
		max <- length(DNAs@listData[[1]])
		pb <- txtProgressBar(min=min, max=max, style=3)

		gNames <- names(DNAs@listData[[1]])
		gNames <- gsub("\\.[LS]\\|.*", "", gNames)
		orthoDNAs <- SimpleList()
		for (i in 1:length(DNAs@listData[[1]])) {
			seqList <- c()
			names <- c()
			for (file in lf) {
				seqList <- c(seqList, toString(DNAs@listData[[file]][i]))
				names <- c(names, names(DNAs@listData[[file]])[i])
			}
			tmp_names <- sapply(lapply(strsplit(names, ""), rev), paste, collapse="")
			tmp_names <- gsub("^.*\\|(.{2,10}).*", "\\1", tmp_names)
			tmp_names <- sapply(lapply(strsplit(tmp_names, ""), rev), paste, collapse="")
			names(seqList) <- tmp_names
			orthoDNAs@listData[[gNames[i]]] <- DNAStringSet(seqList, use.names=TRUE)
			setTxtProgressBar(pb, i)
		}
		cat("\n\n")
		orthoDNAs
})

#	PERFORMING MSA
setGeneric("MSA", function(orthogroups) standardGeneric("MSA"))
setMethod("MSA", signature("SimpleList"),
	function(orthogroups) {
		gNames <- names(orthogroups)
		msa <- SimpleList()
		for (i in 1:length(orthogroups@listData)) {
			msa@listData[[gNames[i]]] <- msa(orthogroups@listData[[i]], substitutionMatrix="default")
			cat(paste(i, " of ", length(orthogroups@listData)))
			cat("\n")
		}
		msa
})

#	TESTING CUTOFFS
setGeneric("testing.cutoff", function(orthogroups) standardGeneric("testing.cutoff"))
setMethod("testing.cutoff", signature("SimpleList"),
	function(orthogroups) {
		cutoffs <- seq(5, 100, by=5)
		accepted_hits_by_cutoff <- SimpleList()
		for (cutoff in cutoffs) {
			accepted <- c()
			for (i in 1:50) {
				musa <- MSA(orthogroups=sample(orthogroups, replace=F, size=10))
				len.accepted <- length(filter.multSA(musa=musa, cutoff=cutoff)$accepted)
				accepted <- c(accepted, (len.accepted / 10) * 100)
			}
			accepted_hits_by_cutoff[[as.character(cutoff)]] <- accepted
		}
		accepted_hits_by_cutoff
})

#	FILTERING multSA object
setGeneric("filter.multSA", function(musa, cutoff) standardGeneric("filter.multSA"))
setMethod("filter.multSA", signature("SimpleList", "numeric"),
	function(musa, cutoff) {
		accepted_musa <- SimpleList()
		rejected_musa <- SimpleList()
		gNames <- names(musa)
		for (i in 1:length(musa@listData)) {
			percentages <- c()
			for (j in 1:length(musa@listData[[i]]@unmasked)) {
				orthoGene <- musa@listData[[i]]@unmasked[[j]]
				t <- table(strsplit(toString(orthoGene), split=""))
				percentages <- c(percentages, ((sum(t[names(t) != "-"])/sum(t)) * 100))
			}
			if (sum(percentages >= cutoff) == length(percentages)) {
				accepted_musa[[gNames[i]]] <- musa@listData[[i]]
			}
			else {
				rejected_musa[[gNames[i]]] <- musa@listData[[i]]
			}
		}
		result <- list(accepted=accepted_musa, rejected=rejected_musa)
		result
})

#	CALLING dnaml OF phylip
setGeneric("perform.dnaml", function(musa) standardGeneric("perform.dnaml"))
setMethod("perform.dnaml", signature("SimpleList"),
	function(musa) {
		
})

