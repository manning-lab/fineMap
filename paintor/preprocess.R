################ preprocess.R ################ 
# R --vanilla --args ${sep = ":" interval} ${gds_file} ${sep="," sample_ids} ${sep="," assoc_files} ${annotation_file} ${sep="," anno_cols} < /finemapping/paintor/preprocess.R

# Load packages
packages <- c("data.table","SeqArray", "GenomicRanges", "SNPRelate", "dplyr", "tidyr", "SeqVarTools")
lapply(packages, library, character.only = TRUE)

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
print(input_args)

interval <- unlist(strsplit(input_args[1], ":"))
gds.file <- input_args[2]
sample.ids.files <- unlist(strsplit(input_args[3],","))
assoc.files <- unlist(strsplit(input_args[4],","))
meta.file <- input_args[5]
annotation.file <- input_args[6]
anno.cols <- unlist(strsplit(input_args[7],","))
mac <- as.numeric(input_args[8])
pval.col <- input_args[9]
effect.col <- input_args[10]
out.pref <- input_args[11]

## # these are from the DCC pipeline, credit -> S. Gogarten 
.variantDF <- function(gds) {
  data.frame(variant.id=seqGetData(gds, "variant.id"),
             chromosome=seqGetData(gds, "chromosome"),
             position=seqGetData(gds, "position"),
             ref=refChar(gds),
             alt=altChar(gds),
             nAlleles=seqNumAllele(gds),
             stringsAsFactors=FALSE)
}
.expandAlleles <- function(gds) {
  .variantDF(gds) %>%
    separate_rows_("alt", sep=",") %>%
    rename_(allele="alt") %>%
    group_by_("variant.id") %>%
    mutate_(allele.index=~1:n()) %>%
    as.data.frame()
}
##

# parse interval
chr <- interval[1]
start <- as.numeric(as.character(interval[2]))
end <- as.numeric(as.character(interval[3]))
print(chr)
print(start)
print(end)

# create output prefix
out.pref <- paste(out.pref, chr, start, end, sep = ".")
print(out.pref)

# Load sample ids
sample.ids <- lapply(sample.ids.files, function(x) fread(x, data.table = F, stringsAsFactors = F, header = F)$V1)
print(names(sample.ids))

# make the output names for zcol and ld
names.suf <- as.character(seq(1,length(sample.ids)))

zcol.names <- paste("ZSCORE", names.suf, sep = ".")
ld.names <- paste("LD", names.suf, sep = ".")

# open gds files
gds.data <- seqOpen(gds.file)

# make sure that interval and gds have the same chr format
gds.chr <- unique(seqGetData(gds.data,"chromosome"))[1]

if (startsWith(gds.chr,"chr") && !(startsWith(chr, "chr"))){
  chr <- sub("^","chr",chr)
} else if (!(startsWith(gds.chr,"chr")) && startsWith(chr, "chr")){
  chr <- sub("chr","",chr)
}

## This is to get the list of common variants among all groups ##
# list for variant ids
var.ids <- list()

# Loop through sample groups getting list of var
for (gind in seq(1,length(sample.ids))){
  
  # filter to samples
  seqSetFilter(gds.data, sample.id = sample.ids[[gind]])
  
  # filter to interval
  seqSetFilterChrom(gds.data, chr, from.bp = start , to.bp = end)
  
  # filter to mac
  seqSetFilterCond(gds.data, mac = mac)
  
  # store variant ids
  var.ids[[gind]] <- seqGetData(gds.data, "variant.id")
  
  # reset filter
  seqResetFilter(gds.data)
}

# get the intersection
var.union = Reduce(intersect, var.ids)
##

## Now other outputs
# prepare marker output
seqResetFilter(gds.data)
seqSetFilter(gds.data, variant.id = var.union)

# gather the data
markers <- .expandAlleles(gds.data)[,c(1,2,3,4,5)]
names(markers) <- c("variant.id","chr","pos","ref","alt")
markers$marker <- paste(markers$chr, markers$pos, markers$ref, markers$alt, sep=":")
##

## Now get the zscores
# calculate z
for (gind in seq(1,length(assoc.files))){
  c.assoc <- fread(assoc.files[gind], data.table = F, stringsAsFactors = F)[,c("chr","pos","ref","alt",pval.col,effect.col)]
  c.assoc$marker <- paste(c.assoc$chr, c.assoc$pos, c.assoc$ref, c.assoc$alt, sep = ":")
  c.assoc <- c.assoc[c.assoc$marker %in% markers$marker, ]
  c.assoc$Z <- abs(qnorm(c.assoc[, pval.col]/2))*sign(c.assoc[,effect.col])
  c.assoc <- c.assoc[,c("marker","Z")]
  colnames(c.assoc) <- c("marker", zcol.names[gind])
  markers <- merge(markers, c.assoc, by.x = "marker", by.y = "marker", all.x = T)
}

# remove any rows with NAs
markers <- na.omit(markers)

# add meta p-value to final output
# load meta results, must have columns named marker (chr:pos:ref:alt) and Zscore
if (meta.file != "NA"){
  meta <- fread(meta.file, data.table = F, stringsAsFactors = F)#[,c("marker","Zscore")] 
  
  # get the right colname
  if ("marker" %in% names(meta)){
    meta.marker <- "marker"
  } else {
    # this is for metal output
    meta.marker <- "MarkerName"
  }
  
  # only take the cols that we need
  meta <- meta[,c(meta.marker,"Zscore")]
  
  # get the right sep character
  for (sepchar in c("-","\\.",";")){
    if (length(unlist(regmatches(meta[1,meta.marker], gregexpr(sepchar, meta[1,meta.marker])))) == 3){
      meta[,meta.marker] <- gsub(sepchar, ":", meta[,meta.marker])
      break
    }
  }
  
  # make sure chr is in the same format
  if (startsWith(meta[1,meta.marker],"chr") && !(startsWith(chr, "chr"))){
    meta[,meta.marker] <- gsub("chr", "", meta[,meta.marker])
  } else if (!(startsWith(meta[1,meta.marker],"chr")) && startsWith(chr, "chr")){
    meta[,meta.marker] <- gsub("^", "chr", meta[,meta.marker])
  }
  
  # get the colname right
  names(meta)[2] <- "ZSCORE.meta"
  
  # merge with markers
  markers <- merge(markers, meta, by.x = "marker", by.y = meta.marker, all.x = T)
  
  # remove any nas
  markers <- na.omit(markers)
  
  # order based on position
  markers <- markers[order(markers$pos),]
  
  # add meta col to zcol names
  write.table("ZSCORE.meta", file = "meta_zcol.txt", row.names = F, col.names = F, sep = "\n", quote = F)
  # zcol.names <- c(zcol.names, "ZSCORE.meta")
  
  # add meta col to ld names
  ld.names <- c(ld.names, "LD.all")
  
} else {
  markers <- markers[order(markers$pos),]
}
##

## This will process the annotation data
# first check if we want to use annotations
if (anno.cols == "NA" || annotation.file == "NA"){
  anno.cols <- c("baseline")
  anno.matrix <- matrix(data = 0, nrow = nrow(markers), ncol = 1)
  markers$state <- "baseline"
} else {
  # Load
  anno.data <- read.table(annotation.file, sep="\t", header=F, as.is=T)
  
  # make sure chr columns are in same format
  if (startsWith(anno.data[1,1], "chr") && !(startsWith(chr, "chr"))){
    anno.data[,1] <- sub("chr","",anno.data[,1])
  } else if (!(startsWith(anno.data[1,1], "chr")) && startsWith(chr, "chr")){
    anno.data[,1] <- sub("^","chr",anno.data[,1])
  } else {}
  
  # make anno into granges object
  anno.data.gr <- GRanges(seqnames = anno.data[,1], ranges = IRanges(start = anno.data[,2], end = anno.data[,3]), state = anno.data[,4])
  
  # make markers into granges object
  markers.gr <- GRanges(seqnames = markers$chr, ranges = IRanges(start = as.numeric(as.character(markers$pos)), end = as.numeric(as.character(markers$pos))), ref = markers$ref, alt = markers$alt)
  
  # get overlap
  markers.ovp.ids <- findOverlaps(markers.gr, anno.data.gr)
  
  # now get the right states for each variant
  markers.ovp <- markers.gr[queryHits(markers.ovp.ids),]
  markers.ovp$state <- anno.data.gr[subjectHits(markers.ovp.ids),]$state
  
  # back to data frame
  markers <- merge(markers, as.data.frame(markers.ovp)[,c("seqnames", "start","ref","alt","state")], by.x = c("chr", "pos","ref","alt"), by.y = c("seqnames", "start","ref","alt"), all.x = T)
  
  # encode states as int positions in |#states| vector
  state.map <- data.frame(state = unique(anno.data[,4]))
  state.map$state.ind <- seq(1,nrow(state.map))
  markers <- merge(markers, state.map, by.x = "state", by.y = "state")
  markers <- markers[order(markers$pos),]
  
  # make matrix of 0s for 1 hot vecs
  anno.matrix <- matrix(data = 0, nrow = nrow(markers), ncol = nrow(state.map))
  
  # make the binary vecs
  for (i in seq(1,nrow(anno.matrix))){
    anno.matrix[i,markers$state.ind[i]] <- 1
  }
}
##

##
# get final list of variants+scores
if (meta.file != "NA"){
  final <- markers[,c("marker","chr","pos",zcol.names,"ZSCORE.meta","state")]
} else {
  final <- markers[,c("marker","chr","pos",zcol.names,"state")]
}


# rename to fit paintor input
names(final)[2] <- "CHR"
##

## Now calculate the LD matricies
# reset filter (not sure if this is necessary)
seqResetFilter(gds.data)

# outmessage for warning if NA in ld
outmessage <- "No NA's in LD matrix detected."

# place to hold flags for na
na.flags <- rep(0,length(assoc.files)+1)

for (gind in seq(1,length(sample.ids))){
  # calculate LD
  ld <- data.frame(snpgdsLDMat(gds.data, method = "corr", slide = 0, sample.id = sample.ids[[gind]], snp.id = markers$variant.id)$LD)
  ld <- ld * ld
  
  # output message if we have some NAs
  if (any(is.na(ld))){
    na.flags[gind] <- 1
    ld[is.na(ld)] <- 0
  }
  
  # save it
  write.table(ld, file = paste0(out.pref,".",ld.names[gind]), row.names = F, col.names = F, sep = " ", quote = F)
}

# make general LD matrix
all.sample.ids <- do.call(c, sample.ids)
ld <- data.frame(snpgdsLDMat(gds.data, method = "corr", slide = 0, sample.id = all.sample.ids, snp.id = markers$variant.id)$LD)
ld <- ld * ld
if (any(is.na(ld))){
  na.flags[length(na.flags)] <- 1
  ld[is.na(ld)] <- 0
}

if (sum(na.flags) > 0){
  inds <- which(na.flags == 1)
  outmessage <- paste("NA's detected in LD matrices: ", paste(inds,collapse = ", "),". Replacing NA values with 0", sep = "")
}

# save it
write.table(ld, file = paste0(out.pref,".LD.all"), row.names = F, col.names = F, sep = " ", quote = F)
##

# write out the markers
write.table(markers[,c("pos","ref","alt","marker")], file = paste0(out.pref,".markers.csv"), row.names = F, col.names = T, sep = ",", quote = F)

# save annotations
if (anno.cols[1] == "baseline"){
  write.table(anno.matrix, file=paste0(out.pref,".annotations"), quote=F, sep=" ", row.names=F, col.names = c("baseline"))
} else {
  write.table(anno.matrix, file=paste0(out.pref,".annotations"), quote=F, sep=" ", row.names=F, col.names = state.map$state)  
}


# save assoc file
write.table(final, file=out.pref, sep=" ", row.names=F, quote=F)

# export col names for zscores
write.table(zcol.names, file = "zcol.txt", row.names = F, col.names = F, sep = "\n", quote = F)

# export the ld names
write.table(ld.names, file = "ld.txt", row.names = F, col.names = F, sep = "\n", quote = F)

# export the annotation names
write.table(anno.cols, file = "anno.txt", row.names = F, col.names = F, sep = "\n", quote = F)

# export output message
write.table(outmessage, file="out_message.txt", sep=" ", row.names=F,col.names=F, quote=F)

# close gds
seqClose(gds.data)

