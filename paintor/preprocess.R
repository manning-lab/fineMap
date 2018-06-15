################ preprocess.R ################ 
# R --vanilla --args ${sep = ":" interval} ${gds_file} ${sep="," sample_ids} ${sep="," assoc_files} ${annotation_file} ${sep="," anno_cols} < /finemapping/paintor/preprocess.R

# Load packages
packages <- c("data.table","SeqArray", "GenomicRanges", "SNPRelate")
lapply(packages, library, character.only = TRUE)

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
interval <- unlist(strsplit(input_args[1], ":"))
gds.file <- input_args[2]
sample.ids.files <- unlist(strsplit(input_args[3],","))
assoc.files <- unlist(strsplit(input_args[4],","))
annotation.file <- input_args[5]
anno.cols <- unlist(strsplit(input_args[6],","))
mac <- as.numeric(inputs_args[7])
pval.col <- input_args[8]
effect.col <- input_args[9]

# parse interval
chr <- interval[1]
start <- as.numeric(as.character(interval[2]))
end <- as.numeric(as.character(interval[3]))

# Load sample ids
sample.ids <- lapply(sample.ids.files, function(x) fread(x, data.table = F, stringsAsFactors = F, header = F)$V1)

# make the output names for zcol and ld
names.suf <- as.character(seq(1,length(sample.ids)))
zcol.names <- paste("ZSCORE", names.suf, sep = ".")
ld.names <- paste("LD", names.suf, sep = ".")

# open gds files
gds.data <- seqOpen(gds.file)

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
}

# get the intersection
var.union = Reduce(intersect, var.ids)
##

## Now calculate the LD matricies
# reset filter (not sure if this is necessary)
seqResetFilter(gds.data)

for (gind in seq(1,length(sample.ids))){
  # calculate LD
  ld <- data.frame(snpgdsLDMat(gds.data, method = "corr", slide = 0, sample.id = sample.ids[[gind]], snp.id = var.union)$LD)
  ld <- ld * ld
  
  # save it
  write.table(ld, file = paste0("Locus1.",ld.names[gind]), row.names = F, col.names = F, sep = " ", quote = F)
}
##

## Now other outputs
# prepare marker output
seqResetFilter(gds.data)
seqSetFilter(gds.data, variant.id = var.union)

# gather the data
chr.v <- seqGetData(gds.data, "chromosome")
pos <- seqGetData(gds.data, "position")
ref <- seqGetData(gds.data, "$ref")
alt <- seqGetData(gds.data, "$alt")

# close gds
seqClose(gds.data)

markers <- data.frame(cbind(chr.v, pos, ref, alt))
names(markers)[1] <- "chr"
markers$marker <- paste(markers$chr, markers$pos, markers$ref, markers$alt, sep=":")

# write out the table
write.table(markers[,c(2,3,4,5)], file = "Locus1.markers.csv", row.names = F, col.names = T, sep = ",", quote = F)
##

## This will process the annotation data
# Load
anno.data <- read.table(annotation.file, sep="\t", header=F, as.is=T)

# make sure chr columns are in same format
if (startsWith(anno.data[1,1], "chr") && !(startsWith(chr, "chr"))){
  anno.data[,1] <- sub("chr","",anno.data[,1])
} else if (!(startsWith(anno.data[1,1], "chr")) && startsWith(chr, "chr")){
  anno.data[,1] <- sub("^","chr",anno.data[,1])
} else {
  pass
}

# make anno into granges object
anno.data.gr <- GRanges(seqnames = anno.data[,1], ranges = IRanges(start = anno.data[,2], end = anno.data[,3]), state = anno.data[,4])

# make markers into granges object
markers.gr <- GRanges(seqnames = markers$chr, ranges = IRanges(start = as.numeric(as.character(markers$pos)), end = as.numeric(as.character(markers$pos))), ref = markers$ref, alt = markers$alt)

# get overlap
markers.ovp <- findOverlaps(markers.gr, anno.data.gr)

# now get the right states for each variant
markers.gr$state <- anno.data.gr[subjectHits(markers.ovp),]$state

# back to data frame
markers <- merge(markers, as.data.frame(markers.gr)[,c("seqnames", "start","ref","alt","state")], by.x = c("chr", "pos","ref","alt"), by.y = c("seqnames", "start","ref","alt"), all.x = T)

# encode states as int positions in |#states| vector
state.map <- data.frame(state = unique(anno.data[,4]))
state.map$state.ind <- seq(1,nrow(state.map))
markers <- merge(markers, state.map, by.x = "state", by.y = "state")

# make matrix of 0s for 1 hot vecs
anno.matrix <- matrix(data = 0, nrow = nrow(markers), ncol = nrow(state.map))

# make the binary vecs
for (i in seq(1,nrow(anno.matrix))){
  anno.matrix[i,markers$state.ind[i]] <- 1
}

# # remove values for annotations that we dont want
# not.inds <- state.map[!(state.map$state %in% anno.cols), "ind"]
# anno.matrix[,not.inds] <- 0

# save annotations
write.table(anno.matrix, file="Locus1.annotations", quote=F, sep=" ", row.names=F, col.names = state.map$state)
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

# final assoc to save
final <- markers[,c("chr","pos",zcol.names)]
names(final)[1] <- "CHR"

# save assoc file
write.table(final, file="Locus1", sep=" ", row.names=F, quote=F)

# export col names for zscores
write.table(zcol.names, file = "zcol.txt", row.names = F, col.names = F, sep = "\n", quote = F)

# export the ld names
write.table(ld.names, file = "ld.txt", row.names = F, col.names = F, sep = "\n", quote = F)

