library(data.table)

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
all.results <- unlist(strsplit(input_args[1], ","))
posterior.threshold <- as.numeric(input_args[2])

all.intervals <- lapply(all.results, function(x) sub(".results","",basename(x)))

all.list <- list()
for (find in seq(1,length(all.results))){
	res <- fread(all.results[find], data.table = F, stringsAsFactors = F)
	res$interval <- all.intervals[find]
	all.list[[find]] <- res
}

all.results <- do.call(rbind, all.list)

all.results <- all.results[all.results$Posterior_Prob > posterior.threshold,]

fwrite(all.results, file = paste0("all.variants.posterior.", posterior.threshold, 
	".csv"), sep = ",")