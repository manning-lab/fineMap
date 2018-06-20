library(data.table)

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
all.results <- unlist(strsplit(input_args[1], ","))
posterior.threshold <- as.numeric(input_args[2])

all.list <- lapply(all.results, function(x) fread(x, data.table = F, stringsAsFactors = F))
all.results <- do.call(rbind, all.list)

all.results <- all.results[all.results$Posterior_Probability < posterior.threshold,]

fwrite(all.results, file = paste0("all.variants.posterior.", posterior_threshold, 
	".csv"), sep = ",")