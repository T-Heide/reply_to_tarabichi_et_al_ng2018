################################################################################
# FILENAME:    '02-create_rds_dataset.R'
################################################################################

# Options:
simResultDir <- "./results/simulations" # Dir containing the simulation results.
datasetDir   <- "./results/datasets"

# Libs: ########################################################################

library(neutralitytestr)
library(dplyr)

# Functions: ###################################################################

paramsFromFileNames <- function(x) {
  
  exprI  <- "([[:digit:]]*)"
  exprF  <- "([[:digit:].]*)"
  expr   <- sprintf("^simulation-mmr_%s-mbr_%s-seed_%s-cst_%s-[[:print:]]*$",
                    exprI, exprF, exprI, exprI)
  base   <- basename(x)
  match  <- regexec(expr, base)
  params <- do.call(rbind, regmatches(base, match))
  
  # Restructure parameter matrix:
  colnames(params) <- c("file","sc_mutation_rate","sc_deltaS","seed","clst")
  rownames(params) <- gsub("-simulated_sequencing[.]tsv$", "", params[,"file"])
  params           <- params[,-1, drop=FALSE]
  
  storage.mode(params) <- "numeric"
  params               <- data.frame(params)
  
  return(params)
}


parseSimFiles <- function(f, ...) {
  # Extract the params from file names:
  cat("- Extracting parameters from file names.\n")
  params <- paramsFromFileNames(f)

  # Neutrality testing:
  cat("- Loading the data.\n")
  data <- lapply(lapply(lapply(f, read.delim), "[", "VAF"), unlist)

  cat("- Performing neutrality tests.\n")
  rsqs       <- list()
  wh         <- sapply(data, function(x) sum(between(x, 0.12, 0.24)) >= 11)
  test       <- lapply(data[wh], neutralitytestr::neutralitytest, ...)
  rsqs[wh]   <- lapply(test, function(x) unlist(x$rsq["metric"])[1])
  params$rsq <- sapply(rsqs, function(x) { if (is.null(x)){return(NA)}
                                           else {return(x)}})
  params$non_neutral <- params$rsq < 0.98
  return(params)
}


parseInBatches <- function(files, batchSiz=200, ...) {
  Nf       <- length(files)
  Nb       <- ceiling(Nf / batchSiz)
  index    <- head(rep(seq_len(Nb), each=batchSiz), Nf)
  splFiles <- split(files, index)
  
  cat(sprintf("Loading %d simulation result files in %d batch(es):\n\n",Nf,Nb))
  
  res_batches <- lapply(seq_along(splFiles), function(i) {
    cat(sprintf("Batch %d/%d:\n", i, Nb))
    res <- parseSimFiles(splFiles[[i]], ...)
    cat("\n")
    return(res)
  })

  return(do.call(rbind, res_batches))
}


# Main: ########################################################################

# Detect cell count files:
countFileMt    <- "^simulation[[:print:]]*-cell_number[.]tsv$"
countFiles     <- list.files(simResultDir, countFileMt, rec=1, full=1)
baseCountFiles <- basename(countFiles)

# Load cell count data:
cellCountData               <- do.call(rbind, lapply(countFiles, read.delim))
cellCountData$subcloneFrac  <- cellCountData$clone2 / cellCountData$total
cellCountData$simID         <- gsub("-cell_number[.]tsv$", "", baseCountFiles)

cellCountData <- cbind(cellCountData, paramsFromFileNames(countFiles))



# Detect sequencing result files:
resFileMt   <- "^simulation[[:print:]]*-simulated_sequencing[.]tsv$"
resFiles    <- list.files(simResultDir, resFileMt, rec=1, full=1)

# Load result data:
resultData     <- parseInBatches(resFiles)
resultDataExt  <- parseInBatches(resFiles,fmin=0.025, fmax=0.45)



# Save as rds files:
dir.create(datasetDir, showWarnings=FALSE, recursive=TRUE)
saveRDS(resultData, file.path(datasetDir, "1f_model_fits.rds"))
saveRDS(resultDataExt, file.path(datasetDir, "1f_model_fits_ext.rds"))
saveRDS(cellCountData, file.path(datasetDir,"cell_counts.rds"))

