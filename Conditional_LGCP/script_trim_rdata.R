#!/usr/bin/env Rscript
# Usage:
#   nohup Rscript script_trim_rdata.R <folder> > output.log 2>&1 &

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript trim_rdata.R <folder>")

folder <- args[1]
files  <- list.files(folder, pattern = "\\.RData$", full.names = TRUE)
cat(sprintf("[INFO] Found %d .RData files in %s\n", length(files), folder))

for (f in files) {
  cat(sprintf("[LOAD] %s\n", f))
  env <- new.env()
  load(f, envir = env)
  obj_name <- ls(env)[1]
  obj <- get(obj_name, envir = env)
  obj[c("step_2b", "step_3")] <- NULL
  assign(obj_name, obj)
  save(list = obj_name, envir = environment(), file = f)
  cat(sprintf("[DONE] %s\n", f))
}