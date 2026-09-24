# ------------------------------------------------------------------------------
#
# Install the R packages this repository needs.
#
#   make install                        (or: Rscript install/install_packages.R)
#
# - with renv: restores the exact versions in renv.lock (the tested environment)
# - otherwise: installs any missing package listed in DESCRIPTION from CRAN
#
# ------------------------------------------------------------------------------

options(repos = c(CRAN = 'https://cloud.r-project.org'))

use_renv <- !identical(Sys.getenv('CGCP_NO_RENV'), 'true')

if (use_renv) {
  if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv')
  renv::restore(lockfile = 'renv.lock', prompt = FALSE)
} else {
  desc <- read.dcf('DESCRIPTION', fields = 'Imports')[1, 'Imports']
  pkgs <- trimws(strsplit(desc, ',')[[1]])
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) install.packages(missing)
}

cat('R', as.character(getRversion()), '\n')
cat('All packages available:', all(vapply(
  trimws(strsplit(read.dcf('DESCRIPTION', fields = 'Imports')[1, 'Imports'], ',')[[1]]),
  requireNamespace, logical(1), quietly = TRUE)), '\n')
