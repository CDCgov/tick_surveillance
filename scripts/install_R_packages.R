#
# This script installs R packages needed for the pipeline
# from downloaded tarballs
#
# Mark Stenglein 10/27/2022
#

#
# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  r_tar_dir=args[1]
  r_lib_dir=args[2]
} else {
  # if running via RStudio
  r_tar_dir = "../lib/R/"
  r_lib_dir = "../lib/R/"
}


# install zip lib - needed for openxlsx
install.packages(paste0(r_tar_dir, "zip_2.2.2.tar.gz"), repos=NULL, lib=r_lib_dir)
# confirm it worked 
library(zip, lib.loc=r_lib_dir)

# install openxlsx
install.packages(paste0(r_tar_dir, "openxlsx_4.2.5.1.tar.gz"), repos=NULL, lib=r_lib_dir)
# confirm it worked 
library(openxlsx, lib.loc=r_lib_dir)
