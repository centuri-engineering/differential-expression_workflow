# ##########################################################################
# This script is used to launch the compilation of the RMarkdown report
# independently of Rstudio interface
# ##########################################################################

WORKING_DIR = "./"

SCRIPT_DIR = file.path( WORKING_DIR, "07_Report")
OUTPUT_DIR = file.path( WORKING_DIR, "07_Report")

rmarkdown::render( input = file.path( SCRIPT_DIR, "diffexp.Rmd"),
                   output_dir = OUTPUT_DIR,
                   output_file  = "diffexp.html",
                   quiet = FALSE)