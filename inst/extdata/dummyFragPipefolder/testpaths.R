# dummyFragPipefolder/testpaths.R
# Path pointg to a dummyFragPipeFolder to test that the package is working


dummyFragPipepath <- "inst/extdata/dummyFragPipeFolder"
dummyhumandatabse <- "inst/extdata/dummyFragPipeFolder/2023-07-12-decoys-reviewed-contam-UP000005640.fas"
dummyexplabels <- c("Treatment", "Vehicle")
dummyconfigpath <- "inst/extdata/dummyFragPipeFolder/1DTPP-TPPR/TPP-TR_config.csv"
dummytpmappic <- "inst/extdata/dummyFragPipeFolder/TP-MAPresults.PNG"

# This should be the last line.
# Note that names are unquoted.
# Using overwrite = T so every time the script runs the
# updated objects are saved, but the default is overwrite = F
usethis::use_data(dummyFragPipepath, dummyhumandatabse, dummyexplabels, dummyconfigpath, dummytpmappic, overwrite = T)


