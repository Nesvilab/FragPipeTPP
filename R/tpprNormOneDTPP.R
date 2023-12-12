
#' Performs TPP-R 1DTPP Normalization in FragPipe outputs
#'
#' @param configurepath: Path to configuration file.
#' @param resultPath: FragPipe output folder.
#'
#' @return void
#' @export
#'
#' @examples tpprNormOneDTPP("C:/FragPipeOutputfolder/TPP-TR_config.xlsx","C:/FragPipeOutputfolder")
tpprNormOneDTPP <- function(configurepath,resultPath){

  #Give data to TPPR
  trData <- TPP::tpptrImport(configTable = configurepath,  fcStr = "rel_fc_",idVar= "Prot_ID")

  #Normalize 1DTPP data
  normResults <- TPP::tpptrNormalize(data=trData)

  #Where to save is set
  savingfolder <- file.path(resultPath, "1DTPP-TPPR")
  setwd(savingfolder)

  #For each experiment extract normalized data
  for (x in normResults$normData){

    # Obtain data
    experimentname <- x@annotation["name"]
    normatrix <- normResults[["normData"]][[experimentname]]@assayData[["exprs"]]

    #Double matrix to data frame to save as csv file
    normdataframe <- as.data.frame(normatrix)
    outputname <- sprintf("TPP-TR_%s.csv",experimentname)
    write.csv(normdataframe, outputname)
  }

}
