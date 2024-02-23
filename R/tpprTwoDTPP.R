tpprTwoDTPP <- function(configfile_path){
  datalist <- c()
  config_dataframe <- read.csv(configfile_path)

  for (filepath in config_dataframe$Path){

    print(filepath)
    #Get the temperature for each experimental folder title
    folderspl <- unlist(strsplit(filepath, "/"))
    folderlabel <- folderspl[length(folderspl)]
    expname <- gsub(".txt", "", folderlabel)

    datalist[[expname]] <- read.csv(filepath, sep = "\t")


  }



}

#Test
twofragpipe <- "Z:/crojaram/TPP_Project/PXD012423/2DTPP/ATP_rep1/FP20-1_build23"
conc_labels <- c(0,0.005,0.05,0.5,2)
labels_exp <- c("42_44","46_48","50_52","54_56","58_60","62_64")
compound <-("ATP")
configtwo <- tmtitotppr_2D(twofragpipe,labels_exp,conc_labels,compound)
tpprTwoDTPP(configtwo)
