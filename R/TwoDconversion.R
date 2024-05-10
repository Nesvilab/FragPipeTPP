library(readr)
library(magrittr)
library(dplyr)
library(tidyr)
library(TPP2D)

tmtiheader_to_tpprheaders2D <- function(directory_of_interest, configtemperatures){

  #Find annotation file: Connection between temperature/concentration and TMT label
  annotationfile <- ""
  files_in_specific_directory <- list.files(directory_of_interest)
  for (file in files_in_specific_directory){
    if (grepl("annotation.txt", file)) {
      annotationfile <-  file.path(directory_of_interest, file)

    }
  }
  # Initiate a list to rename column headers
  output_renaming_dict <- list('Protein ID' = 'Prot_ID', "Gene" = "clustername", "Unique Peptides" = "qupm",
                               'Unique Spectral Count' = "qusm")

  # Dictionary to translate TMT labels to TPP-R labels
  tppr_labels_dict <- list(
    "126" = "sumionarea_126", "127N" = "sumionarea_127L", "127C" = "sumionarea_127H",
    "128N" = "sumionarea_128L", "128C" = "sumionarea_128H", "129N" = "sumionarea_129L",
    "129C" = "sumionarea_129H", "130N" = "sumionarea_130L", "130C" = "sumionarea_130H",
    "131N" = "sumionarea_131L"
  )

  #Concentration lables to use fro TPP-R

  #Initiate a list to stores TMT lable with temperature+concentration connection
  tmt_tempconc_dict <- list()
  # Read annotation file
  lines <- readLines(annotationfile)
  for (line in lines) {
    line <- gsub("\n", "", line)
    splits <- strsplit(line, "\t")[[1]]

    # Obtain TMT label and experimental label (temperature_concentration)
    tmt_label <- splits[1]
    temp_conc <- splits[2]

    #TMT label as key, temperature+concentration as value
    tmt_tempconc_dict[[tmt_label]] <- temp_conc

    # Add the experimental label as key and tmt label as value
    output_renaming_dict[[temp_conc]] <- tppr_labels_dict[[tmt_label]]


  }

  #"Dictionary" to output TMT lable and concnetration relations per temperature
  output_dict <- list()

  #Get the temperature for each experimental folder title
  folderspl <- unlist(strsplit(directory_of_interest, "/"))
  folderlabel <- folderspl[length(folderspl)]
  configtemperatures <- unlist(strsplit(folderlabel, "_"))

  #For the temperatures in the experiment
  for (temp_val in configtemperatures[1:2]){
    #Create a sublist for concentration values
    conc_vals_for_tppr <- list()
    #For each TMT label
    for (labeltmt in names(tmt_tempconc_dict)){

      # Change formatting to match the one needed by TPP-R config file
      config_tmt_label <- gsub("sumionarea_", "", tppr_labels_dict[[labeltmt]])
      label_spl <- unlist(strsplit(as.character(tmt_tempconc_dict[labeltmt]), "_"))
      #Get concentration
      config_conc_label <- label_spl[2]
      #Get temperature
      temp_flag <- label_spl[1]
      #If the concentration happens to be at the current temperature being looked at save the concentration
      if (temp_val==temp_flag){
        conc_vals_for_tppr [[config_tmt_label]] <- as.numeric(config_conc_label)
        output_dict[[temp_val]] <- conc_vals_for_tppr
      }else{
        conc_vals_for_tppr [[config_tmt_label]] <- "-"
        output_dict[[temp_val]] <- conc_vals_for_tppr

      }

    }

  }

  return(list(output_renaming_dict, output_dict))
}



fragpipe_to_TPPR <- function(expfolder, configtemperatures) {
  # Function to transform protein.tsv to TPPR input .txt file

  # How to rename column names, assuming renaming_dict is a custom function yet to be defined
  list_oldcoltonewcol_tmt_to_tempdict <- tmtiheader_to_tpprheaders2D(expfolder, configtemperatures)
  oldcoltonewcol <- list_oldcoltonewcol_tmt_to_tempdict[[1]]
  tmt_to_tempdict <- list_oldcoltonewcol_tmt_to_tempdict[[2]]

  # Create data frame by reading the file
  protfile <- "protein.tsv"
  protsv_path <-  file.path(expfolder, protfile)
  print("protsv_path")
  print(protsv_path)
  #This creates a tibble, incompatible with the normalization downstream, bult makes it easier to do other stuff
  DT <- readr::read_tsv(protsv_path)


  # Copy original dataframe to prevent changes to the main DT
  outputDT <- DT

  # What columns to extract from protein.tsv. After "Indistinguishable Proteins" there are the
  # columns with abundance info
  thecolumns <- colnames(outputDT)
  referenceindex <- match("Indistinguishable Proteins", thecolumns)
  protindex <- match("Protein ID", thecolumns)
  geneindex <- match("Gene", thecolumns)
  uniquepepindex <- match("Unique Peptides", thecolumns)
  uniquespectralindex <- match("Unique Spectral Count", thecolumns)
  outputcolumns <- thecolumns[(referenceindex + 1):length(thecolumns)]

  # Add as first column the one that contains the protein ID (the second column in the protein.tsv file), then the unique peptides per protein
  outputcolumns <- c(thecolumns[protindex], thecolumns[geneindex], thecolumns[uniquepepindex], thecolumns[uniquespectralindex], outputcolumns)

  # Data frame with desired columns
  finalDT <- dplyr::select(outputDT, all_of(outputcolumns))

  # Rename to TPP-R compatible columns
  names(finalDT) <- oldcoltonewcol[names(finalDT)]

  #128H and 131L are used as ref columns
  newcolumns <- colnames(finalDT)

  #TPP2D compatability


  finalDT["ref_fc_126"] <- finalDT["sumionarea_126"]/finalDT["sumionarea_128H"]
  finalDT["ref_fc_127L"] <- finalDT["sumionarea_127L"]/finalDT["sumionarea_128H"]
  finalDT["ref_fc_127H"] <- finalDT["sumionarea_127H"]/finalDT["sumionarea_128H"]
  finalDT["ref_fc_128L"] <- finalDT["sumionarea_128L"]/finalDT["sumionarea_128H"]
  finalDT["ref_fc_128H"] <- finalDT["sumionarea_128H"]/finalDT["sumionarea_128H"]
  finalDT["ref_fc_129L"] <- finalDT["sumionarea_129L"]/finalDT["sumionarea_131L"]
  finalDT["ref_fc_129H"] <- finalDT["sumionarea_129H"]/finalDT["sumionarea_131L"]
  finalDT["ref_fc_130L"] <- finalDT["sumionarea_130L"]/finalDT["sumionarea_131L"]
  finalDT["ref_fc_130H"] <- finalDT["sumionarea_130H"]/finalDT["sumionarea_131L"]
  finalDT["ref_fc_131L"] <- finalDT["sumionarea_131L"]/finalDT["sumionarea_131L"]


   #finalDT <- finalDT %>% tibble::add_column(ref_fc_127L = finalDT$sumionarea_127L/finalDT$sumionarea_128H, .before = "sumionarea_126")
  #finalDT <- finalDT %>% tibble::add_column(ref_fc_127H = finalDT$sumionarea_127H/finalDT$sumionarea_128H, .before = "sumionarea_126")
  #finalDT <- finalDT %>% tibble::add_column(ref_fc_128L = finalDT$sumionarea_128L/finalDT$sumionarea_128H, .before = "sumionarea_126")
  #finalDT <- finalDT %>% tibble::add_column(ref_fc_128H = finalDT$sumionarea_128H/finalDT$sumionarea_128H, .before = "sumionarea_126")
  #finalDT <- finalDT %>% tibble::add_column(ref_fc_129L = finalDT$sumionarea_129L/finalDT$sumionarea_131L, .before = "sumionarea_126")
  #finalDT <- finalDT %>% tibble::add_column(ref_fc_129H = finalDT$sumionarea_129H/finalDT$sumionarea_131L, .before = "sumionarea_126")
  #finalDT <- finalDT %>% tibble::add_column(ref_fc_130L = finalDT$sumionarea_130L/finalDT$sumionarea_131L, .before = "sumionarea_126")
  #finalDT <- finalDT %>% tibble::add_column(ref_fc_130H = finalDT$sumionarea_130H/finalDT$sumionarea_131L, .before = "sumionarea_126")
  #finalDT <- finalDT %>% tibble::add_column(ref_fc_131L = finalDT$sumionarea_131L/finalDT$sumionarea_131L, .before = "sumionarea_126")


  print(finalDT)
  return(list(finalDT, tmt_to_tempdict))
}

# This code assumes that there is a function called renaming_dict already defined.
# The renaming_dict function should return a list with two elements: oldcoltonewcol (a named vector or list for renaming columns) and tmt_to_tempdict.

#' twoDConversion: Function to convert protein.tsv filesinto TPPR input for 2DTPP analysis
#'
#' @param fragpipefolder: Path to the folder where FragPipe saves results
#' @param experimentlables: A vector of strings to determine what labels are use for different experiment types
#' @param concentrationlabels: A vector of strings indicatin the temperatures used in each experiment (plex)
#' @param compound_val: A string indicating the compounf of interest (e.g. ATP)
#'
#' @return configsavepath: Path to the newly created configuration file
#' @export
#'
#' @examples twofragpipe <- "Z:/crojaram/TPP_Project/PXD012423/2DTPP/ATP_rep1/FP20-1_build23"
#' conc_labels <- c(0,0.005,0.05,0.5,2)
#' labels_exp <- c("42_44","46_48","50_52","54_56","58_60","62_64")
#' compound <-("ATP")
#' configtwo <- twoDConversion(twofragpipe,labels_exp,conc_labels,compound)

twoDConversion <- function(fragpipefolder, experimentlabels, concentrationlabels, compound_val){


  # List all folders in the specified directory
  folders <- list.dirs(fragpipefolder, recursive = FALSE)

  #print(folders)

  #Create TPP-TPPR folder if it does not exist already
  tpprfolder <- file.path(fragpipefolder, "2DTPP-TPPR")

  if (dir.exists(tpprfolder)) {
    print("Saving results in 2DTPP-TPPR folder")
  } else {
    print("Creating 2DTPP-TPPR folder")
    dir.create(tpprfolder)
  }

  #Vectors to store congif file columns
  Compound <- c()
  Experiment <- c()
  Temperature <- c()
  TMTsix <- c()
  TMTsevenL <- c()
  TMTsevenH <- c()
  TMTeightL <- c()
  TMTeightH <- c()
  TMTnineL <- c()
  TMTnineH <- c()
  TMTtenL <- c()
  TMTtenH <- c()
  TMTelevenL <- c()
  refencecol <- c()
  pathcol <- c()

  #list to create comparison columns
  comparison <- list()

  #Find the experimental names
  for (folder in folders) {


    # Select experiment folders to extract the number for replicates and
    #look for those labels among the columns of the tsv file
    for (label in experimentlabels) {

      lscolumnstoget <- list("Index")
      if (grepl(label, folder)) {

        # Splitting the file path using strsplit
        path_components <- unlist(strsplit(folder, "/"))

        #Adding temperature and replicate for future use
        temperatures_replabel <- path_components[length(path_components)]

        #print(paste("temperatures_replabel",temperatures_replabel))

        #For each temperature set
        conditionlabel <- unlist(strsplit(temperatures_replabel, "_"))

        temp_conce_labels_ls <- list()

        for (item in conditionlabel[1:2]) {

          #Temperature
          Experiment <- append(Experiment, temperatures_replabel)
          #Add compound
          Compound <- append(Compound, compound_val)
          #Create new path to file with data
          outputfilename <- paste(temperatures_replabel, "txt", sep = ".")
          Outputname <- file.path(tpprfolder, outputfilename)
          pathcol <- append(pathcol, Outputname )
          #For each concentration value
          for (conce_val in concentrationlabels){
            #Add the columns to extract from TMTI file
            temp_conc_label <- paste(item,conce_val, sep = "_")
            lscolumnstoget[[length(lscolumnstoget) + 1]] <- paste("X",temp_conc_label, sep="")

          }
        }

        #protein.tsv to TPPR

        convertedfile <- fragpipe_to_TPPR(folder, experimentlabels)


        #Use annotation file to replace the temperature_experiment with appropriate TMT -labels
        #tpprheaders <- tmtiheader_to_tpprheaders2D(folder, Temperature)

        #Extract the TMT-label to Celcius conversion
        #configtempvals <- tpprheaders[2][[1]]
        configtempvals <- convertedfile[2][[1]]

        #For each temperature extract the appropriate TMT lables use for the available concentration
        for (temperature_val in names(configtempvals)){

          Temperature <- append(Temperature, temperature_val)

          TMTsix <- append(TMTsix, configtempvals[[temperature_val]]$`126` )
          TMTsevenL <- append(TMTsevenL, configtempvals[[temperature_val]]$`127L` )
          TMTsevenH <- append(TMTsevenH, configtempvals[[temperature_val]]$`127H` )
          TMTeightL <- append(TMTeightL, configtempvals[[temperature_val]]$`128L` )
          TMTeightH <- append(TMTeightH, configtempvals[[temperature_val]]$`128H` )
          TMTnineL <- append(TMTnineL, configtempvals[[temperature_val]]$`129L` )
          TMTnineH <- append(TMTnineH, configtempvals[[temperature_val]]$`129H` )
          TMTtenL <- append(TMTtenL, configtempvals[[temperature_val]]$`130L` )
          TMTtenH <- append(TMTtenH, configtempvals[[temperature_val]]$`130H` )
          TMTelevenL <- append(TMTelevenL, configtempvals[[temperature_val]]$`131L` )

          if (configtempvals[[temperature_val]]$`131L` == "0"){
            refencecol <- append(refencecol, "131L")
          }else{
            refencecol <- append(refencecol, "128H")
          }

        }

        #print(convertedfile[1])
        #print(class(convertedfile[1]))

        # Save data frame as a tab-delimited text file
        #write.table(newdataframe, file = Outputname, sep = "\t", row.names = FALSE, quote = FALSE) - use base R (slower)
        #data.table::fwrite(convertedfile[1], file = Outputname, sep = "\t")
        write.table(convertedfile[1], file = Outputname, sep = "\t", row.names = FALSE, quote = FALSE)


      }}
  }

  print("Writing Configuration File...")
  #Create configuration file
  configurationdf <- data.frame(Compound, Experiment, Temperature, TMTsix, TMTsevenL, TMTsevenH, TMTeightL, TMTeightH, TMTnineL, TMTnineH, TMTtenL, TMTtenH, TMTelevenL, refencecol, pathcol)
  names(configurationdf) <- c("Compound", "Experiment", "Temperature", "126","127L", "127H", "128L","128H",	"129L", "129H",	"130L",	"130H",	"131L", "RefCol", "Path")
  configsavepath <- file.path(fragpipefolder, "2DTPP-TPPR", "TPP-TR_config.csv")
  write.csv(configurationdf, configsavepath, row.names = FALSE)

  return(configsavepath)
}


