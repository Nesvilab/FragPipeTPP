
#' Extract protein description for TP-MAP compatability
#'
#' @param filepathdb: Path to protein database (FASTA format)
#'
#' @return list with protein uniprot accession numbers and corresponding description line
#' @export
#'
#' @examples fastaparser("C:/protein.fas")
fastaparser <- function(filepathdb) {
  # Initialize a list
  identifiers <- list()

  # Read the FASTA file
  con <- file(filepathdb, "r")
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    # Check for protein headers
    if (substr(line, 1, 1) == ">") {
      # Split the line using "|" as a delimiter
      parts <- unlist(strsplit(line, "\\|"))
      # Extract the protein ID string
      keyproteinid <- parts[2]
      # Add protein ID as key and whole description line as value
      identifiers[[keyproteinid]] <- line
    }
  }
  close(con)

  return(identifiers)
}



#' TPP-R output to TP-MAP input (for 1DTPP analysis)
#'
#' @param configurationfile_path: Path to configuration file.
#' @param databasepath: Path to protein database (FASTA format)
#'
#' @return void
#' @export
#'
#' @examples tpprTotpmap("C:/FragPipeOutputfolder/TPP-TR_config.xlsx","C:/protein.fas")
tpprTotpmap <- function(configurationfile_path, databasepath){

  #Read configuration file
  config_dataframe <- data.table::fread(configurationfile_path)


  #obtain the path to the config file folder
  tpprfolder <- dirname(configurationfile_path)

  framestojoin <- list()

  #Access the experiments in the config file
  for (exp in config_dataframe$"Experiment"){

    #vector for renaming header
    newnames <- c("Accession")

    #Where in the column
    rownum <- which(config_dataframe$"Experiment" == exp, arr.ind=TRUE)


    #Find each experimental TPP-R result file
    experiment_tpprdf_path <- file.path(tpprfolder,  sprintf("TPP-TR_%s.csv",exp))
    experiment_tpprdf <- data.table::fread(experiment_tpprdf_path)

    #Converting NA to zeroes
    experiment_tpprdf[is.na(experiment_tpprdf)] <- 0


    #Obtain column names from experimental TPP-R result file
    tpprheaders <- colnames(experiment_tpprdf)
    #print(tpprheaders)
    for (label in tpprheaders){
      if (grepl("rel_fc", label)){
        splitlabel <- unlist(strsplit(label, "_"))
        tmtlabel <- splitlabel[3]

        #Need to add the .. before tmtlabel otherwise it thinks I am literally looking for tmtlabel column
        temperaturelabel <- config_dataframe[rownum, tmtlabel]
        mod_exp <- gsub("_", "", exp)
        completetemplabel <- sprintf("Ref_%.1f_%s", temperaturelabel , mod_exp)
        newnames <- append(newnames, completetemplabel)

      }
    }

    #Appropriate TP-MAP header
    names(experiment_tpprdf) <- newnames
    #print(colnames(experiment_tpprdf))
    framestojoin[[length(framestojoin) + 1]] <- experiment_tpprdf
  }
  #print(framestojoin)

  #Merge dataframes
  masterframe <- data.frame()
  for (dataframe in framestojoin){
    if (nrow(masterframe) != 0 && ncol(masterframe) != 0) {
      #print("The data frame is not empty")
      masterframe <- tidyft::inner_join(masterframe, dataframe, by = "Accession")

    } else {
      #print("The data frame is empty")
      masterframe <- dataframe
    }

  }

  #Creating description columnfrom database file
  descriptions_ls <- fastaparser(databasepath)

  descriptions <- c()

  for (protein in masterframe$"Accession"){
    #print(protein)
    #print(descriptions_ls[[protein]])

    descriptions <- append(descriptions, descriptions_ls[[protein]] )
  }

  #print(length(descriptions))
  #print(length(masterframe$"Accession"))

  #Added description column
  masterframe$Description <- descriptions


  #Creating path to save results
  Outputname <- file.path(tpprfolder, "1DTPP.tpmap.txt")

  data.table::fwrite(masterframe, file = Outputname, sep = "\t")

  print(Outputname)

}
