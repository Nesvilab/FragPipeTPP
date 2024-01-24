#' Extract protein description for TP-MAP compatability
#'
#' @param filepathdb
#'
#' @returnlist with protein uniprot accession numbers and corresponding description line
#' @noRd
#'
#' @examples fastaparser("C:/protein.fas")
fastaparser <- function(filepathdb) {
  # Function to obtain protein information string whole from FASTA file
  # @param filepathdb: path, FASTA file location
  # @return: list, key: Protein Accession; value: string

  # Intialize list
  identifiers <- list()

  # Parse FASTA file
  con <- file(filepathdb, "r")
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    # Only pay attention to protein headers
    if (startsWith(line, ">")) {
      # Look only at the string after "|"
      modregex <- regexpr("\\|", line, ignore.case = TRUE)
      proteinidindeces <- regmatches(line, modregex)[[1]]

      # What is the protein id string
      keyproteinid <- substr(line, proteinidindeces[1] + 1, proteinidindeces[2] - 1)

      # Add protein id as key and whole description line as value
      identifiers[[keyproteinid]] <- line
    }
  }
  close(con)

  return(identifiers)
}

#' Organizing dataframe column headers for TP-MAP compatibility
#'
#' @param colstoarrenge: what columns to change
#'
#' @return list of lists with the connections between old and new names.
#' @noRd
#'
#' @examples column_titlesformatter(list)
column_titlesformatter <- function(colstoarrenge) {

  renaming_dict <- list()

  # List of the column names to arrange
  arrangedcols <- c(colstoarrenge[1], colstoarrenge[length(colstoarrenge)], colstoarrenge[-c(1, length(colstoarrenge))])

  # Change TMTI temperature column headers to TP-MAP compatible ones

  #Change TMTI headers to TPMAP compatible ones
  for (item in arrangedcols) {

    #print(item)
    newcolumntitle <- ""

    if (item == "Description") {

      newcolumntitle <- "Description"
    } else if (item == "Index") {

      newcolumntitle <- "Accession"
    } else {

      #newcolumntitle <- substr(item, 2, nchar(item))
      newcolumntitle <- sub("^.", "Ref_", item)

    }

    #print(newcolumntitle)


    renaming_dict[[item]] <- newcolumntitle
  }

  return(list(arrangedcols, renaming_dict))
}

#' Putting together the output data frame to be input into TP-MAP for 2DTPP data
#'
#' @param proteindb: path to FASTA database
#' @param TMTIPath: path to TMTI file
#'
#' @return data.frame
#' @noRd
#'
#' @examples dataframecreator("C:/protein.fas", "C:/FragPipeOutputfolder/tmt-report")
dataframecreator <- function(proteindb, TMTIPath) {
  # Putting together the output data frame to be input into TP-MAP
  # @param proteindb: path to FASTA database
  # @param TMTIPath: path to TMTI file
  # @return: data.frame

  fastafile <- fastaparser(proteindb)
  abundanceDT <- read.table(TMTIPath, sep = "\t", header = TRUE)
  setwd(dirname(TMTIPath))

  # Copy original to prevent changes to the main DT
  outputDT <- abundanceDT

  # What columns to extract from the tmt-report ratio file
  thecolumns <- colnames(outputDT)
  referenceindex <- which(thecolumns == "ReferenceIntensity")

  # First column and the ones that contain abundance information, excluding reference
  outputcolumns <- c(thecolumns[(referenceindex + 1):length(thecolumns)], thecolumns[1])

  descriptionoutlist <- character(0)

  # Gathering description lines from FASTA file based on protein index from TMTI file
  for (item in outputDT$Index) {
    if (item %in% names(fastafile)) {
      tpmapline <- fastafile[[item]]
      tpmapline <- gsub("\n", "", tpmapline)
      tpmapline <- gsub(">", "", tpmapline)
      descriptionoutlist <- c(descriptionoutlist, tpmapline)
    }
  }

  # Selecting the columns to output and adding the description info for each protein
  finalDT <- outputDT[, outputcolumns]
  finalDT$Description <- descriptionoutlist

  # Format selected columns to be TP-MAP compatible
  oldcols <- names(finalDT)
  newcols <- column_titlesformatter(oldcols)[[1]]
  finalDT <- finalDT[, newcols]

  # Rename Index header to Accession
  names(finalDT)[names(finalDT) == "Index"] <- "Accession"

  colnames(finalDT) <- column_titlesformatter(oldcols)[[2]]

  # Replace NA with zero
  finalDT[is.na(finalDT)] <- 0

  return(finalDT)
}


#' Performs 2DTPP data file conversion FragPipe outputs to TP-MAP input
#'
#' @param resultPath: FragPipe output folder.
#' @param fastadatabsefile: .fas file of the protein database
#'
#' @return void
#' @export
#'
#' @examples
tmtiTotpmap_2D <- function(resultPath, fastadatabsefile){

  #Variables to find the correct TMTI .tsv file
  tmtifile <- "tmt-report/ratio_protein_None.tsv"
  tmtifilepath <-  file.path(resultPath, tmtifile)


  # Calling dataframecreator function
  tpmapDT <- dataframecreator(fastadatabsefile, tmtifilepath)
  #print(tpmapDT)

  # Create folder to save 2DTPP Results
  twodtpp_output <- file.path(resultPath, "2DTPP")
  if (!file.exists(twodtpp_output)) {
    dir.create(twodtpp_output)
  }
  setwd(twodtpp_output)

  # Save tpmapDT to a CSV file
  write.table(tpmapDT, file = "2DTPP.tpmap.txt", sep = "\t", row.names = FALSE)
}
