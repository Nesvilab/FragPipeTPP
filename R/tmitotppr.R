
#' Translate TMTI headers to TPPR headers
#'
#' @param directory_of_interest" Experiment folder inside FragPipe Output folder
#'
#' @return a list of a dictionary of old and new lables, as well as a list of temperature and TMT conversion
#' @export
#'
#' @examples tmtiheader_to_tpprheaders("C:/FragPipeOutputfolder"/Treatment)
tmtiheader_to_tpprheaders <- function(directory_of_interest){

  annotationfile <- ""
  files_in_specific_directory <- list.files(directory_of_interest)
  for (file in files_in_specific_directory){
    if (grepl("annotation.txt", file)) {
      annotationfile <-  file.path(directory_of_interest, file)

    }
  }
  # Initiate a list to rename column headers
  output_renaming_dict <- list('Index' = 'Prot_ID', "NumberPSM" = "qssm")

  # Dictionary to translate TMT labels to TPP-R labels
  tppr_labels_dict <- list(
    "126" = "rel_fc_126", "127N" = "rel_fc_127L", "127C" = "rel_fc_127H",
    "128N" = "rel_fc_128L", "128C" = "rel_fc_128H", "129N" = "rel_fc_129L",
    "129C" = "rel_fc_129H", "130N" = "rel_fc_130L", "130C" = "rel_fc_130H",
    "131N" = "rel_fc_131L"
  )

  #Temperature lables to use fro TPP-R
  temp_vals_for_tppr <- list()

  # Read annotation file
  lines <- readLines(annotationfile)
  for (line in lines) {
    line <- gsub("\n", "", line)
    splits <- strsplit(line, "\t")[[1]]

    # Obtain TMT label and experimental label (temperature_experimenttype-replicatw)
    tmt_label <- splits[1]
    exp_label <- splits[2]

    # Add the experimental label as key and tmt label as value
    output_renaming_dict[[exp_label]] <- tppr_labels_dict[[tmt_label]]

    # Change formatting to match the one needed by TPP-R config file
    config_tmt_label <- gsub("rel_fc_", "", tppr_labels_dict[[tmt_label]])
    exp_label_spl <- unlist(strsplit(exp_label, "_"))
    config_temp_label <- exp_label_spl[1]
    temp_vals_for_tppr[[config_tmt_label]] <- as.numeric(config_temp_label)
  }
  return(list(output_renaming_dict, temp_vals_for_tppr))
}




#' Convert FragPipe Output to TPP-R input (for 1DTPP analysis)
#'
#' @param fragpipefolder: Path to the folder where FragPipe saves results
#' @param experimentlables: A vector of stringd to determine what labels are use for different experiment types
#'
#' @return Path to the created TPP-R configuration file
#' @export
#'
#' @examples tmtitpttpr("C:/FragPipeOutputfolder", c("Treatment", "Vehicle"))
tmtitotppr <- function(fragpipefolder, experimentlables){

  #Variables to find the correct TMTI .tsv file
  tmtifile <- "tmt-report/ratio_protein_None.tsv"
  inputfile <-  file.path(fragpipefolder, tmtifile)

  # Open TSV file
  tmtidata <- read.delim(inputfile)

  #Extract column names
  headers <- colnames(tmtidata)

  # List all folders in the specified directory
  folders <- list.dirs(fragpipefolder, recursive = FALSE)

  #vectors to store congif file columns
  Experiment <- c()
  Condition <- c()
  Replicate <- c()
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
  pathcol <- c()

  #list to create comparison columns
  comparison <- list()

  #Find the experimental names
  for (folder in folders) {
    foundexplabel <- ""

    # Select experiment folders to extract the number for replicates and
    #look for those labels among the columns of the tsv file
    for (label in experimentlables) {
      lscolumnstoget <- list("Index", "NumberPSM")
      if (grepl(label, folder)) {
        combined_string <- paste("Folder: ", folder)

        # Splitting the file path using strsplit
        path_components <- unlist(strsplit(folder, "/"))
        #print(path_components[length(path_components)])

        #Obtain exp label
        foundexplabel <- path_components[length(path_components)]
        Experiment <- append(Experiment, foundexplabel)
        #remove replicate to set condition
        conditionlabel <- unlist(strsplit(foundexplabel, "_"))
        Condition <- append(Condition, conditionlabel[1])
        Replicate <- append(Replicate, conditionlabel[2])

        #Slice data frame according to the experiment and replicate
        for (item in headers) {
          if (grepl(foundexplabel, item)){
            #print(item)
            # Remove the first 'X' character from the string, which was added when using read.delim()
            lscolumnstoget[[length(lscolumnstoget) + 1]] <- item

          }
        }


        # Selecting columns by name
        newdataframe <- tmtidata[, unlist(lscolumnstoget)]
        print(class(newdataframe))
        #print(newdataframe$"Index")
        #Normalization to lowest temperature (first two columns are protein id and qssm filtering criteria)
        newdataframe[, 3:ncol(newdataframe)] <-   newdataframe[, 3:ncol(newdataframe)] / newdataframe[,3]

        #print(newdataframe[, -1])

        #TODO: Use annotation file to replace the temperature_experiment with appropriate TMT -labels
        tpprheaders <- tmtiheader_to_tpprheaders(folder)
        print(tpprheaders)
        #Extract the TMT-label to Celcius conversion
        configtempvals <- tpprheaders[2][[1]]


        TMTsix <- append(TMTsix, configtempvals$`126` )
        TMTsevenL <- append(TMTsevenL, configtempvals$`127L` )
        TMTsevenH <- append(TMTsevenH, configtempvals$`127H` )
        TMTeightL <- append(TMTeightL, configtempvals$`128L` )
        TMTeightH <- append(TMTeightH, configtempvals$`128H` )
        TMTnineL <- append(TMTnineL, configtempvals$`129L` )
        TMTnineH <- append(TMTnineH, configtempvals$`129H` )
        TMTtenL <- append(TMTtenL, configtempvals$`130L` )
        TMTtenH <- append(TMTtenH, configtempvals$`130H` )
        TMTelevenL <- append(TMTelevenL, configtempvals$`131L` )

        #Place correct headers
        names(newdataframe) <- unlist(tpprheaders[1])
        #print(newdataframe$"Prot_ID")

        outputfilename <- paste(foundexplabel, "txt", sep = ".")

        #Create TPP-TPPR folder if it does not exist already
        tpprfolder <- file.path(fragpipefolder, "1DTPP-TPPR")

        if (dir.exists(tpprfolder)) {
          print("Saving results in 1DTPP-TPPR folder")
        } else {
          print("Creating 1DTPP-TPPR folder")
          dir.create(tpprfolder)
        }

        #Create new path
        Outputname <- file.path(tpprfolder, outputfilename)
        pathcol <- append(pathcol, Outputname )


        # Save data frame as a tab-delimited text file
        #write.table(newdataframe, file = Outputname, sep = "\t", row.names = FALSE, quote = FALSE) - use base R (slower)
        data.table::fwrite(newdataframe, file = Outputname, sep = "\t")

      }
    }

  }

  #Create configuration file
  configurationdf <- data.frame(Experiment, Condition, Replicate, TMTsix, TMTsevenL, TMTsevenH, TMTeightL, TMTeightH, TMTnineL, TMTnineH, TMTtenL, TMTtenH, TMTelevenL, pathcol)
  names(configurationdf) <- c("Experiment", "Condition", "Replicate", "126","127L", "127H", "128L","128H",	"129L", "129H",	"130L",	"130H",	"131L", "Path")
  configsavepath <- file.path(fragpipefolder, "1DTPP-TPPR", "TPP-TR_config.csv")
  write.csv(configurationdf, configsavepath, row.names = FALSE)

  return(configsavepath)


}
