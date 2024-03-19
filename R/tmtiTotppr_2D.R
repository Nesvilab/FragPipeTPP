#' Translate TMTI headers to TPPR headers for 2DTPP data
#'
#' @param directory_of_interest: fragPipe directory
#' @param configtemperatures: A list fo temperatures
#'
#' @return a list of a dictionary of old and new labels, as well as a list of temperatures with lists of TMT anc Concentration conversion
#' @noRd
#'
#' @examples tmtiheader_to_tpprheaders("C:/FragPipeOutputfolder")
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
  output_renaming_dict <- list('Index' = 'Prot_ID', "NumberPSM" = "qusm")

  # Dictionary to translate TMT labels to TPP-R labels
  tppr_labels_dict <- list(
    "126" = "rel_fc_126", "127N" = "rel_fc_127L", "127C" = "rel_fc_127H",
    "128N" = "rel_fc_128L", "128C" = "rel_fc_128H", "129N" = "rel_fc_129L",
    "129C" = "rel_fc_129H", "130N" = "rel_fc_130L", "130C" = "rel_fc_130H",
    "131N" = "rel_fc_131L"
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
      config_tmt_label <- gsub("rel_fc_", "", tppr_labels_dict[[labeltmt]])
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


#' Convert FragPipe Output to TPP-R input (for 2DTPP analysis)
#'
#' @param fragpipefolder: Path to the folder where FragPipe saves results
#' @param experimentlables: A vector of stringd to determine what labels are use for different experiment types
#'
#' @return Path to the created TPP-R configuration file
#' @export
#'
#' @examples tmtitpttpr("C:/FragPipeOutputfolder", c("42_44","46_48","50_52","54_56","58_60","62_64"), c(0,0.005,0.05,0.5,2), c(2mg ATP))
tmtitotppr_2D <- function(fragpipefolder, experimentlabels, concentrationlabels, compound_val){

  #Variables to find the correct TMTI .tsv file
  tmtifile <- "tmt-report/ratio_protein_None.tsv"
  inputfile <-  file.path(fragpipefolder, tmtifile)

  # Open TSV file
  tmtidata <- read.delim(inputfile)

  #Extract column names
  headers <- colnames(tmtidata)

  #print(headers)

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

      lscolumnstoget <- list("Index", "NumberPSM")
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

        # Selecting columns by name
        print(paste("Extracting columns needed for experiment", temperatures_replabel))
        newdataframe <- tmtidata[, unlist(lscolumnstoget)]

        #Normalization to lowest temperature (first two columns are protein id and qssm filtering criteria)
        #newdataframe[, 3:ncol(newdataframe)] <-   newdataframe[, 3:ncol(newdataframe)] / newdataframe[,3]

        #Use annotation file to replace the temperature_experiment with appropriate TMT -labels
        tpprheaders <- tmtiheader_to_tpprheaders2D(folder, Temperature)

        #Extract the TMT-label to Celcius conversion
        configtempvals <- tpprheaders[2][[1]]

        #For each temperature extract the appropriate TMT lables use for the avialble concentration
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

        #Place correct headers
        names(newdataframe) <- unlist(tpprheaders[1])
        #print(newdataframe)


        # Save data frame as a tab-delimited text file
        #write.table(newdataframe, file = Outputname, sep = "\t", row.names = FALSE, quote = FALSE) - use base R (slower)
        data.table::fwrite(newdataframe, file = Outputname, sep = "\t")


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


#Test
twofragpipe <- "Z:/crojaram/TPP_Project/PXD012423/2DTPP/ATP_rep1/FP20-1_build23"
conc_labels <- c(0,0.005,0.05,0.5,2)
labels_exp <- c("42_44","46_48","50_52","54_56","58_60","62_64")
compound <-("ATP")
configtwo <- tmtitotppr_2D(twofragpipe,labels_exp,conc_labels,compound)

tpp2dResults <- TPP::analyze2DTPP(configTable = configtwo,
                                  compFc = TRUE,
                                  idVar = "Prot_ID",
                                  intensityStr = "rel_fc_",
                                  nonZeroCols = "qusm",
                                  methods = "doseResponse",
                                  createReport = "none",
                                  resultPath = list_of_data[2]
)
