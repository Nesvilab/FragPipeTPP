#' Translate TMTI headers to TPPR headers for 2DTPP data
#'
#' @param directory_of_interest
#'
#' @return a list of a dictionary of old and new lables, as well as a list of temperature and TMT conversion
#' @noRd
#'
#' @examples tmtiheader_to_tpprheaders("C:/FragPipeOutputfolder")
tmtiheader_to_tpprheaders2D <- function(directory_of_interest, configtemperatures){

  #print(configtemperatures)
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

  #Concentration lables to use fro TPP-R
  conc_vals_for_tppr <- list()
  tmt_tempconc_dict <- list()

  # Read annotation file
  lines <- readLines(annotationfile)
  for (line in lines) {
    line <- gsub("\n", "", line)
    splits <- strsplit(line, "\t")[[1]]

    # Obtain TMT label and experimental label (temperature_concentration)
    tmt_label <- splits[1]
    temp_conc <- splits[2]

    print(tmt_label)
    print(temp_conc)

    tmt_tempconc_dict[[tmt_label]] <- temp_conc

    # Add the experimental label as key and tmt label as value
    output_renaming_dict[[temp_conc]] <- tppr_labels_dict[[tmt_label]]

    # Change formatting to match the one needed by TPP-R config file
    config_tmt_label <- gsub("rel_fc_", "", tppr_labels_dict[[tmt_label]])
    exp_label_spl <- unlist(strsplit(temp_conc, "_"))
    config_temp_label <- exp_label_spl[2]
    conc_vals_for_tppr [[config_tmt_label]] <- as.numeric(config_temp_label)
  }

  print(tmt_tempconc_dict)
  print(names(tmt_tempconc_dict))

  for (labeltmt in tmt_tempconc_dict){
    print(labeltmt)
  }


  return(list(output_renaming_dict, conc_vals_for_tppr))
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


  #vectors to store congif file columns

  #Feb 16th: Replace first three with Compound, Experiment and Temperature
  #Create TPP-TPPR folder if it does not exist already
  tpprfolder <- file.path(fragpipefolder, "2DTPP-TPPR")

  if (dir.exists(tpprfolder)) {
    print("Saving results in 2DTPP-TPPR folder")
  } else {
    print("Creating 2DTPP-TPPR folder")
    dir.create(tpprfolder)
  }


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

        #print(path_components)

        #Adding to the columns of the future config file
        temperatures_replabel <- path_components[length(path_components)]

        #print(paste("temperatures_replabel",temperatures_replabel))




        conditionlabel <- unlist(strsplit(temperatures_replabel, "_"))

        temp_conce_labels_ls <- list()

        for (item in conditionlabel[1:2]) {
          #Temperature
          Experiment <- append(Experiment, temperatures_replabel)
          Temperature <- append(Temperature, item)
          Compound <- append(Compound, compound_val)
          outputfilename <- paste(temperatures_replabel, "txt", sep = ".")
          #Create new path
          Outputname <- file.path(tpprfolder, outputfilename)
          pathcol <- append(pathcol, Outputname )

          for (conce_val in concentrationlabels){
            #print(conce_val)
            temp_conc_label <- paste(item,conce_val, sep = "_")
            lscolumnstoget[[length(lscolumnstoget) + 1]] <- paste("X",temp_conc_label, sep="")

            #print(paste("temp_conc_label",temp_conc_label))

            #TODO: Gather temperatures and concentrations for config file




            }
          }




        #print(paste("lscolumnstoget",lscolumnstoget))
        # Selecting columns by name
        newdataframe <- tmtidata[, unlist(lscolumnstoget)]
        #print(class(newdataframe))
        #print(newdataframe$"Index")

       # print(colnames(newdataframe))

        #Normalization to lowest temperature (first two columns are protein id and qssm filtering criteria)
        newdataframe[, 3:ncol(newdataframe)] <-   newdataframe[, 3:ncol(newdataframe)] / newdataframe[,3]

        #print(newdataframe[, -1])

        #TODO: Use annotation file to replace the temperature_experiment with appropriate TMT -labels
        tpprheaders <- tmtiheader_to_tpprheaders2D(folder, Temperature)
        #print(tpprheaders)

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

        #(colnames(newdataframe))
        #print(tpprheaders[1])
        #Place correct headers
        names(newdataframe) <- unlist(tpprheaders[1])
        #print(newdataframe)







        # Save data frame as a tab-delimited text file
        #write.table(newdataframe, file = Outputname, sep = "\t", row.names = FALSE, quote = FALSE) - use base R (slower)
        data.table::fwrite(newdataframe, file = Outputname, sep = "\t")


      }}
  }

  print(Compound)
  print(Experiment)
  print(Temperature)

  #Create configuration file
  configurationdf <- data.frame(Compound, Experiment, Temperature, TMTsix, TMTsevenL, TMTsevenH, TMTeightL, TMTeightH, TMTnineL, TMTnineH, TMTtenL, TMTtenH, TMTelevenL, pathcol)
  names(configurationdf) <- c("Compound", "Experiment", "Temperature", "126","127L", "127H", "128L","128H",	"129L", "129H",	"130L",	"130H",	"131L", "Path")
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
