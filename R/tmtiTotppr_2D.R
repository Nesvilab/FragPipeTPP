#' Convert FragPipe Output to TPP-R input (for 2DTPP analysis)
#'
#' @param fragpipefolder: Path to the folder where FragPipe saves results
#' @param experimentlables: A vector of stringd to determine what labels are use for different experiment types
#'
#' @return Path to the created TPP-R configuration file
#' @export
#'
#' @examples tmtitpttpr("C:/FragPipeOutputfolder", c("42_44","46_48","50_52","54_56","58_60","62_64"), c(0,0.005,0.05,0.5,2), c(2mg ATP))
tmtitotppr_2D <- function(fragpipefolder, experimentlabels, concentrationlabels, compound){

  #Variables to find the correct TMTI .tsv file
  tmtifile <- "tmt-report/ratio_protein_None.tsv"
  inputfile <-  file.path(fragpipefolder, tmtifile)

  # Open TSV file
  tmtidata <- read.delim(inputfile)

  #Extract column names
  headers <- colnames(tmtidata)

  print(headers)

  # List all folders in the specified directory
  folders <- list.dirs(fragpipefolder, recursive = FALSE)

  #print(folders)


  #vectors to store congif file columns

  #Feb 16th: Replace first three with Compound, Experiment and Temperature

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
        Compound <- append(compound, compound)
        temperatures_replabel <- path_components[length(path_components)]
        Experiment <- append(Experiment, temperatures_replabel)
        conditionlabel <- unlist(strsplit(temperatures_replabel, "_"))

        #print(conditionlabel)
        temp_conce_labels_ls <- list()
        for (item in conditionlabel[1:2]) {
          #Temperature
          #print(item)
          for (conce_val in concentrationlabels){
            #print(conce_val)
            temp_conc_label <- paste(item,conce_val, sep = "_")
            lscolumnstoget[[length(lscolumnstoget) + 1]] <- paste("X",temp_conc_label, sep="")
            #print(temp_conc_label)
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
        tpprheaders <- tmtiheader_to_tpprheaders(folder)
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
        print(newdataframe)

        outputfilename <- paste(temperatures_replabel, "txt", sep = ".")


      }}
  }
}



