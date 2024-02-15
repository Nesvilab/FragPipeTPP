#' Convert FragPipe Output to TPP-R input (for 2DTPP analysis)
#'
#' @param fragpipefolder: Path to the folder where FragPipe saves results
#' @param experimentlables: A vector of stringd to determine what labels are use for different experiment types
#'
#' @return Path to the created TPP-R configuration file
#' @export
#'
#' @examples tmtitpttpr("C:/FragPipeOutputfolder", c("42_44","46_48","50_52","54_56","58_60","62_64"), c(0,0.005,0.05,0.5,2))
tmtitotppr_2D <- function(fragpipefolder, experimentlabels, concentrationlabels){

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
    for (label in experimentlabels) {

      lscolumnstoget <- list("Index", "NumberPSM")
      if (grepl(label, folder)) {

        # Splitting the file path using strsplit
        path_components <- unlist(strsplit(folder, "/"))

        #print(path_components)

        temperatures_replabel <- path_components[length(path_components)]

        Experiment <- append(Experiment, temperatures_replabel)

        #print(Experiment)

        conditionlabel <- unlist(strsplit(temperatures_replabel, "_"))

        #Adding to the columns of the future config file
        Condition <- append(Condition, conditionlabel[1])
        Condition <- append(Condition, conditionlabel[2])
        #Double up due to the precense of two temperatures
        Replicate <- append(Replicate, conditionlabel[3])
        Replicate <- append(Replicate, conditionlabel[3])

        #print(conditionlabel)

        for (item in conditionlabel[1:2]) {
          print(item)
          for (conce_val in concentrationlabels){
            temp_conc_label <- paste(item,conce_val, sep = "_")


            #Slice data frame according to the experiment and replicate
            for (item in headers) {
              if (grepl(temp_conc_label, item)){
                #print(item)
                # Remove the first 'X' character from the string, which was added when using read.delim()
                lscolumnstoget[[length(lscolumnstoget) + 1]] <- item

              }
            }

          }}

        print(lscolumnstoget)


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

      }}
  }
}



