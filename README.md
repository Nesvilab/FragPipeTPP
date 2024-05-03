# Thermal Protein Profiling (TPP) Data Analysis with FragPipe: FragPipeTPP
###### Update: May 3rd, 2024

## What is TPP?
TPP is an unbiased method to search for drug targets by monitoring protein stability in the presence and absence of the drug compound of interest. Different samples at different temperatures are prepared and them label with [Tandem Mass Tags (TMT)](https://pubs.acs.org/doi/10.1021/ac0262560). Next, samples are pooled together and analyzed by typical mass spectrometry-based proteomics. In the links below users ar encouraged to review the Nature Protocol method and a recent overview of TPP to understand how the data is analyzed and how protein melting curves are obtained. 

  - 2015 [Nature Methods TPP protocol](https://www.nature.com/articles/nprot.2015.101)
  - Recent [small overview](https://pubmed.ncbi.nlm.nih.gov/36368297/)
  - Data use in this package was analyzed by FragPipe. Raw data was obtained from [a TPP study](https://www.nature.com/articles/s41467-019-09107-y), where the effects of ATP in protein solubility were determined.     
  
Comparing a protein melting curve with and without drug is important to determine if the compound had a destabilizing effect on the protein (melting temperature decrease) or the opposite stabilizing effect (melting temperature increase) or none.

## What is FragPipe?
- It is a set of computation tools capable of protein identification and quantification from bottom-up proteomics data. One of its modules is capable of analyzing TMT datasets, such as the one produced in TPP experiments. FragPipe is a very fast and robust pipeline, but there are no FragPipe-compatible TPP tools (until now). To learn how to use FragPipe see the link below:
  
  - [FragPipe website](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html])
  
- For downstream TPP analysis of FragPipe data, this package makes it possible to use both TPP-R and TP-MAP.
  
## What is TPP-R?
It is an R package available through Bioconductor capable of analyzing both 1D (temperature change) and 2D (temperate and drug concentration changes) TPP datasets. It produces melting curves for each protein as well as metrics on the fittings of the curves. For more details see the link below:

  - [TPP-R at Bioconductor](https://bioconductor.org/packages/release/bioc/html/TPP.html)

## What is TP-MAP?
It is a GUI-based TPP analysis program also available for both 1DTPP and 2DTPP data. Melting curves (1DTPP) or matrices (2DTPP) can be visualized in TP-MAP. Protein-protein interactions analysis using STRING are also possible. for more information visit:

  - [2022 bioRxiv](https://www.biorxiv.org/content/10.1101/2021.02.22.432361v2)
  
  
### Package Installation

To analyze 1DTPP data produced by FragPipe, output files can be converted to input files compatible with both the R-Package TPP-R and the Java GUI program TP-MAP. Briefly, the paths to the FragPipe output folder, the protein database and labels for control and treatment samples are the only variables needed
 by the package. If there is no previous experience using R or RStudio, it is recommended to visit [R for beginners](https://education.rstudio.com/learn/beginner/).

In RStudio Console:

### 1. Make sure renv package is installed. This package makes sure all needed packages are installed. To install run: 

#### install.packages("renv")

If successful a message will appeared saying: "The downloaded binary packages are in", and it will indicate a path in the user's local machine.

### 2. Install TPP and FragPipeTPP packages

#### renv::install("bioc::TPP")

Below is a successful FragPipeTPP package installation:

The following package(s) will be installed:
- TPP [3.28.0]
These packages will be installed into "C:/Users/caror/AppData/Local/R/win-library/4.3".

Do you want to proceed? [Y/n]: Y

##### Installing packages --------------------------------------------------------
- Installing TPP ...                            OK [copied from cache in 0.37s]
Successfully installed 1 package in 0.39 seconds.

#### renv::install("Nesvilab/FragPipeTPP")

Below is a successful FragPipeTPP package installation:

###### Downloading packages -------------------------------------------------------
- Downloading FragPipeTPP from GitHub ...    OK [427.5 Kb in 0.67s]
Successfully downloaded 1 package in 1.1 seconds.

The following package(s) will be installed:
- FragPipeTPP [Nesvilab/FragPipeTPP@main]
These packages will be installed into "C:/Users/caror/AppData/Local/R/win-library/4.3".

Do you want to proceed? [Y/n]: Y

###### Installing packages --------------------------------------------------------
- Installing FragPipeTPP ...                 OK [built from source and cached in 2.1s]
Successfully installed 1 package in 2.4 seconds.

### 3. Load FragPipeTPP

#### library(FragPipeTPP)



**The Package is ready to use!**

There are three functions to use:

- **tmitotppr** (FragPipe to TPP-R)

- **tpprNormOneDTPP** (TPPR all melting curve normalization only - required to input 1DTPP data into TP-MAP)

- **tpprTotpmap** (from TPP-R normalized FragPipe data to TP-MAP input file)
  
- **tmtiTotpmap_2D** (FragPipe to TP-MAP conversion of 2DTPP data)

- **TwoDConversion** (FragPipe to TP-MAP conversion of 2DTPP data)

## For information and to see sample code results use browseVignettes("FragPipeTPP") in the RStudio Console.
