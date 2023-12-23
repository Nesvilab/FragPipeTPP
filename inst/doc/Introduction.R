## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
library(FragPipeToTPPR, TPP)

## ----Path to local Frag Pipe folder-------------------------------------------
FPpath <- FragPipeToTPPR::dummyFragPipepath

## ----FragPipePath for building , include=FALSE--------------------------------
FPpath<- "C:/Users/crojaram/Documents/FragPipeToTPPR/inst/extdata/dummyFragPipeFolder"

## ----Experimental lables------------------------------------------------------
explabels <- FragPipeToTPPR::dummyexplabels

## ----TMT-Integrator (in FragPipe) ouput to TPP-R input------------------------
configpath <- FragPipeToTPPR::tmtitotppr(FPpath, explabels)

## ----ConfigPath for building , include=FALSE----------------------------------
configpath <- "C:/Users/crojaram/Documents/FragPipeToTPPR/inst/extdata/dummyFragPipeFolder/1DTPP-TPPR/TPP-TR_config.csv"

## ----Check path of previous output--------------------------------------------
print(configpath)

## ----Database for building , include=FALSE------------------------------------
hoomandb <- "C:/Users/crojaram/Documents/FragPipeToTPPR/inst/extdata/dummyFragPipeFolder/2023-07-12-decoys-reviewed-contam-UP000005640.fas"

## ----Melting Curve Normalization----------------------------------------------
FragPipeToTPPR::tpprNormOneDTPP(configpath,FPpath)

## ----TPP-R to TP-MAP----------------------------------------------------------

FragPipeToTPPR::tpprTotpmap(configpath,hoomandb)

