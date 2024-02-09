## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
library(FragPipeTPP, TPP)

## ----Path to local Frag Pipe folder-------------------------------------------
FPpath <- FragPipeTPP::dummyFragPipepath

## ----FragPipePath for building , include=FALSE--------------------------------
FPpath<- "C:/Users/crojaram/Documents/FragPipeTPP/inst/extdata/dummyFragPipeFolder"

## ----Experimental lables------------------------------------------------------
explabels <- FragPipeTPP::dummyexplabels

## ----TMT-Integrator (in FragPipe) ouput to TPP-R input------------------------
configpath <- FragPipeTPP::tmtitotppr(FPpath, explabels)

## ----ConfigPath for building , include=FALSE----------------------------------
configpath <- "C:/Users/crojaram/Documents/FragPipeTPP/inst/extdata/dummyFragPipeFolder/1DTPP-TPPR/TPP-TR_config.csv"

## ----Check path of previous output--------------------------------------------
print(configpath)

## ----Database for building , include=FALSE------------------------------------
hoomandb <- "C:/Users/crojaram/Documents/FragPipeTPP/inst/extdata/dummyFragPipeFolder/2023-07-12-decoys-reviewed-contam-UP000005640.fas"

## ----Melting Curve Normalization----------------------------------------------
FragPipeTPP::tpprNormOneDTPP(configpath,FPpath)

## ----TPP-R to TP-MAP----------------------------------------------------------

FragPipeTPP::tpprTotpmap(configpath,hoomandb)

