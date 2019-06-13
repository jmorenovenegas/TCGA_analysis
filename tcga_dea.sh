#!/bin/bash

echo TCGA TCNB vs Luminal A Differential Expression Analysis workflow

if [ $# == 0 ]; then

echo 'You must include a parameter:'
echo 'All if you want to get the data and perform the analysis'
echo 'Data if you only want to get the data'
echo 'Analysis if you already have the data available and you only want to perform the analysis'

fi

if [ $1 == 'All' ]  || [ $1 == 'Data' ]; then
	
echo Getting TCGA data

Rscript modulesR/1_getData.R  

echo Done

echo Subsetting data

Rscript modulesR/2_subsettingData.R

echo Done

fi

if [ $1 == 'All' ]  || [ $1 == 'Analysis' ] ; then

Rscript modulesR/3_DEA.R TNBC

echo Done

echo Clustering and survival analysis

Rscript modulesR/5_Clustering_and_survival_analysis.R

echo Done

fi

exit

