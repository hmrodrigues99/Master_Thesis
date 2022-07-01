#!/bin/bash

#This script receives a trimmomatic_quality.csv, a STAR_log.csv and a featureCounts_log.csv of a given BioProject and computes its BioProject_Statistics.csv

export PATH=$PATH:$HOME/data/$1

paste trimmomatic_log.csv STAR_log.csv featureCounts_log.csv > ${1}_Statistics.csv
#rm trimmomatic_log.csv
#rm STAR_log.csv
#rm featureCounts_log.csv
