###R code to read in the data for tumor volume and cell infiltration
#  the analysis parts are in analysis_vX.R
#
#reading the data about Tumor volume and cell infiltration
####################################################
#   we now move these below code to the second version v2.0
#	please see those there
#
#now let's starting connecting to the other data!!!
#  1. tumor volume
#  2. cell infiltration
#	
#	there are two ways to do this
#		1) using the extra information (variables)
#			as the supplementary variable to see where 
#			the cyotkine modes correlate with the supplementary variables
#		2) add these ones as predicting variables, to see how to 
#			combined together to separate/classify the groups (treatments)

library(tidyverse)
library(readxl)
library(dplyr)
library(magrittr)
library(here)

#load data...
setwd("/mnt/sb_ub1/share/Feng/VISTA/MC38/")
#tumor volume
volume.tumor<-read_excel("48 mice picked from EVH_011_MC38 061422.xlsx", sheet="Sheet2")

#clean up, removing some unnecessary rows/columns
names(volume.tumor)[c(1,5)]<-c("Mouse_ID", "Volume")
volume.tumor$Group<-sub(x=volume.tumor$Mouse_ID, 
		pattern="\\-[0-9]{1,2}", replacement="", fixed=F)
volume.tumor$Mouse_Num<-sub(x=volume.tumor$Mouse_ID, 
		pattern="[1-6]{1}\\-", replacement="", fixed=F)

#cells
cells.infiltration<-read_excel("TME FINAL_CD45_CD8_CD4.xlsx")

##clean up, removing some unnecessary rows/columns
#   first get rid of #1-9
cells.infiltration<-cells.infiltration[-c(1:10, 59:60),]
cells.infiltration<-cells.infiltration[,-c(4:7, 10:14)]

names(cells.infiltration)<-c("Mouse_ID", "Sample_CD45","Percent_CD45", 
		"Sample_CD8", "Percent_CD8","Sample_CD4", "Percent_CD4" )
cells.infiltration  %<>% mutate_at(vars(starts_with("Percent_")),as.numeric)

cells.infiltration %<>% mutate(Mouse_ID=str_replace(Mouse_ID,"#","")) %>%
	 mutate(Group_ID=str_replace(Mouse_ID,"\\-[0-9]{1,2}","")) %>%
	 mutate(Mouse_Num=str_replace(Mouse_ID,"[0-9]{1,2}\\-","")) 

README<-"data file containing infiltrating cells in tumor and tumor volume data\n"
README<-paste0(README,"the data were generated in data_tumor_cells.v1.0.R\n")

#save it and pass them to the other modules
save(file="data_tumor_cellInfil.RData", 
		cells.infiltration, volume.tumor, README)