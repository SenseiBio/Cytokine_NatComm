#R code, accessory functions

#this is the function to detect CCL names and non-CCL cytokine names
# then switch the non-CCL names to CCL ones.
# here we rely on the input and only determine its CCL names if they show in the input 
# names. No extra data base for now.
#
#'@param cnames a vector of cytokine names 
#'@return the list of switched cytokine names.
#'@todo need to further testing about the duplicated cytokines how to avoid mis-match
#	.......
UnifyCytokineNames<-function(cnames){
###take care the names first, getting ccl names when they are available.
#	get rid of the duplicated whenever they are.
# take care the name duplication, and also need to get ccl name
#cnames<-colnames(all_dat)[-c(1,2)]
# cnames<-stats.mc38$cytokine # for testing 

#remembering the duplicates, denoted by ".2" or ".1" (anything with .#)
#
duplicated.names<-sub(x=cnames,"^[a-zA-Z0-9\\-]+","")
duplicated.names<-sub(x=duplicated.names,"^\\/[a-zA-Z0-9]+","")
duplicated.names<-sub(x=duplicated.names,"^\\([a-zA-Z0-9\\-]+\\)","")

cnames.ccl<-sub(x=cnames,"^[a-zA-Z0-9\\-]+","")
cnames.ccl<-sub(x=cnames.ccl,"^\\([a-zA-Z0-9\\-]+\\)","")
cnames.ccl<-sub(x=cnames.ccl,"^\\/","")
cnames.ccl<-sub(x=cnames.ccl,"^\\.[0-9]+","")

cnames.nccl<-sub(x=cnames, "\\/[a-zA-Z0-9\\-]+$","")
cnames.nccl<-sub(x=cnames.nccl, "\\.[0-9]+$","")

for(i in 1:length(cnames.ccl))
{
	if(cnames.ccl[i] !="")
	{
		#replace with ccl names
		cnames[cnames.nccl==cnames.nccl[i]]=cnames.ccl[i]
	}
}
#colnames(all_dat)[-c(1:2)]<-cnames
return(cnames);
}