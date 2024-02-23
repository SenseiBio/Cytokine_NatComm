#processing msd/luminex data
#

##   to download gmt pathways from gsea broad institute
## http://www.gsea-msigdb.org/gsea/msigdb/index.jsp
### to download http://www.gsea-msigdb.org/gsea/downloads.jsp
# 
#====================\
#https://bioinformaticsbreakdown.com/how-to-gsea/
#
####note: the cytokine symbols (entrez) was referenced to 
## a cytokine registry database (flat file) at 
#			https://www.immport.org/resources
#install.packages("readxl")

library(readxl)
library(dplyr)
library(here)

#########section to estimate the detection limits for "OOR<" values
###Note: we estimate by eyes for the detetion limits
##
dl_Th17<-c(6, 105,1.5,15,
			7,6.5,65,10,
			10,1.7)
dl_Th17_cytokines<-c("IL-17F", "IL-21","IL-22","IL-23",
					"IL-25","IL-27","IL-31","IL-33",
					"CD40L","MIP-3a")
					
dl_23plex<-c(2.1, 4.5,5.2,1.2,
			0.55,2.1,1.5,0.5,
			0.65,1.0,0.85,3.0,
			4.5,15,4.5,11,
			0.55,2.1,17,0.5,
			3,3.1,3.5)
dl_23plex_cytokines<-c("Eotaxin", "G-CSF","GM-CSF","IFN-g",
						"IL-1a","IL-1b","IL-2","IL-3",
						"IL-4","IL-5","IL-6","IL-9",
						"IL-10","IL-12(p40)","IL-12(p70)","IL-13",
						"IL-17A","KC","MCP-1","MIP-1a",
						"MIP-1b","RANTES","TNF-a")					
						
dl_31plex<-c(0.55,16, 4.5,0.11,
				1.5,3,0.15,0.17,
				4.5,1.0,4.5,0.22,
				0.4,0.6,0.85,2.1,
				7.5,1.7,5.2,0.15,
				0.2,0.19,0.2,3.0,
				2.2,3.0,5.5,0.17,
				4.5,2.7,1.7)
dl_31plex_cytokines<-c("BCA-1/CXCL13","CTACK/CCL27","ENA-78/CXCL5","Eotaxin/CCL11",
					"Eotaxin2/CCL24","Fractalkine/CX3CL1","GM-CSF","I-309/CCL1",
					"I-TAC/CXCL11","IFN-g","IL-1b","IL-2",
					"IL-4","IL-6","IL-10","IL-16",
					"IP-10/CXCL10","KC/CXCL1","MCP-1/CCL2","MCP-3/CCL7",
					"MCP-5/CCL12","MDC/CCL22","MIP-1a/CCL3","MIP-1b/CCL4",
					"MIP-3a/CCL20","MIP-3b/CCL19","RANTES/CCL5","SCYB16/CXCL16",
					"SDF-1a/CXCL12","TARC/CCL17","TNF-a")

#setwd("E:\\feng\\personal_record\\cv\\jobSearch2\\SenseiBio\\Cytokine_Tim")
#setwd("C:\\Feng\\Vista\\Cytokines")
cdata_Th17 <- read_xlsx(here("Data","EVH_011 Mu Cytokine Th17 10-plex Concentrations 01Sep22 ASC.xlsx"),
		skip=6, sheet=1,col_names=T,progress=T)
		
		#get rid of extra suffix for the cytokine numbers on header
colnames(cdata_Th17)<-sub(x=colnames(cdata_Th17),"\\s\\([0-9A-Za-z]+\\)","")
cdata_Th17$Group<-paste0(sub(x=cdata_Th17$"Mouse ID", "\\-[0-9a-zA-Z]+$",""), 
			sub(x=cdata_Th17$"Mouse ID", "^[0-9]+\\-[0-9]+",""))
#now doing the numbers
#first get rid of the * in the front of the numbers
cdata_Th17_temp<- cdata_Th17 %>% apply( 2, FUN=function(x){sub(x=x,"*","",fixed=T)}) #%>%
				 #apply(2, FUN=function(x){1})
scale_factor<-0.667		#<-used to get rid of OOR data, under the detection limit 
cdata_Th17_temp <-as.data.frame(cdata_Th17_temp)		 
for(i in 3:dim(cdata_Th17_temp)[2])
{
	cid<-colnames(cdata_Th17)[i]
	cdata_Th17_temp[,i]<-sub(x=cdata_Th17_temp[,i],"OOR <", replacement=dl_Th17[which(dl_Th17_cytokines==cid)]*scale_factor,fixed=T) 
	cdata_Th17_temp[,i]<-as.numeric(cdata_Th17_temp[,i])
}
cdata_Th17<-cdata_Th17_temp
#done 
	

# 31 plex panel	
cdata_Cyto31 <- read_xlsx(here("Data","EVH_011 Mu Cytokine 31plex Concentrations 13Sep22 ASC.xlsx"),
		skip=6, sheet=1,col_names=T,progress=T, n_max=74)
colnames(cdata_Cyto31)[c(1,2)]=cdata_Cyto31[1,c(1,2)]
colnames(cdata_Cyto31)[2]="Mouse ID"
cdata_Cyto31 =cdata_Cyto31[-1,]
colnames(cdata_Cyto31) <- colnames(cdata_Cyto31) %>% sub(pattern="\\s\\([0-9A-Za-z]+\\)",replacement="") %>%
								sub(pattern="Mo ",replacement="", fixed=T)
cdata_Cyto31$Group<-paste0(sub(x=cdata_Cyto31$"Mouse ID", "\\-[0-9a-zA-Z]+$",""), 
			sub(x=cdata_Cyto31$"Mouse ID", "^[0-9]+\\-[0-9]+",""))

#now doing the numbers
#first get rid of the * in the front of the numbers
cdata_Cyto31_temp<- cdata_Cyto31 %>% apply( 2, FUN=function(x){sub(x=x,"*","",fixed=T)}) #%>%
				 #apply(2, FUN=function(x){1})
	
cdata_Cyto31_temp <-as.data.frame(cdata_Cyto31_temp)		 
for(i in 3:dim(cdata_Cyto31_temp)[2])
{
	cat ("doing index ", i, "..\n")
	cid<-colnames(cdata_Cyto31)[i]
	cdata_Cyto31_temp[,i]<-sub(x=cdata_Cyto31_temp[,i],"OOR <", replacement=dl_31plex[which(dl_31plex_cytokines==cid)]*scale_factor,fixed=T) 
	cdata_Cyto31_temp[,i]<-as.numeric(cdata_Cyto31_temp[,i])
}
cdata_Cyto31<-cdata_Cyto31_temp








		
cdata_Cyto23 <- read_xlsx(here("Data","EVH_011 Mu Cytokine 23plex Concentrations 07Sep22 ASC.xlsx"),
		skip=7, sheet=1,col_names=T,progress=T#, n_max=74
		)
#colnames(cdata_Cyto23)[c(1,2)]=cdata_Cyto23[1,c(1,2)]
colnames(cdata_Cyto23)[c(1:2)]=c("Group","Mouse ID")
#cdata_Cyto23 =cdata_Cyto23[-1,]
colnames(cdata_Cyto23) <- colnames(cdata_Cyto23) %>% sub(pattern="\\s\\([0-9A-Za-z]+\\)",replacement="") %>%
								sub(pattern="Mo ",replacement="", fixed=T)
cdata_Cyto23$Group<-paste0(sub(x=cdata_Cyto23$"Mouse ID", "\\-[0-9a-zA-Z]+$",""), 
			sub(x=cdata_Cyto23$"Mouse ID", "^[0-9]+\\-[0-9]+",""))
#now doing the numbers
#first get rid of the * in the front of the numbers
cdata_Cyto23_temp<- cdata_Cyto23 %>% apply( 2, FUN=function(x){sub(x=x,"*","",fixed=T)}) #%>%
				 #apply(2, FUN=function(x){1})
	
cdata_Cyto23_temp <-as.data.frame(cdata_Cyto23_temp)		 
for(i in 3:dim(cdata_Cyto23_temp)[2])
{
	cat ("doing index ", i, "..\n")
	cid<-colnames(cdata_Cyto23)[i]
	cdata_Cyto23_temp[,i]<-sub(x=cdata_Cyto23_temp[,i],"OOR <", 
				replacement=dl_23plex[which(dl_23plex_cytokines==cid)]*scale_factor,fixed=T) 
	cdata_Cyto23_temp[,i]<-as.numeric(cdata_Cyto23_temp[,i])
}
cdata_Cyto23<-cdata_Cyto23_temp

#now save the data raw data
save(file=here("Data","cytokine_data_raw_MC38.RData"), cdata_Cyto23, cdata_Cyto31, cdata_Th17)

#now we need to clean up the data
#




#now doing the pheatmap for
library(pheatmap)
#rotate the matrix to have row as cytokines and col as mouse
cdata_Th17_mat<-as.matrix(cdata_Th17[,c(3:dim(cdata_Th17)[2])])
cdata_Th17_mat<-log(t(cdata_Th17_mat))
colnames(cdata_Th17_mat)<-cdata_Th17$"Mouse ID"
#annotation
cdata_Th17$Group<-factor(cdata_Th17$Group, levels=c("1","2","5","6","8","9","1b","2b","6b","9b"))
mouse_group<-as.data.frame(cdata_Th17[,1])

rownames(mouse_group)<-cdata_Th17$"Mouse ID"
colnames(mouse_group)<-"Mouse Group"

gaps_col<-aggregate(cdata_Th17$"Mouse ID", by=list(cdata_Th17$Group), length)

h17<-pheatmap(cdata_Th17_mat,scale="none",  cluster_cols = F, cluster_rows=T,
			show_colnames = F, annotation_col=mouse_group, gaps_col=c(48)
			)
png(file="Th17_cytokine_panel.png", width=800, height=300)
h17
dev.off()

#rotate the matrix to have row as cytokines and col as mouse
cdata_Cyto31_mat<-as.matrix(cdata_Cyto31[,c(3:dim(cdata_Cyto31)[2])])
cdata_Cyto31_mat<-log(t(cdata_Cyto31_mat))
colnames(cdata_Cyto31_mat)<-cdata_Cyto31$"Mouse ID"
#annotation
cdata_Cyto31$Group<-factor(cdata_Cyto31$Group, levels=c("1","2","5","6","8","9","2b","8b","9b"))
mouse_group_31<-as.data.frame(cdata_Cyto31[,1])

rownames(mouse_group_31)<-cdata_Cyto31$"Mouse ID"
colnames(mouse_group_31)<-"Group"

#gaps_col<-aggregate(cdata_Cyto31$"Mouse ID", by=list(cdata_Cyto31$Group), length)

c31<-pheatmap(cdata_Cyto31_mat,scale="none",  cluster_cols = F, cluster_rows=T,
			show_colnames = F, #border_color = NA,
			annotation_col=mouse_group_31
			, gaps_col=c(48)
			)
png(file="Cyto31_cytokine_panel.png", width=800, height=900)
c31
dev.off()

#rotate the matrix to have row as cytokines and col as mouse
cdata_Cyto23_mat<-as.matrix(cdata_Cyto23[,c(3:dim(cdata_Cyto23)[2])])
cdata_Cyto23_mat<-log(t(cdata_Cyto23_mat))
colnames(cdata_Cyto23_mat)<-cdata_Cyto23$"Mouse ID"
#annotation
cdata_Cyto23$Group<-factor(cdata_Cyto23$Group, levels=c("1","2","5","6","8","9"))#,"2b","8b","9b"))
mouse_group_23<-as.data.frame(cdata_Cyto23[,1])

rownames(mouse_group_23)<-cdata_Cyto23$"Mouse ID"
colnames(mouse_group_23)<-"Group"

#gaps_col<-aggregate(cdata_Cyto31$"Mouse ID", by=list(cdata_Cyto31$Group), length)

c23<-pheatmap(cdata_Cyto23_mat,scale="none",  cluster_cols = F, cluster_rows=T,
			show_colnames = F,  border_color = "grey60",annotation_col=mouse_group_23
			#, gaps_col=c(48)
			)
png(file="Cyto23_cytokine_panel.png", width=800, height=630)
c23
dev.off()


#put things together
# take care of *B data, get average of b and original sample (technical replicates)
cdata_Cyto31_log<-log(cdata_Cyto31[,-c(1,2)])
cdata_Cyto31_log<-cbind(cdata_Cyto31[,c(1,2)], cdata_Cyto31_log)

temp<-cdata_Cyto31_log
temp<-temp[grep(temp$Group, pattern="^[289]"),]
temp$"Mouse ID"<-sub(x=temp$"Mouse ID", pattern="b","")

#loop through columns to do averaging.
means<-data.frame()
for( i in 3:dim(temp)[2])
{
	if(i==3)
	{
		means=aggregate(temp[,i], by=list(temp$"Mouse ID"), mean)
		rownames(means)<-means[,1]
		colnames(means)[2]<-colnames(temp)[i]
	}
	else
	{	a<-aggregate(temp[,i], by=list(temp$"Mouse ID"), mean)
		rownames(a)<-a[,1]
		colnames(a)[2]<-colnames(temp)[i]
		means<-cbind(means,a[rownames(means),2])
	}
}

colnames(means)[-1]<-colnames(temp)[-c(1:2)]
colnames(means)[1]<-"Mouse ID"
means$"Group"<-sub(x=means$"Mouse ID", pattern="\\-[0-9]+$","")
cdata_Cyto31_clean<- cdata_Cyto31_log[grep(cdata_Cyto31_log$Group, pattern="[156]"),]
cdata_Cyto31_clean<- rbind(cdata_Cyto31_clean,means)
cdata_Cyto31_clean$Group<-factor(cdata_Cyto31_clean$Group, levels=unique(cdata_Cyto31_clean$Group))
rownames(cdata_Cyto31_clean)<-cdata_Cyto31_clean$"Mouse ID"
# Th17 
cdata_Th17_log<-log(cdata_Th17[,-c(1,2)])
cdata_Th17_log<-cbind(cdata_Th17[,c(1,2)], cdata_Th17_log)
temp<-cdata_Th17_log 
temp<-temp[-grep(temp$Group, pattern="b$"),]  #it seems in this case Th17 is not technical replicates, we get rid of them and only use the first set of data.
cdata_Th17_clean<-temp
cdata_Th17_clean$Group<-factor(cdata_Th17_clean$Group, levels=unique(cdata_Th17_clean$Group))
rownames(cdata_Th17_clean)<-cdata_Th17_clean$"Mouse ID"

#nothing to do with cyot23. since there is * bsamples.
cdata_Cyto23_log<-log(cdata_Cyto23[,-c(1,2)])
cdata_Cyto23_log<-cbind(cdata_Cyto23[,c(1,2)], cdata_Cyto23_log)
cdata_Cyto23_clean<-cdata_Cyto23_log
rownames(cdata_Cyto23_clean)<-cdata_Cyto23_clean$"Mouse ID"

#put them together.
all_dat<-cbind(cdata_Cyto23_clean, cdata_Cyto31_clean[cdata_Cyto23$"Mouse ID",-c(1:2)], cdata_Th17_clean[cdata_Cyto23$"Mouse ID",-c(1:2)])


###take care the names first, getting ccl names when they are available.
#	get rid of the duplicated whenever they are.
# take care the name duplication, and also need to get ccl name
cnames<-colnames(all_dat)[-c(1,2)]
cnames.ccl<-sub(x=cnames,"^[a-zA-Z0-9\\-]+","")
cnames.ccl<-sub(x=cnames.ccl,"^\\([a-zA-Z0-9\\-]+\\)","")
cnames.ccl<-sub(x=cnames.ccl,"^\\/","")
cnames.nccl<-sub(x=cnames, "\\/[a-zA-Z0-9\\-]+$","")

for(i in 1:length(cnames.ccl))
{
	if(cnames.ccl[i] !="")
	{
		#replace with ccl names
		cnames[cnames.nccl==cnames.nccl[i]]=cnames.ccl[i]
	}
}
colnames(all_dat)[-c(1:2)]<-cnames

index.dup<-duplicated(colnames(all_dat))
#to make things easy we get rid of the duplicated ones
all_dat.original<-all_dat

all_dat<-all_dat[,!index.dup]

#now let's plot them together into the heatmap
#all_dat$Group<-factor(cdata_Cyto31$Group, levels=c("1","2","5","6","8","9","2b","8b","9b"))
mouse_group_all<-as.data.frame(all_dat[,1])

rownames(mouse_group_all)<-all_dat$"Mouse ID"
colnames(mouse_group_all)<-"Group"
all_mat<-t(all_dat[,-c(1:2)])
colnames(all_mat)<-all_dat$"Mouse ID"
#gaps_col<-aggregate(cdata_Cyto31$"Mouse ID", by=list(cdata_Cyto31$Group), length)
ReadMe="all_mat is the clean data (removing duplicated cytokines. all_dat.original contains all cytokines data (including duplicated and log transformed)"
save(file=here("Data","cdata_all_MC38.RData"), all_dat.original, all_dat, ReadMe 
			)


c_all.scaled<-pheatmap(all_mat,scale="row",  cluster_cols = F, cluster_rows=T,
			show_colnames = F, #border_color = NA,
			annotation_col=mouse_group_all
			#, gaps_col=c(48)
			)
c_all<-pheatmap(all_mat,scale="none",  cluster_cols = F, cluster_rows=T,
			show_colnames = F, #border_color = NA,
			annotation_col=mouse_group_all
			#, gaps_col=c(48)
			)
png(file="cytokines_all.png", width=600, height=900)
c_all
dev.off()

png(file="cytokines_all_scaled.png", width=600, height=900)
c_all.scaled
dev.off()

#now do group heatmap.
# 
temp<-all_dat
#loop through columns to do averaging.
g.means<-data.frame()
for( i in 3:dim(temp)[2])
{
	if(i==3)
	{
		g.means=aggregate(temp[,i], by=list(temp$"Group"), mean)
		rownames(g.means)<-g.means[,1]
		colnames(g.means)[2]<-colnames(temp)[i]
	}
	else
	{	a<-aggregate(temp[,i], by=list(temp$"Group"), mean)
		rownames(a)<-a[,1]
		colnames(a)[2]<-colnames(temp)[i]
		g.means<-cbind(g.means,a[rownames(g.means),2])
	}
}
colnames(g.means)[-1]<-colnames(temp)[-c(1:2)]
colnames(g.means)[1]<-"Group"

#
g_mat<-t(g.means[,-c(1)])
colnames(g_mat)<-g.means$"Group"
#gaps_col<-aggregate(cdata_Cyto31$"Mouse ID", by=list(cdata_Cyto31$Group), length)
x<-rownames(g_mat)
			index.tnf<-which(x=="TNF-a")
			index.ifn<-which(x=="IFN-g")
			x[index.tnf]<-expression(paste("TNF",alpha,""))
			x[index.ifn]<-expression(paste("IFN",gamma,""))
c_g.scaled<-pheatmap(g_mat,scale="row",  cluster_cols = F, cluster_rows=T,
			show_colnames = T, 
			labels_row=x,#border_color = NA,
			#annotation_col=mouse_group_all
			#, gaps_col=c(48)
			)
			
			
c_g<-pheatmap(g_mat,scale="none",  cluster_cols = F, cluster_rows=T,
			show_colnames = T, #border_color = NA,
			labels_row=x,
			#annotation_col=mouse_group_all
			#, gaps_col=c(48)
			)
			
png(file="cytokines_group_scaled.png", width=400, height=900)
c_g.scaled
dev.off()
png(file="cytokines_group.png", width=400, height=900)
c_g
dev.off()

pdf(file="cytokines_group_scaled.pdf", width=4, height=9)
c_g.scaled
dev.off()

#now we need to modify the names to make them the entrez symbols
cname<-colnames(all_dat)

symbol<-cname
symbol[symbol=="G-CSF"]<-"CSF3"
symbol[symbol=="GM-CSF"]<-"CSF2"
symbol[symbol=="IL-12(p40)"]<-"IL12B"
symbol[symbol=="IL-12(p70)"]<-"IL12A"
symbol[symbol=="TNF-a"]<-"TNF"
symbol[symbol=="CD40L"]<-"CD40LG"
symbol[symbol=="IFN-g"]<-"IFNG"
#symbol[symbol=="IL-1a"]<-"IL1A"
#symbol[symbol=="IL-1b"]<-"IL1B"
colnames(all_dat)<-symbol

#now let's do anova and then fsd
#
anova.table<-data.frame()
temp<-all_dat

for( i in 3:dim(temp)[2])
{
		#doing anova
		x<-temp[,i]
		grp<-temp$Group
		temp.lm<-lm(x~grp)
		#follow up
		pd<-TukeyHSD(aov(temp.lm))
		
	if(i==3)
	{
		
		anova.table<-data.frame(cytokine=colnames(temp)[i], 
				"anova.F_Value"=anova(temp.lm)[1,"F value"],"anova.Pr"=anova(temp.lm)[1,"Pr(>F)"])
		anova.table<-cbind(anova.table,t(pd$grp[,4] ))
	}
	else
	{	
		a<-data.frame(cytokine=colnames(temp)[i], 
				"anova.F_Value"=anova(temp.lm)[1,"F value"],"anova.Pr"=anova(temp.lm)[1,"Pr(>F)"])
		a<-cbind(a,t(pd$grp[,4] ))
		anova.table<-rbind(anova.table, a)
	}
	
}
anova.table$p.adj<-p.adjust(anova.table[,3])

anova.table<-anova.table[,c(1:3,dim(anova.table)[2], 4:(dim(anova.table)[2]-1))]
anova.table.NA<-anova.table
#remove not significant ones
anova.table.NA[anova.table.NA$p.adj>0.05, c(5: (dim(anova.table.NA)[2])) ]<-NA
colnames(anova.table.NA)<-sub(x=colnames(anova.table.NA), "-","v",fixed=T)
write.csv(anova.table.NA,file="cytoke_analysis_anova_Tukey.cvs")
#write_xlsx(anova.table, file="cytokine_anova_Tukey.xlsx", col_names=T)

#now start doing the gsea
library(fgsea)
library(tibble)
#first read the data base gmt
pathways<-gmtPathways("msigdb_v2022.1.Hs_files_to_download_locally\\msigdb_v2022.1.Hs_GMTs\\msigdb.v2022.1.Hs.symbols.gmt")
# prepare the rank list file

anova.table$cytokine<-sub(x=anova.table$cytokine, "-","", fixed=T)
anova.table$cytokine<-toupper(anova.table$cytokine)
rownames(anova.table)<-anova.table$cytokine
library(tidyr)
#get name done first
colnames(all_dat)<-sub(x=colnames(all_dat), pattern="-", "", fixed=T)
colnames(all_dat)[-c(1,2)]<-toupper(colnames(all_dat)[-c(1,2)])
#reshape the data table to long format
all_dat_long<- all_dat %>% pivot_longer(where(is.numeric), names_to="Cytokine", values_to="Conc")
# first group 2 vs 1.
groupMeans<-aggregate(all_dat_long$Conc, by=list(group=all_dat_long$Group,Cytokine=all_dat_long$Cytokine),FUN=mean)

g1.means<-groupMeans[groupMeans$group=="1",]
rownames(g1.means)<-g1.means$Cytokine

g2.means<-groupMeans[groupMeans$group=="2",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"2-1")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways, rnk, minSize=5) %>% as_tibble() %>% arrange(pval)


res_2v1<-res
			  
#6 v 1
g2.means<-groupMeans[groupMeans$group=="6",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"6-1")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways, rnk, minSize = 5) %>% as_tibble() %>% arrange(pval)
res_6v1<-res		  


#9 v 1
g2.means<-groupMeans[groupMeans$group=="9",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"9-1")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways, rnk, minSize = 5) %>% as_tibble() %>% arrange(pval)
res_9v1<-res	

#8 v 1
g2.means<-groupMeans[groupMeans$group=="8",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"8-1")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways, rnk, minSize = 5) %>% as_tibble() %>% arrange(pval)
res_8v1<-res		  
	  
#5 v 1
g2.means<-groupMeans[groupMeans$group=="5",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"5-1")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways, rnk, minSize = 5) %>% as_tibble() %>% arrange(pval)
res_5v1<-res
 library(readr)

write_csv(res_2v1,file="cytoke_gsea_2v1_mSigDB.cvs")
write_csv(res_5v1,file="cytoke_gsea_5v1_mSigDB.cvs" )
write_csv(res_6v1,file="cytoke_gsea_6v1_mSigDB.cvs" )
write_csv(res_8v1,file="cytoke_gsea_8v1_mSigDB.cvs" )
write_csv(res_9v1,file="cytoke_gsea_9v1_mSigDB.cvs" )

#plot msigdb results
res_2v1$Group="2v1"
res_5v1$Group="5v1"
res_6v1$Group="6v1"
res_8v1$Group="8v1"
res_9v1$Group="9v1"

all_msigdb<-rbind(res_2v1, res_5v1, res_6v1, res_8v1, res_9v1)
num<-10
	temp_2v1<-res_2v1[order(res_2v1$NES, decreasing=T)[c(1:num, (dim(res_2v1)[1]-num):dim(res_2v1)[1])],]
	temp_5v1<-res_5v1[order(res_5v1$NES, decreasing=T)[c(1:num, (dim(res_5v1)[1]-num):dim(res_5v1)[1])],]
	temp_6v1<-res_6v1[order(res_6v1$NES, decreasing=T)[c(1:num, (dim(res_6v1)[1]-num):dim(res_6v1)[1])],]
	temp_8v1<-res_8v1[order(res_8v1$NES, decreasing=T)[c(1:num, (dim(res_8v1)[1]-num):dim(res_8v1)[1])],]
	temp_9v1<-res_9v1[order(res_9v1$NES, decreasing=T)[c(1:num, (dim(res_9v1)[1]-num):dim(res_9v1)[1])],]
	all_msigdb<-rbind(temp_2v1, temp_5v1, temp_6v1, temp_8v1, temp_9v1)

	#temp<-res_2v1[order(res_2v1$NES, decreasing=F)[1:20],]
all_msigDB_plot<-ggplot(all_msigdb, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA (mSigDB)") +	   
  facet_wrap(~Group, scales = "free_y", nrow=5)
  
  png("Cytokine_gsea_mSigDB.png", width=800, height=1200)
	all_msigDB_plot
	dev.off()

###### combo vs pd1 (2)####
#        8v2
#########################
g1.means<-groupMeans[groupMeans$group=="2",]
rownames(g1.means)<-g1.means$Cytokine

g2.means<-groupMeans[groupMeans$group=="8",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"8-2")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways, rnk, minSize=5) %>% as_tibble() %>% arrange(pval)
#plot
res_8v2_msigDB_plot<-ggplot(res, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA") + 
  theme_minimal()
res_8v2_msigDB<-res
##############################
g1.means<-groupMeans[groupMeans$group=="2",]
rownames(g1.means)<-g1.means$Cytokine

g2.means<-groupMeans[groupMeans$group=="9",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"9-2")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways, rnk, minSize=5) %>% as_tibble() %>% arrange(pval)
#plot
res_9v2_msigDB_plot<-ggplot(res, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA") + 
  theme_minimal()
res_9v2_msigDB<-res
#putting them together 
res_8v2_msigDB$Group<-"8v2"
res_9v2_msigDB$Group<-"9v2"

temp_8v2<-res_8v2_msigDB[order(res_8v2_msigDB$NES, decreasing=T)[c(1:num, (dim(res_8v2_msigDB)[1]-num):dim(res_8v2_msigDB)[1])],]
	temp_9v2<-res_9v2_msigDB[order(res_9v2_msigDB$NES, decreasing=T)[c(1:num, (dim(res_9v2_msigDB)[1]-num):dim(res_9v2_msigDB)[1])],]
	all_msigdb_combo<-rbind(temp_8v2, temp_9v2)

	#temp<-res_2v1[order(res_2v1$NES, decreasing=F)[1:20],]
all_msigDB_combo_plot<-ggplot(all_msigdb_combo, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA (mSigDB)") +	   
  facet_wrap(~Group, scales = "free_y", nrow=2)
  
  png("Cytokine_gsea_mSigDB_8or9vs2.png", width=800, height=600)
	all_msigDB_combo_plot
	dev.off()
	
	write_csv(res_8v2_msigDB,file="cytoke_gsea_8v2_msigDB.xls")
	write_csv(res_9v2_msigDB,file="cytoke_gsea_9v2_msigDB.cvs" )
####################
### doing the hallmark pathways###
###################################

  
#  all_msigDB_plot<-ggplot(all_msigdb[all_msigdb$pval<=0.05,], aes(reorder(pathway, NES), NES)) +
#  geom_col(aes(fill=pval<0.05)) +
#  coord_flip() +
#  labs(x="Pathway", y="Normalized Enrichment Score",
#       title="Hallmark pathways NES from GSEA") +
#	facet_wrap(~Group, scales = "free_y")

####doing the hallmark pathways.
pathways.hallmark<-gmtPathways("msigdb_v2022.1.Hs_files_to_download_locally\\msigdb_v2022.1.Hs_GMTs\\h.all.v2022.1.Hs.symbols.gmt")
# prepare the rank list file

#2v1
g2.means<-groupMeans[groupMeans$group=="2",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"2-1")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways.hallmark, rnk, minSize=5) %>% as_tibble() %>% arrange(pval)
#plot
res_2v1_plot<-ggplot(res, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
res_2v1_hallmark<-res

#5 v 1
g2.means<-groupMeans[groupMeans$group=="5",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"5-1")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways.hallmark, rnk, minSize = 5) %>% as_tibble() %>% arrange(pval)

#plot
res_5v1_plot<-ggplot(res, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
res_5v1_hallmark<-res

#6 v 1
g2.means<-groupMeans[groupMeans$group=="6",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"6-1")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways.hallmark, rnk, minSize = 5) %>% as_tibble() %>% arrange(pval)
#plot
res_6v1_plot<-ggplot(res, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
res_6v1_hallmark<-res	

#8 v 1
g2.means<-groupMeans[groupMeans$group=="8",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"8-1")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways.hallmark, rnk, minSize = 5) %>% as_tibble() %>% arrange(pval)
res_8v1_plot<-ggplot(res, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
res_8v1_hallmark<-res

#9 v 1
g2.means<-groupMeans[groupMeans$group=="9",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"9-1")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways.hallmark, rnk, minSize = 5) %>% as_tibble() %>% arrange(pval)

res_9v1_plot<-ggplot(res, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
res_9v1_hallmark<-res	

###### combo vs pd1 (2)####
#        8v2
#########################
g1.means<-groupMeans[groupMeans$group=="2",]
rownames(g1.means)<-g1.means$Cytokine

g2.means<-groupMeans[groupMeans$group=="8",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"8-2")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways.hallmark, rnk, minSize=5) %>% as_tibble() %>% arrange(pval)
#plot
res_8v2_plot<-ggplot(res, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
res_8v2_hallmark<-res

g1.means<-groupMeans[groupMeans$group=="2",]
rownames(g1.means)<-g1.means$Cytokine

g2.means<-groupMeans[groupMeans$group=="9",]
rownames(g2.means)<-g2.means$Cytokine

fc2<-g1.means[rownames(anova.table),"x"]-g2.means[rownames(anova.table),"x"]
anova.table$rank<- -log10(anova.table$"9-2")*(-1*fc2)
rnk<-anova.table[,c("cytokine","rank")]
rnk<-deframe(rnk)

res<- fgsea(pathways=pathways.hallmark, rnk, minSize=5) %>% as_tibble() %>% arrange(pval)
#plot
res_9v2_plot<-ggplot(res, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
res_9v2_hallmark<-res

plotEnrichment(pathways.hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], rnk)

#put them together of all vs 1. to draw a figure.
res_2v1_hallmark$Group<-"2 vs 1"
all_hallmark_gsea<-res_2v1_hallmark

res_5v1_hallmark$Group<-"5 vs 1"
all_hallmark_gsea<-rbind(all_hallmark_gsea, res_5v1_hallmark)
res_6v1_hallmark$Group<-"6 vs 1"
all_hallmark_gsea<-rbind(all_hallmark_gsea, res_6v1_hallmark)

res_8v1_hallmark$Group<-"8 vs 1"
all_hallmark_gsea<-rbind(all_hallmark_gsea, res_8v1_hallmark)

res_9v1_hallmark$Group<-"9 vs 1"
all_hallmark_gsea<-rbind(all_hallmark_gsea, res_9v1_hallmark)

all_plot<-ggplot(all_hallmark_gsea, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  #theme_minimal()+
  theme(legend.position = c(0.8, 0.2)) +
  facet_wrap(~Group)
png(file="gsea_hallmark_result.png", width=800, height=500)
all_plot
dev.off()  
#write to the file

write_csv(res_2v1_hallmark,file="cytoke_gsea_2v1_hallmark.xls")
write_csv(res_5v1_hallmark,file="cytoke_gsea_5v1_hallmark.cvs" )
write_csv(res_6v1_hallmark,file="cytoke_gsea_6v1_hallmark.cvs" )
write_csv(res_8v1_hallmark,file="cytoke_gsea_8v1_hallmark.cvs" )
write_csv(res_9v1_hallmark,file="cytoke_gsea_9v1_hallmark.cvs" )

########put together 8/9 v 2
res_8v2_hallmark$Group<-"8v2"
res_9v2_hallmark$Group<-"9v2"

combo_hallmarm_gsea<-res_8v2_hallmark
combo_hallmarm_gsea<-rbind(combo_hallmarm_gsea, res_9v2_hallmark)
combo_all_plot<-ggplot(combo_hallmarm_gsea, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  #theme_minimal()+
  theme(legend.position = c(0.8, 0.2)) +
  facet_wrap(~Group, nrow=2)
  png(file="cytokine_gsea_hallmark_8or9v2.png", width=600, height=500)
	combo_all_plot 
	dev.off()
## write to the file
write_csv(res_8v2_hallmark,file="cytoke_gsea_8v2_hallmark.xls")
write_csv(res_9v2_hallmark,file="cytoke_gsea_9v2_hallmark.cvs" )

save(file="dataset_gsea.RData", all_hallmark_gsea, all_msigdb, all_msigdb_combo, combo_hallmarm_gsea)
#save(

####start doing for the complete msigdb
	

	  
#################left over don't run

topPathwaysUp <- res[res$ES > 0,][head(order(res$pval), n=10), "pathway"]
topPathwaysDown <- res[res$ES < 0,][head(order(res$pval), n=10), "pathway"]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], rnk, res, 
              gseaParam = 0.5)
