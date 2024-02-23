# R code to do further analysis to see patterns
#=== Start this module 
#		so far only work on MC38
#			1/19/2023
####################################
#		- PCA
#      - MFA
#  for now. 12/15/2022
# 
#
#-----------add more more modes of data
#	1. tumor volume
#	2. tumor cell infiltration
#		1/17/2023
#     
#-----------------------------------

# Using the data generated from cytokine_process.R in both experiments MC38 (cytokine) and MB49 (cytokine2) 
#

### load libraries
library(FactoMineR)
library(factoextra)
library(ggpubr)
library(ggalt)
library(tidyverse)
library(readxl)
library(dplyr)
library(magrittr)
library(usethis)
library(here)
library(corrplot)
library(randomForest)
library(caret)

####now read the data, whenever possible, let's save/read data from shared folder
#setwd("/mnt/sb_ub1/share/Feng/VISTA/Cytokines/") #linux machine
#setwd("E:\\feng\\personal_record\\cv\\jobSearch2\\SenseiBio\\Cytokines")
#both data sets were saved by cytokine_process.R module in separate experiments.
load(here("Data","cdata_all_MC38.RData")) 
mc38<-all_dat.original
mc38$group<-rep(c("Ctrl", "aPD1","aSNS101_10mg/Kg", "aSNS101_30mg/Kg", "aPD1+aSNS101_10mg/Kg","aPD1+aSNS101_30mg/Kg"), 
											times=rep(8,times=6))
mc38$group<-factor(mc38$group, levels=c("Ctrl", "aPD1","aSNS101_10mg/Kg", "aSNS101_30mg/Kg", "aPD1+aSNS101_10mg/Kg","aPD1+aSNS101_30mg/Kg"))
											
mc38<-mc38[,-1]
mc38$tumor<-"MC38"


write_csv(file=here("Data","cdata_all_MC38_log.csv"),mc38[,c(66,1:65)])

#read the tumore and cell infiltration data.
load(here("Data","data_tumor_cellInfil.RData"))
volume.tumor<-as.data.frame(volume.tumor)
rownames(volume.tumor)<-volume.tumor$Mouse_ID

temp<-cbind(volume.tumor,mc38[rownames(volume.tumor),"group"])
colnames(temp)[8]<-"group"
temp <-temp[!is.na(temp$group),]

write_csv(file=here("Data","tumorVolume_all_MC38.csv"),
		temp[,c(7,8,1:5)])


	####mouse number errors!!!! Don't run for now.
		cells.infiltration<-as.data.frame(cells.infiltration)
		rownames(cells.infiltration)<-cells.infiltration$Mouse_ID
	#####################################

#now put the data together

all.data<-cbind(mc38, volume.tumor[rownames(mc38),])
		#cells.infiltration[rownames(mc38),])
all.data<-all.data[,-c(1,68:71,73:74)]

#plot the distribution of tumor size
x<-lm(Volume~group,data=all.data)
anova(x)
x.av<-aov(x)
TukeyHSD(x.av)

sum.mean<-aggregate(all.data$Volume, by=list(all.data$group), FUN=mean)
names(sum.mean)<-c("Group", "Volume")
sum.se<-aggregate(all.data$Volume, by=list(all.data$group), 
	FUN=function(x){sd(x)/sqrt(length(x))})
sum.data<-cbind(sum.mean, sum.se)
names(sum.data)<-c("Group", "M","Group2","SE")

#all.data$group<-factor(all.data, levels=c(""))
#plot the group-wise tumor volume changes
dp<-ggplot(all.data, aes(x=group, y=Volume, fill=group))+
	#geom_bar(sum.data,aes(x=Group, y=Volume, fill=Group),stat="identity")+
	# geom_errorbar(sum.data,aes(ymin=M-SE, ymax=M+SE), width=.2) 
	geom_boxplot()+
	 geom_dotplot(
	 	binaxis = "y", stackdir = "center",binwidth=50)

bp<-ggplot(sum.data, aes(x=Group, y=M, fill=Group))+
	geom_bar(stat="identity")+
	ylab(expression(paste("Tumor Volume (",um^3,")")))+
	geom_errorbar(aes(ymin=M-SE, ymax=M+SE), width=.2) +
	theme(legend.position="none",
			axis.text.x=element_blank(),
			axis.title.x=element_blank(),
			axis.title.y=element_text(size=18))
png("tumorVol_barplot_nolabel.png", width=750,height=500)
bp
dev.off()

## first we try pca (not MFA yet)
#try correlation
mat<-all.data[,-c(65,66)]#-c("Group","tumor")]
mat <- mat %>% 
	select(-CCL2.1)

mat.cor<-(cor(mat)) 
#the last row is the correlation to volume
#we sort to get the order
volume.cor<-mat.cor[c("Volume"),] 
index<-order(volume.cor)

#remove the CCL2.1



num_cytokines<-8
mat.sort<-mat[,index]

mat.sort <- mat.sort %>% 
	rename("CCL3 "="CCL3.1", " IL-6"="IL-6.1", 
				"Tumor Vol"="Volume")

#ind<-grep(x=colnames(mat.sort),pattern="IL-1a") 
#colnames(mat.sort)[ind]<-expression(paste("IL-1",alpha,""))
mat.cor<-cor(mat.sort[,
		c(1:num_cytokines,(dim(mat.cor)[2]-num_cytokines):dim(mat.cor)[2])])
write.csv(file=here("Data","correlationMatrix_all.csv"),
			as.data.frame(cor(mat.sort)))
colnames(mat.cor)[5]=":IL-1~alpha"
colnames(mat.cor)[6]=":TNF-alpha"
colnames(mat.cor)[13]=":IL-1~beta"

rownames(mat.cor)[5]=":IL-1~alpha"
rownames(mat.cor)[6]=":TNF-alpha"
rownames(mat.cor)[13]=":IL-1~beta"

#now let's do correlation matrix
png(file="correlation.png", width=600, height=600)
corrplot(mat.cor, type="upper", #order="hclust", 
         tl.col="black", tl.srt=45,tl.cex=1.3)
dev.off()

pdf(file="correlation.pdf", width=6, height=6)
corrplot(mat.cor, type="upper", #order="hclust", 
         tl.col="black", tl.srt=45,tl.cex=1.1)
dev.off()
png(file="volume_comparison_bar.png", width=700, height=500)
bp
dev.off()
png(file="volume_comparison_box.png", width=700, height=500)
dp
dev.off()
#do pca first 
mat<-all.data[,-c(65,66)]#-c("Group","tumor")]
results<-prcomp(mat, scale.=T, center=T)
#results<-PCA(mat, scale.=T, center=T)
#prcomp(mb49[,-c(65,66],scale.=T, center=F)
var_explained <-results$sdev^2 / sum(results$sdev^2)

#create scree plot
num<-10
plot(c(1:num),var_explained[1:num], xlab = "principal component", ylab = "Proportion of Variance Explained",
		ylim = c(0, 1), type = "b", main = "Scree Plot")
		
plot(results$x[,1], results$x[,2], col=as.integer(factor(all.data$tumor)))

	
plot(results$x[,1], results$x[,2], col=as.integer(factor(all.data$group)), pch=as.integer(factor(all.data$group)))

#doiing PCA without supplementary
#doing separate analysis
#mc38first
res.mc38<-prcomp(mc38[,-c(1,66,67)],scale.=T, center=T)

fig_eig.mc38<-fviz_eig(res.mc38,addlabels = TRUE, ylim = c(0, 50))
	
	png(file="scree_mc38.png",width=550, height=400)
	fig_eig.mc38
	dev.off()
fviz_pca_ind(res.mc38, axes=c(1,2),
             habillage = as.factor(mc38$group), # color by groups 
			 geom=c("point"), #choice=c("group"),
             #palette = c("#00AFBB", "#E7B800","red","blue","pink","),
             addEllipses = T, #ellipse.type = "confidence", 
			 ellipse.level=0.95, pointsize=4,
             repel = TRUE # Avoid text overlapping
             ) 

res.mc38$contrib<-(res.mc38$sdev^2)/sum(res.mc38$sdev^2)
			 
fig_mc38_pc12<-ggplot(as.data.frame(res.mc38$x), aes(x=PC1, y=PC2, 
			group=as.factor(mc38$group), colour=as.factor(mc38$group)))+
	geom_point(aes(size=all.data$Volume))+
	geom_encircle()+ xlim(c(-16,15))+ylim(c(-8,8))+
	geom_vline(xintercept=0,color="grey",linetype=2,linewidth=1.6)+
	geom_hline(yintercept=0,color="grey",linetype=2,linewidth=1.6)+
	xlab(paste0("PC1 (",format(res.mc38$contrib[1]*100, nsmall=1, digit=2),"%)"))+
	ylab(paste0("PC2 (",format(res.mc38$contrib[2]*100, nsmall=1, digit=2),"%)"))+
	guides(color = guide_legend(title = ""),
		 size = "none"#guide_legend(title = "Tumor Vol")
		 )+
	#scale_size_continuous(name="Tumor Vol", range = c(2,15))+
	scale_colour_discrete(labels=
		c("Control",expression(alpha*"PD1"),
			expression(alpha*"SNS101+"), expression(alpha*"SNS101++"),
			expression(alpha*"PD1&SNS101+"),
			expression(alpha*"PD1&SNS101++")))+
	theme(legend.position=c(0.15,0.11),legend.text = element_text(size=20))
png(file="PCA_mc38_tumorVolumeSize.png", width=750,height=750)
fig_mc38_pc12
dev.off()
 
#	fig_mc38_pc34<-ggplot(as.data.frame(res.mc38$x), aes(x=PC3, y=PC4, 
#				group=as.factor(mc38$group), colour=as.factor(mc38$group)))+
#		geom_point(size=5)+geom_encircle()+xlim(c(-16,15))+ylim(c(-8,8))+
#		geom_vline(xintercept=0,color="grey",linetype=2,size=1.6)+
#		geom_hline(yintercept=0,color="grey",linetype=2,size=1.6)+
#		xlab(paste0("PC3 (",format(res.mc38$contrib[3]*100, nsmall=1, digit=2),"%)"))+
#		ylab(paste0("PC4 (",format(res.mc38$contrib[4]*100, nsmall=1, digit=2),"%)"))+
#		guides(color = guide_legend(title = "Groups"))+
#		theme()	

#png(file="fig_ind_pc12.mc38.png", width=800, height=700)
#fig_mc38_pc12 #, fig_ind_pc34, nrow=1, ncol=2)
#dev.off()
#	png(file="fig_ind_pc34.mc38.png", width=800, height=700)
#	fig_mc38_pc34 #, fig_ind_pc34, nrow=1, ncol=2)
#	dev.off()
	
#	fig_var_contrib_pc1<-fviz_contrib(res.mc38, choice = "var", axes = 1, top = 40,
#	             #palette = "jco"
#				 )+theme(legend.position="none",axis.text=element_text(size=12))
#	fig_var_contrib_pc2<-fviz_contrib(res.mc38, choice = "var", axes = 2, top = 40,
#	             #palette = "jco"
#				 )+theme(legend.position="none",axis.text=element_text(size=12))
#
#	fig_var_contrib_pc3<-fviz_contrib(res.mc38, choice = "var", axes = 3, top = 40,
#	             #palette = "jco"
#				 )+theme(legend.position="none",axis.text=element_text(size=12))
#	png(file="fig_var_contrib.mc38.png", width=900, height=900)
#	ggarrange(fig_var_contrib_pc1, fig_var_contrib_pc2,fig_var_contrib_pc3, nrow=3, ncol=1)
#	dev.off()
			 
	#individual
#	fig_ind_contrib_pc1<-fviz_contrib(res.mc38, choice = "ind", axes = 1, top = 40,
#	             #palette = "jco"
#				 )+theme(legend.position="none",axis.text=element_text(size=12))
#	fig_ind_contrib_pc2<-fviz_contrib(res.mc38, choice = "ind", axes = 2, top = 40,
#	             #palette = "jco"
#				 )+theme(legend.position="none",axis.text=element_text(size=12))
#
#	fig_ind_contrib_pc3<-fviz_contrib(res.mc38, choice = "ind", axes = 3, top = 40,
#	             #palette = "jco"
#				 )+theme(legend.position="none",axis.text=element_text(size=12))
#	png(file="fig_ind_contrib.mc38.png", width=900, height=900)
#	ggarrange(fig_ind_contrib_pc1, fig_ind_contrib_pc2,fig_ind_contrib_pc3, nrow=3, ncol=1)
#	dev.off()

#doing PCA with supplementary
pca.2<-PCA(mat, scale.unit = TRUE, ncp = 15, 
		quanti.sup =65, graph = FALSE)
png(file="pca_variable_sup.png",width=650, height=650)
plot(pca.2, choix = "var")
dev.off()

#make a new df with pcs and volume

dfm<-as.data.frame(res.mc38$x[,1:20])

dfm<-cbind(dfm, mat[rownames(dfm),"Volume"])
colnames(dfm)[dim(dfm)[2]]<-"Volume"
#dfm<-as.data.frame(dfm)
lm.pcs<-lm(Volume~PC1+PC2+PC3+PC4+#PC5+PC6+PC7+PC8+
		PC9+#PC10+PC11+PC12+PC13+PC14+PC15+
		PC16+PC17#+PC18+PC19+PC20
		, data=dfm)
summary(lm.pcs)

lm.pcs.all<-lm(Volume~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+
		PC9+PC10+PC11+PC12+PC13+PC14+PC15+
		PC16+PC17+PC18+PC19+PC20
		, data=dfm)
summary(lm.pcs.all)


##doing variable selection using machine learning
#using random forrest algorithm
control <- rfeControl(functions = rfFuncs, # random forest
                      method = "repeatedcv", # repeated cv
                      repeats = 5, # number of repeats
                      number = 10) # number of folds

#prepare the input, predicting vs independent varaibles             
x <- mat %>%
  select(-Volume) %>%
  as.data.frame()
x<-scale(x, center=T, scale=T)

# Target variable
y <- mat$Volume

# Training: 80%; Test: 20%
set.seed(2021)
inTrain <- createDataPartition(y, p = .99, list = FALSE)[,1]

x_train <- x[ inTrain, ]
x_test  <- x[-inTrain, ]

y_train <- y[ inTrain]
y_test  <- y[-inTrain]

result_rfe1 <- rfe(x = x_train, 
                   y = y_train, 
                   sizes = c(1:40),#size is number of features in the model.
                   rfeControl = control)

# Print the results
result_rfe1

# Print the selected features
predictors(result_rfe1)

# Print the results visually
#x1<-ggplot(data = result_rfe1, metric = "Accuracy") + theme_bw()
#x2<-ggplot(data = result_rfe1, metric = "Kappa") + theme_bw()

x3<-ggplot(data = result_rfe1, metric = "RMSE") + theme_bw()

#ggarrange(x1,x2,x3, nrow=2)


#importance
varimp_data <- data.frame(feature = row.names(varImp(result_rfe1))[1:10],
                          importance = varImp(result_rfe1)[1:10, 1])
varimp_data$feature[7]<-"CCL2"

varimp_data$feature[10]<-"IL-2"

varimp_data <- varimp_data %>% filter(!is.element (feature,  c("IL-27", "CCL11","IL-2")))

imp<-ggplot(data = varimp_data, 
       aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
  geom_bar(stat="identity") + labs(x = "Cytokines", y = "Variable Importance") + 
  geom_text(aes(label = round(importance, 2)), vjust=1.6, color="white", size=4) + 
  theme_bw() + theme(legend.position = "none",
  		axis.text.y = element_text(size = 15),
  		axis.text.x =element_text(size=15),
  		axis.title.x =element_text(size=15),
  		axis.title.y =element_text(size=15),
  		)   

  	

png(file="rfe_RMSE_imp.png",width=750, height=1000)
ggarrange(x3,imp, nrow=2)
dev.off()
pdf(file="rfe_RMSE_imp.pdf", width=7.5, height=10)
ggarrange(x3,imp, nrow=2)
dev.off()
png(file="rfe_RMSE_imp_bar.png",width=750, height=550)
imp
dev.off()

#lastly, do linear regression
lm.rfe<-lm(Volume~CCL7+mat$"IL-3"+CXCL12+CCL24+CCL12+CXCL10+CCL2.1,data=mat)





	#################DON"T RUN##########
	##   left over from previous analysis
	##############3

###doing MFA with tumor volume as supplementary for now
group<-rep(1, dim(all.data)[2])
type<-c(rep("s", 64),"n","n","s")
ind.sup<-c(65,66)
name.group<-names(all.data)

res.mfa <- MFA(all.data, 
               group = group, 
               type = type,#ind.sup=ind.sup,
               name.group = name.group,
               num.group.sup = c(65, 66), #what this is?
               graph = FALSE)

#contribution or eigen
fig.cont<-fviz_eig(res.mfa, addlabels = TRUE, ylim = c(0, 60))
pdf(file="scree.pdf",width=7, height=4)
fig.cont
dev.off()

png(file="scree_both.png",width=550, height=400)
fig.cont
dev.off()
#individual			   
fig_ind_pc12<-fviz_pca_ind(res.mfa, axes=c(1,2),
             habillage = as.factor(all.data$tumor), # color by groups 
			 geom=c("point"), #choice=c("group"),
             #palette = c("#00AFBB", "#E7B800","red","blue","pink","),
             addEllipses = T, #ellipse.type = "confidence", 
			 ellipse.level=0.95, pointsize=4,
             repel = TRUE # Avoid text overlapping
             ) 
fig_ind_pc34<-fviz_pca_ind(res.mfa, axes=c(3,4),
             habillage = as.factor(all.data$tumor), # color by groups 
			 geom=c("point"), #choice=c("group"),
             #palette = c("#00AFBB", "#E7B800","red","blue","pink","),
             addEllipses = T, #ellipse.type = "confidence", 
			 ellipse.level=0.95, pointsize=4, alpha.var=0,
             repel = TRUE # Avoid text overlapping
             ) +xlim(c(-6,6.5))+ylim(c(-6,6.5))
png(file="fig_ind_pc12.png", width=800, height=700)
fig_ind_pc12 #, fig_ind_pc34, nrow=1, ncol=2)
dev.off()
png(file="fig_ind_pc34.png", width=800, height=700)
fig_ind_pc34 #, fig_ind_pc34, nrow=1, ncol=2)
dev.off()

#variable (cytokines)
fig_var<-fviz_mfa_var(res.mfa, choice=c("quanti.var"), #palette = "jco", 
             col.var.sup = "violet", repel = TRUE, select.var=list(contrib=10),
			 geom = c("point", "text"), legend = "bottom")

# Contributions to dimension 1
fig_var_contrib<-fviz_contrib(res.mfa, choice = "quanti.var", axes = 1, top = 40,
             #palette = "jco"
			 )+theme(legend.position="none")
png(file="fig_var_both.png", width=1000, height=700)
		ggarrange(fig_var, fig_var_contrib, nrow=2, ncol=1)
dev.off()
		




##### start plotting cytokine data (heatmaps) to back up the findings in MFA for **MC38** data.
#get the list of cytokines contributing to the PC1 changes.

contrib.pc1.mc38<-res.mc38$rotation[,"PC1"]^2
m.contrib<-mean(contrib.pc1.mc38)  #this is the 1/numb of cytokines
std.contrib<-sd(contrib.pc1.mc38)

#first get cytokines that are contributing above the average
list.cytokines.ba.mc38<-names(contrib.pc1.mc38[contrib.pc1.mc38>m.contrib])
#draw heatmaps using these cytokines.
dat.ba.mc38<-mc38[, list.cytokines.ba.mc38]
dat.ba.mc38<-cbind(mc38[,c("Mouse ID", "group", "tumor")], dat.ba.mc38)
temp<-dat.ba.mc38
#loop through columns to do averaging.
g.means<-data.frame()
for( i in 4:dim(temp)[2])
{
	if(i==4)
	{
		g.means=aggregate(temp[,i], by=list(temp$"group"), mean)
		rownames(g.means)<-g.means[,1]
		colnames(g.means)[2]<-colnames(temp)[i]
	}
	else
	{	a<-aggregate(temp[,i], by=list(temp$"group"), mean)
		rownames(a)<-a[,1]
		colnames(a)[2]<-colnames(temp)[i]
		g.means<-cbind(g.means,a[rownames(g.means),2])
	}
}
colnames(g.means)[-1]<-colnames(temp)[-c(1:3)]
colnames(g.means)[1]<-"Group"

#need to read the statistics in order to label the changes.
#
stats.mc38<-read.table(file="cytoke_analysis_anova_Tukey.cvs", header=T, sep=",")

#
source("/mnt/sb_ub1/Feng/VISTA/Cytokines/functions.R")
stats.mc38$cytokine<-UnifyCytokineNames(stats.mc38$cytokine)
stats.mc38$cytokine<-sub(x=stats.mc38$cytokine, pattern=".2", ".1",fixed=T)

index<-which(duplicated(stats.mc38$cytokine))
stats.mc38$cytokine[index]<-paste0(stats.mc38$cytokine[index], ".1")
rownames(stats.mc38)<-stats.mc38$cytokine

g_mat<-t(g.means[,-c(1)])
colnames(g_mat)<-g.means$"Group"

#gaps_col<-aggregate(cdata_Cyto31$"Mouse ID", by=list(cdata_Cyto31$Group), length)
x<-rownames(g_mat)
index.nonsig<-c()
#now for each item in x, we need to consult the stat table to get significance.
for( i in 1:length(x))
{
	#check the stat table 
	#omnibus sig
	temp<-x[i]	
	if(stats.mc38[x[i],"p.adj"]<0.05)
	{
		#check the followup analysis for 
		if(stats.mc38[x[i],"X8v1"]<0.05 |stats.mc38[x[i],"X9v1"]<0.05)
		{
			temp<-paste0(temp, " *")
		}
		if(stats.mc38[x[i],"X8v2"]<0.05 |stats.mc38[x[i],"X9v2"]<0.05)
		{
			temp<-paste0(temp, " #")
		}
	}
	else{
		index.nonsig<-c(index.nonsig, i)
	}
	x[i]<-temp
}

x<-x[-index.nonsig]
g_mat<-g_mat[-index.nonsig,]

#			index.tnf<-which(x=="TNF-a")
#			index.ifn<-which(x=="IFN-g")
#			x[index.tnf]<-expression(paste("TNF",alpha,""))
#			x[index.ifn]<-expression(paste("IFN",gamma,""))


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
			
png("heatmap_cytokines_pc1Sig_mc38.png", width=400, height=900)
c_g.scaled
dev.off()



