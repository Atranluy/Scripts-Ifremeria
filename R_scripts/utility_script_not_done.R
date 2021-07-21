############" utility script from A.tran LU Y ##############
### Using mainly snprelate packages
###

library(SNPRelate)
library(gdsfmt)


########### calculating multi locus heterozygosity on individual #############


Multilocus_Het<-function(genofile,popmap){
genofile<-genofile
popmap<-popmap
sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
ID=NULL
for(i in popmap$V1){
  a<-which(i==sample_id)
  ID<-append(ID,a)
}
genomatrix<-snpgdsGetGeno(genofile,sample.id =sample_id[ID])
MLH_all_SNPs<-NULL
MLH_all_SNPs_noNA<-NULL
for(i in 1:nrow(genomatrix)){
  temp_ind<-genomatrix[i,]
  het_is<-which(temp_ind=="1")
  MLH_temp<-(length(het_is)/ncol(genomatrix))
  MLH_temp1<-(length(het_is)/sum(table(temp_ind,exclude =NA)))
  MLH_all_SNPs<-append(MLH_all_SNPs,MLH_temp)
  MLH_all_SNPs_noNA<-append(MLH_all_SNPs_noNA,MLH_temp1)
}

MLH_all<-cbind(sample_id[ID],as.numeric(MLH_all_SNPs))
MLH_all<-cbind(MLH_all,as.numeric(MLH_all_SNPs_noNA))
MLH_all<-as.data.frame(MLH_all)
MLH_all$V2<-as.numeric(MLH_all$V2)
MLH_all$V3<-as.numeric(MLH_all$V3)
MLH_all$POP<-as.factor(popmap$V2[ID])
colnames(MLH_all)<-c('ID_ind',"MLH_all","MLH_wo_NA","POP")

return(MLH_all)
}

######################  Site frequency spectrum ######################

##################### Folded SFS  by using allele with maf
maf_SFS_1D<-function(genofile,ID_ind){
sfs=NULL
sfs_maf=NULL
geno_temp<-snpgdsGetGeno(genofile,sample.id = ID_ind)
for(j in 1:ncol(geno_temp)){
  sfs_temp<-factor(geno_temp[,j],levels =c(0,1,2,NA))
  sfs_temp1<-table(sfs_temp)
  A<-((sfs_temp1[1]*2)+sfs_temp1[2])
  B<-((sfs_temp1[2])+sfs_temp1[3]*2)
  sfs_temp2<-append(A,B)
  sfs<-cbind(sfs,sfs_temp2)
  sfs_maf<-append(sfs_maf,min(sfs_temp2))}
return(sfs_maf)}

#########################################  Treemix and baypass format ################################

Treemix_Baypass_format<-function(genofile,popmap){
genofile<-genofile
popmap<-popmap
lvls<-c(0,1,2,NA)
format_Baypass=NULL
format_Treemix=NULL
sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
popmap$V2<-as.factor(popmap$V2)
for(pop in levels(popmap$V2)){
  pop_temp<-which(pop==popmap$V2)
  print(pop)
  geno_temp<-snpgdsGetGeno(genofile,sample.id =sample_id[pop_temp])
  freq_tab<-NULL
  for(i in 1:ncol(geno_temp)){
    col_temp<-geno_temp[,i]
    col_temp<-factor(col_temp,levels = lvls)
    tab_temp<-table(col_temp)
    freq_tab<-rbind(freq_tab,tab_temp)}
  pop1<-NULL
  for(j in 1:nrow(freq_tab)){
    t1<-( (as.numeric(freq_tab[j,1])*2)+as.numeric(freq_tab[j,2]))
    t2<-( (as.numeric(freq_tab[j,3])*2)+as.numeric(freq_tab[j,2]))
    temp<-cbind(t1,t2)
    pop1<-rbind(pop1,temp)
    assign(paste0(pop),pop1)}
  pop2<-paste0(pop1[,1],",",pop1[,2])
  format_Baypass<-cbind(format_Baypass,pop1)
  format_Treemix<-cbind(format_Treemix,pop2)
}
format_final<-list(format_Baypass,format_Treemix)
return(format_final)
}










