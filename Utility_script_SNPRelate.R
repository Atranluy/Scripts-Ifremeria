############" utility script from A.tran LU Y ##############
### Using mainly snprelate packages
###

library(SNPRelate)
library(gdsfmt)
library(inbreedR)


########### calculating multi locus heterozygosity on individual #############
Multilocus_Het_test<-function(genofile,popmap){
genofile<-genofile
popmap<-popmap
sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
ID=NULL
for(i in popmap$V1){
  a<-which(i==sample_id)
  ID<-append(ID,a)
}

number_het<-function(x){
  length(which(x==1))/length(x)}
number_het_noNA<-function(x){
  length(which(x==1))/length(x[!is.na(x)])}

genomatrix<-snpgdsGetGeno(genofile,sample.id =sample_id[ID])
aa<-apply(genomatrix,1,number_het)
ab<-apply(genomatrix,1,number_het_noNA)
MLH_all<-cbind(sample_id[ID],as.numeric(aa))
MLH_all<-cbind(MLH_all,as.numeric(ab))
MLH_all<-as.data.frame(MLH_all)
MLH_all$V2<-as.numeric(MLH_all$V2)
MLH_all$V3<-as.numeric(MLH_all$V3)
MLH_all$POP<-as.factor(popmap$V2[ID])
colnames(MLH_all)<-c('ID_ind',"MLH_all","MLH_wo_NA","POP")
MLH_all<-as.data.frame(MLH_all)
MLH_all<-MLH_all[order(MLH_all$POP),]
return(MLH_all)
}

Multilocus_Het<-function(genofile){
genofile<-genofile
popmap<-popmap
sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
genomatrix<-snpgdsGetGeno(genofile,sample.id =sample_id)
genomatrix[genomatrix==2]<-0
MLH_inbreedR<-MLH(genomatrix)
MLH_inbreedR_final<-cbind(sample.id,MLH_inbreedR)
MLH_inbreedR_final<-as.data.frame(MLH_inbreedR_final)
colnames(MLH_inbreedR_final)<-c("Ind","MLH")
return(MLH_inbreedR_final)
}
######################  Site frequency spectrum ######################

##################### Folded SFS  by using allele by using the maf
maf_SFS_1D<-function(genofile,ID_ind){
sfs=NULL
sfs_maf=NULL
geno_temp<-snpgdsGetGeno(genofile,sample.id = ID_ind)
for(j in 1:ncol(geno_temp)){
  sfs_temp<-factor(geno_temp[,j],levels =c(0,1,2,NA))
  sfs_temp1<-table(sfs_temp)
  B<-((sfs_temp1[1]*2)+sfs_temp1[2])
  A<-((sfs_temp1[2])+sfs_temp1[3]*2)
  sfs_temp2<-append(A,B)
  sfs<-cbind(sfs,sfs_temp2)
  sfs_maf<-append(sfs_maf,min(sfs_temp2))}
return(sfs_maf)}


#########################################  Treemix and baypass format ################################

Treemix_Baypass_format<-function(genofile,popmap){
genofile<-genofile
popmap<-as.data.frame(popmap)
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
    t2<-( (as.numeric(freq_tab[j,1])*2)+as.numeric(freq_tab[j,2]))
    t1<-( (as.numeric(freq_tab[j,3])*2)+as.numeric(freq_tab[j,2]))
    temp<-cbind(t1,t2)
    pop1<-rbind(pop1,temp)
    assign(paste0(pop),pop1)}
  pop2<-paste0(pop1[,1],",",pop1[,2])
  format_Baypass<-cbind(format_Baypass,pop1)
  format_Treemix<-cbind(format_Treemix,pop2)
}
colnames(format_Treemix)<-levels(popmap$V2)
format_final<-list(Baypass=format_Baypass,Treemix=format_Treemix)
return(format_final)
}


##############  heterozygosity per SNP ##################
#!!!! pas fini !!!!! #


Het_per_SNPs_<-function(genofile,popmap){
if (is.null(popmap)){
genomatrix<-snpgdsGetGeno(genofile)
number_het<-function(x){
  length(which(x==1))/length(x)}
number_het_noNA<-function(x){
  length(which(x==1))/length(x[!is.na(x)])}
aa<-apply(genomatrix,2,number_het)
ab<-apply(genomatrix,2,number_het_noNA)
write.table(aa,'het_per_snp_with_na.txt',col.names = F,row.names = F,quote = F)
write.table(ab,'het_per_snp_with_NOna.txt',col.names = F,row.names = F,quote = F)}
  else{
    genofile<-genofile
    popmap<-popmap
    sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
    ID=NULL
    for(i in popmap$V1){
      a<-which(i==sample_id)
      ID<-append(ID,a)
    }
    genomatrix<-snpgdsGetGeno(genofile,sample.id =sample_id[ID])
    number_het<-function(x){
      length(which(x==1))/length(x)}
    number_het_noNA<-function(x){
      length(which(x==1))/length(x[!is.na(x)])}
    aa<-apply(genomatrix,2,number_het)
    ab<-apply(genomatrix,2,number_het_noNA)
    write.table(aa,'het_per_snp_with_na.txt',col.names = F,row.names = F,quote = F)
    write.table(ab,'het_per_snp_with_NOna.txt',col.names = F,row.names = F,quote = F)
  }
  }








