file_locat<-getwd()
setwd(file_locat)
argv <- commandArgs(TRUE)
#script Adrien Tran Lu Y 

VCF<-argv[1]  # argument one is the vcf
rep_file <- argv[2] # argument two is the file (tab delimited) that contain replicate with their number of couple 

library(SNPRelate)
library(gdsfmt)

########### Rscript which calculate genotyping error from VCF with 4 differents consideration of missing data.

File_repl<-read.table(rep_file, sep="\t", header=T)
File_repl$Couple<-as.factor(File_repl$Couple)

filenames <- Sys.glob("*snps.vcf")
nom_file<-NULL
t1<-NULL
t2<-NULL
t3<-NULL
t4<-NULL
  vcf.fn<-filenames
  snpgdsVCF2GDS(vcf.fn,paste0(filenames,'.gds'), method="biallelic.only")
  genofile <- snpgdsOpen(paste0(filenames,'.gds'))
  geno_ma<-snpgdsGetGeno(genofile,with.id = T)
  geno_ma1<-t(geno_ma$genotype)
  nombre_de_row<-nrow(geno_ma1)
  colnames(geno_ma1)<-geno_ma$sample.id
  colnames(geno_ma1)
  
  

all_taux_geno<-NULL
all_taux_no_NA<-NULL
all_taux_onlyNA<-NULL
list_nom_rep<-c() 
  
for( i in levels(File_repl$Couple)){
  aa<-which(i==File_repl$Couple)
  print(length(File_repl$Nom_replicas[aa]))
  long<-length(File_repl$Nom_replicas[aa])
  list_ind<-File_repl$Nom_replicas[aa]
  if(long==2){
    #########
    ind1<-which(colnames(geno_ma1)==list_ind[1])
    ind2<-which(colnames(geno_ma1)==list_ind[2])
    print(paste0("pair of replicate  ",ind1," and ",ind2))
    
    
    k<-geno_ma1[,ind1]==geno_ma1[,ind2]
    kl<-table(k,exclude = F )
    tot_nodiff<-sum(kl)
    nb_diff_withno_NA<-(nombre_de_row-tot_nodiff)
    taux_err_genotyped<-((nb_diff_withno_NA/(kl[1]+nb_diff_withno_NA))*100)
    taux_err<-((nb_diff_withno_NA/nombre_de_row)*100)
    taux_err_with_NA<-(((kl[2]+nb_diff_withno_NA)/nombre_de_row)*100)
    taux_err_only_NA<-((kl[2]/nombre_de_row)*100)
    all_taux<-rbind(all_taux,taux_err_with_NA)
    all_taux_no_NA<-rbind(all_taux_no_NA,taux_err)
    all_taux_geno<-rbind(all_taux_geno,taux_err_genotyped)
    all_taux_onlyNA<-rbind(all_taux_onlyNA,taux_err_only_NA)
    
    list_nom_rep<-c(list_nom_rep,paste0(list_ind[1],'_',list_ind[2]))}
  if(long==3){
    #########
    ind1<-which(colnames(geno_ma1)==list_ind[1])
    ind2<-which(colnames(geno_ma1)==list_ind[2])
    print(paste0("pair of replicate ",ind1," and ",ind2))
    k<-geno_ma1[,ind1]==geno_ma1[,ind2]
    kl<-table(k,exclude = F )
    tot_nodiff<-sum(kl)
    nb_diff_withno_NA<-(nombre_de_row-tot_nodiff)
    taux_err_genotyped<-((nb_diff_withno_NA/(kl[1]+nb_diff_withno_NA))*100)
    taux_err<-((nb_diff_withno_NA/nombre_de_row)*100)
    taux_err_with_NA<-(((kl[2]+nb_diff_withno_NA)/nombre_de_row)*100)
    taux_err_only_NA<-((kl[2]/nombre_de_row)*100)
    all_taux<-rbind(all_taux,taux_err_with_NA)
    all_taux_no_NA<-rbind(all_taux_no_NA,taux_err)
    all_taux_geno<-rbind(all_taux_geno,taux_err_genotyped)
    all_taux_onlyNA<-rbind(all_taux_onlyNA,taux_err_only_NA)
    list_nom_rep<-c(list_nom_rep,paste0(list_ind[1],'_',list_ind[2]))
    ind1<-which(colnames(geno_ma1)==list_ind[1])
    ind2<-which(colnames(geno_ma1)==list_ind[3])
    print(paste0("pair of replicate ",ind1," and ",ind2))
    k<-geno_ma1[,ind1]==geno_ma1[,ind2]
    kl<-table(k,exclude = F )
    tot_nodiff<-sum(kl)
    nb_diff_withno_NA<-(nombre_de_row-tot_nodiff)
    taux_err_genotyped<-((nb_diff_withno_NA/(kl[1]+nb_diff_withno_NA))*100)
    taux_err<-((nb_diff_withno_NA/nombre_de_row)*100)
    taux_err_with_NA<-(((kl[2]+nb_diff_withno_NA)/nombre_de_row)*100)
    taux_err_only_NA<-((kl[2]/nombre_de_row)*100)
    all_taux<-rbind(all_taux,taux_err_with_NA)
    all_taux_no_NA<-rbind(all_taux_no_NA,taux_err)
    all_taux_geno<-rbind(all_taux_geno,taux_err_genotyped)
    all_taux_onlyNA<-rbind(all_taux_onlyNA,taux_err_only_NA)
    list_nom_rep<-c(list_nom_rep,paste0(list_ind[1],'_',list_ind[3]))                 
    ind1<-which(colnames(geno_ma1)==list_ind[2])
    ind2<-which(colnames(geno_ma1)==list_ind[3])
    print(paste0("pair of replicate ",ind1," and ",ind2))
    k<-geno_ma1[,ind1]==geno_ma1[,ind2]
    kl<-table(k,exclude = F )
    tot_nodiff<-sum(kl)
    nb_diff_withno_NA<-(nombre_de_row-tot_nodiff)
    taux_err_genotyped<-((nb_diff_withno_NA/(kl[1]+nb_diff_withno_NA))*100)
    taux_err<-((nb_diff_withno_NA/nombre_de_row)*100)
    taux_err_with_NA<-(((kl[2]+nb_diff_withno_NA)/nombre_de_row)*100)
    taux_err_only_NA<-((kl[2]/nombre_de_row)*100)
    all_taux<-rbind(all_taux,taux_err_with_NA)
    all_taux_no_NA<-rbind(all_taux_no_NA,taux_err)
    all_taux_geno<-rbind(all_taux_geno,taux_err_genotyped)
    all_taux_onlyNA<-rbind(all_taux_onlyNA,taux_err_only_NA)
    
    list_nom_rep<-c(list_nom_rep,paste0(list_ind[2],'_',list_ind[3]))
    }
  
}
  row.names(all_taux)<-list_nom_rep
  row.names(all_taux_geno)<-list_nom_rep
  row.names(all_taux_no_NA)<-list_nom_rep
  row.names(all_taux_onlyNA)<-list_nom_rep
  ##########
  
  tt1<-t(all_taux)
  tt2<-t(all_taux_geno)
  tt3<-t(all_taux_no_NA)
  tt4<-t(all_taux_onlyNA)
  
  t1<-rbind(t1,tt1)
  t2<-rbind(t2,tt2)
  t3<-rbind(t3,tt3)
  t4<-rbind(t4,tt4)
  #write.table(nom_file,"nom_file.csv",sep=";", row.names = F, col.names = T)
  
  
  
  ####### write genotyping error rate only for genotyped SNPs for both replicate over all SNPs (no  missing data)
  write.table(t2,"all_taux_erreur_geno.csv",sep=";", row.names = F, col.names = T,quote=F)
  
  ####### write genotyping error rate only for genotyped SNPs for both replicate over all SNPs (consideringg)
  write.table(t1,"all_taux_erreur.csv",sep=";", row.names = F, col.names = F)
  
  write.table(t3,"all_taux_no_NA_erreur.csv",sep=";", row.names = F, col.names = F)
  
  write.table(t4,"all_taux_onlyNAerreur.csv",sep=";", row.names = F, col.names = F)



  