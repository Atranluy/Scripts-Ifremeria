file_locat<-getwd()
setwd(file_locat)
argv <- commandArgs(TRUE)
#script Adrien

VCF<-argv[1]
rep_file <- argv[2]

library(SNPRelate)
library(gdsfmt)





File_repl<-read.table("Replicas_Txt_all.txt", sep="\t", header=T)
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
  
  
  
  
  
  
  #list_ind<-c("In188-2a","In188-2b","In218","In218_Rinter_1","In218_Rinter_2","In227","In227_Rinter_1","In232","In232_Rinter_1","In286-2","In286-2_Rinter_1","In366-2","In366-2_Rinter_1","In381-2","In381-2_Rinter_1","In390-2","In390-2_Rinter","In391-2","In391-2_Rinter","In392-2","In392-2_Rinter","In399-2","In399-2_Rinter_1","In002-3","In002-3_Rintra_1","In093-2","In093-2_Rinter","In154-2","In154-2_Rinter_1",'In155-2',"In155_2_Rinter_1","In422-2","In422-2_Rinter","In429-2","In429-2_Rinter")
  

  
  

  list_ind<-c("In093-2", "In093-2_Rinter", "In154-2", "In154-2_Rinter_1", "In155_2_Rinter_1", "In155-2", "In188-2a", "In188-2b", "In218", "In218_Rinter_1", "In218_Rinter_2", "In227", "In227_Rinter_1", "In232", "In232_Rinter_1", "In366-2", "In366-2_Rinter_1", "In390-2", "In390-2_Rinter", "In391-2", "In391-2_Rinter", "In392-2", "In392-2_Rinter", "In399-2", "In399-2_Rinter_1", "In422-2", "In422-2_Rinter", "In429-2", "In429-2_Rinter", "In561", "In561_Rinter_1", "In561_Rintra_1", "In661", "In661_Rinter_1", "In661_Rintra_1","In002-3","In002-3_Rintra_1","In286-2","In286-2_Rinter_1","In381-2","In381-2_Rinter_1","In435", "In435_Rinter_1","In487-2_Rinter_1","In487","In647","In647_Rinter_1", "In400-2","In400-3","In390-3","In391-3","In392-3")
all_taux<-NULL
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
    print(paste0("couple de replicas ",ind1," et ",ind2))
    
    
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
    print(paste0("couple de replicas ",ind1," et ",ind2))
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
    print(paste0("couple de replicas ",ind1," et ",ind2))
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
    print(paste0("couple de replicas ",ind1," et ",ind2))
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
  
  write.table(t2,"all_taux_erreur_geno.csv",sep=";", row.names = F, col.names = T,quote=F)
  write.table(t1,"all_taux_erreur.csv",sep=";", row.names = F, col.names = F)
  write.table(t3,"all_taux_no_NA_erreur.csv",sep=";", row.names = F, col.names = F)
  write.table(t4,"all_taux_onlyNAerreur.csv",sep=";", row.names = F, col.names = F)


mddd<-(kl[1]+nb_diff_withno_NA)
mddd/nombre_de_row
  
  