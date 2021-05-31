file_locat<-getwd()
setwd(file_locat)

library(SNPRelate)
library(gdsfmt)


argv <- commandArgs(TRUE)
x <- argv[1] ### argument two is the file (tab delimited) that contain replicate with their number of couple 
###### with header "Nom_replicates Couple"
###### ind1 1 
###### ind1_r1 1
###### ind2 2
###### ind2_r2 2
vcf_file<- argv[2] ##### VCF 
rep_table<-read.table(x,header = T)
rep_table$Couple<-as.factor(rep_table$Couple)

####### create GDS file from the VCF with snprelate
snpgdsVCF2GDS(vcf_file,paste0(vcf_file,'.gds'), method="biallelic.only")
genofile <- snpgdsOpen(paste0(vcf_file,'.gds'))

####### calculate the missing data 

miss_data<-snpgdsSampMissRate(genofile,with.id = T)
miss_data<-as.data.frame(miss_data)
miss_data$ind<-row.names(miss_data)

id_to_del=NULL
for( i in levels(rep_table$Couple)){
  aa<-which(i==rep_table$Couple)
  #print(length(rep_table$Nom_replicates[aa]))
  long<-length(rep_table$Nom_replicates[aa])
  list_ind<-rep_table$Nom_replicates[aa]
if(long==2){
  list_temp<-NULL
  ind1<-which(miss_data$ind==list_ind[1])
  ind2<-which(miss_data$ind==list_ind[2])
  list_temp<-append(miss_data$miss_data[ind1],miss_data$miss_data[ind2])
  minval<-min(list_temp)
  id_min<-which(list_temp==minval)
  id_to_del_temp<-list_ind[-id_min]
  id_to_del<-append(id_to_del,id_to_del_temp)}
  if(long==3){
    list_temp<-NULL
    ind1<-which(miss_data$ind==list_ind[1])
    ind2<-which(miss_data$ind==list_ind[2])
    ind3<-which(miss_data$ind==list_ind[3])
    list_temp<-append(miss_data$miss_data[ind1],miss_data$miss_data[ind2])
    list_temp<-append(list_temp,miss_data$miss_data[ind3])
    minval<-min(list_temp)
    id_min<-which(list_temp==minval)
    id_to_del_temp<-list_ind[-id_min]
    id_to_del<-append(id_to_del,id_to_del_temp)}}
	
######### output a txt_file with ind name of all replicate with the higher value of missing data. 
write.table(id_to_del,"replicate_with_higher_err_rate.txt",row.names = F, col.names = F, quote = F)

 
 
  
  
  








