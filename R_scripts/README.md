# ***Readme for RScripts***


##  **Calculating genotype error rate**
the Rscript **Error_rate_replicate_VCF.R** will calculate the genotyped error rate between replicate from a VCF file.  
It will use the package R SNPrelate.  
This will simply calculate the number of non-concordant genotype between replicate in percent with 4 differents way.  

### 4 Differents table of genotyped error rate 
1. **Table_geno_error.csv** is a table of the percent of genotyped error over all SNPs (including missing data *NA* but not considered as an error)
2. **Table_geno_error_with_NA.csv** is a table of the percent of genotyped error over all SNPs (including missing data *NA* and considered as an error)
3. **Table_geno_error_with_noNA.csv** is a table of the percent of genotyped error over all tagged SNPs for both replicate (missing data is not taken in account)
4. **Table_geno_error_only_NA.csv** is a table of the percent of genotyped error while considering missing data as the only error.
_____

To run this script on a HPC, you just have to write :
```
Rscript Error_rate_replicate_VCF.R my.vcf file_replicate.txt
```
* my.vcf just correspond to the name of the vcf (e.g my vcf is called  *populations.snps.vcf*, you will just write *populations.snps.vcf*)
* file_replicate.txt is just a text file with 2 tab-delimited column.
*  The first one correspond of the same sample id as used in the VCF file and the second correspond to a integer of pair of replicate. an exemple below
### Exemple
Replicates|Pairs
----------|----------
name_ind1|1
name_ind1_rep|1
name_ind2|2
name_ind2_rep|2
name_ind3|3
name_ind3_rep1|3
name_ind3_rep2|3

## **Rscript to remove replicate with higher rate of missing data**
the Rscrip **Remove_replicate_missdata.R** will calculate the percent of missing data for all replicate from a VCF file and return the list of replicate with the highest value.
This script run simply on a HPC by using :
```
Rscript Remove_replicate_missdata.R my.vcf file_replicate.txt
 ```
This script will need the same input as described in the previous script.



