# **Readme for RScripts**


## **Calculating genotype error rate**
the Rscript **Error_rate_replicate_VCF.R** will calculate the genotyped error rate between replicate from a VCF file.

It will use the package R SNPrelate.

This will simply calculate the number of non-concordant genotype between replicate in percent with 4 differents way.

### 4 Differents table of genotyped error rate 
1. **Table_geno_error.csv** is a table of the percent of genotyped error over all SNPs (including missing data *NA* but not considered as an error)
2. **Table_geno_error_with_NA.csv** is a table of the percent of genotyped error over all SNPs (including missing data *NA* and considered as an error)
3. **Table_geno_error_with_noNA.csv** is a table of the percent of genotyped error over all tagged SNPs for both replicate (missing data is not taken in account)
4. **Table_geno_error_only_NA.csv** is a table of the percent of genotyped error while considering missing data as the only error. 
_____

To run this script on a HPC, you just have to write "Rscript Error_rate_replicate_VCF.R {
