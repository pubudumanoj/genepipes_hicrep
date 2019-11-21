# Evaluates reproducibility of Hi-C data using hicrep
# Inputs: two hic-interaction matrices to compare (defined by the design file)
# Written by Pubudu Nawarathna - Nov 2019
# Usage : Rscript hicrep.R -d path_design -c path_interaction_matrices -o output_dir

library(data.table)
library(hicrep)
#library(dplyr)

# Usage

usage = function(errM) {
  cat("\nUsage : Rscript hicrep.R [option] <Value>\n")
  cat("       -d        : design file\n")
  cat("       -c        : raw count file\n")
  cat("       -o        : output directory\n")
  cat("       -h        : this help\n\n")
  
  stop(errM)
}

set.seed(123456789)

#The inputs for perform_hicrep function are simply two Hi-C matrices to be compared.
#The Hi-C matrices should have the dimension N×(3+N).
#The three additional initial columns are the chromosome name
#and mid-point coordinates of the two contact bins

load_file <- function(fname=NULL){
  
  return(fread(fname))
}

perform_hicrep <- function(mat1=NULL, mat2=NULL, bin=1000000, smooth=1, boundary=5000000){
  
  #mat1 = mat1[-1,]
  #mat2 = mat2[-1,]
  mat1 <- load_file(mat1)
  mat2 <- load_file(mat2)
  
  #Rename the first column
  colnames(mat1)[1] <- "var1"
  colnames(mat2)[1] <- "var1"
  
  #create a new column with chromosome name
  mat1$chr <- sub("\\-.*", "", mat1$var1)
  mat2$chr <- sub("\\-.*", "", mat2$var1)
  
  #create a new column with bin start value
  mat1$bin_start <- as.numeric(sub(".*\\-", "", mat1$var1))
  mat2$bin_start <- as.numeric(sub(".*\\-", "", mat2$var1))
  
  #create a new column with bin end value
  mat1$bin_end <- as.numeric(mat1$bin_start + bin)
  mat2$bin_end <- as.numeric(mat2$bin_start + bin)
  
  mat1 <- as.data.frame(mat1)
  mat2 <- as.data.frame(mat2)
  
  #Restructure the data frame acoording to the hicrep input format
  #Required structure of the data set
  ##       V1      V2       V3 V4 V5 V6 V7 V8 V9 V10
  ## 1  chr22       0  1000000  0  0  0  0  0  0   0
  ## 2  chr22 1000000  2000000  0  0  0  0  0  0   0
  ## 3  chr22 2000000  3000000  0  0  0  0  0  0   0
  ## 4  chr22 3000000  4000000  0  0  0  0  0  0   0
  ## 5  chr22 4000000  5000000  0  0  0  0  0  0   0
  ## 6  chr22 5000000  6000000  0  0  0  0  0  0   0
  ## 7  chr22 6000000  7000000  0  0  0  0  0  0   0
  ## 8  chr22 7000000  8000000  0  0  0  0  0  0   0
  ## 9  chr22 8000000  9000000  0  0  0  0  0  0   0
  ## 10 chr22 9000000 10000000  0  0  0  0  0  0   0
  
  
  mat1 <- mat1[,c(dim(mat1)[2]-2,dim(mat1)[2]-1,dim(mat1)[2], 3:(dim(mat1)[2]-3))]
  mat2 <- mat2[,c(dim(mat2)[2]-2,dim(mat2)[2]-1,dim(mat2)[2], 3:(dim(mat2)[2]-3))]
  
  #mat1 <- as.data.frame(mat1)
  #mat2 <- as.data.frame(mat2)


  #rownames(mat1) = make.names(1:dim(mat1)[1], unique=TRUE)
  #rownames(mat2) = make.names(1:dim(mat2)[1], unique=TRUE)
  Pre_HiC <- prep(mat1, mat2, bin, smooth, boundary)
  
  SCC.out = get.scc(Pre_HiC, bin, boundary)
  return(SCC.out[[3]])
}


hicrep_all_chr <- function(){
  
  #read all files from two different samples and send two files at a time to the function perform_hicrep
}
#mat1 <- fread("interaction_matrices/SRR1658581/genomewideMatrices/HTD_SRR1658581_MboI_genomewide_Res_1000000_raw.txt", header = F)
#mat1 <- fread("/project/6002326/pubudu/job_outputs/Fiona_hicseq/chromosomeMatrices_MSC16_AD_R1_HiC/HTD_MSC16_AD_R1_HiC_Arima_chr10_50000_raw.txt")
#mat2 <- fread("/project/6002326/pubudu/job_outputs/Fiona_hicseq/chromosomeMatrices_MSC72_AD_R1_HiC/HTD_MSC72_AD_R1_HiC_Arima_chr10_50000_raw.txt")
#perform_hicrep(mat1, mat2, bin=50000, boundary = 500000)
sample1="MSC16_AD_R1_HiC"
sample2="MSC72_AD_R1_HiC"
#select files inside the "interaction_matrices" which are located in "chromosomeMatrices" for given sample1 and 2 (for the comparison)
file_list_sample1 <- list.files(path=paste0("interaction_matrices/",sample1), pattern = "_raw.txt", recursive = TRUE)
file_list_sample2 <- list.files(path=paste0("interaction_matrices/",sample2), pattern = "_raw.txt", recursive = TRUE)
#only select files which have phrase "chromosomeMatrices"

file_list_sample1 <- file_list_sample1[ grepl("chromosomeMatrices", file_list_sample1) ]
file_list_sample2 <- file_list_sample2[ grepl("chromosomeMatrices", file_list_sample2) ]

#ignore all files from "genomewideMatrices" folders
file_list_sample1 <- file_list_sample1[ !grepl("genomewideMatrices", file_list_sample1) ]
file_list_sample2 <- file_list_sample2[ !grepl("genomewideMatrices", file_list_sample2) ]

length(file_list_sample1)




#get after
#body=`echo "$entry" | awk -F_ '{print $NF}'`
#rm_end=echo "$body" | cut -d_r1.fastq -f1`
