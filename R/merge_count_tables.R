## Obtaining common count table for TElocal output

library(tidyverse)
library(data.table)
library(stringr)
library(purrr)
library(zeallot)


## The following functions work appropriately provided that TE annotation in count tables has a form of "TE_ID:TE_family:TEclass" (standard for TElocal output)
## And the annotations for genes DO NOT contain ":" (e.g. "ENSG01101100","LOL1","ENSMUS01101100").


## Extract counts for genes (entries with annotation that DOES NOT contain ":" [common gene names, ENSEMBL IDs etc])
genesep<-function(df){
  ge_cnt<-df[df$gene.TE %in% grep(":",df$gene.TE,value=TRUE, invert=TRUE),]
  return(ge_cnt)
}

# Extract counts for TEs only (entries with annotation that DOES contain ":", e.g. "hAT-N1_Mam_dup9:hAT-Tip100:DNA")
tesep<-function(df){
  te_cnt<-df[df$gene.TE %in% grep(":",df$gene.TE,value=TRUE),]
  return(te_cnt)
}


# telocal_load is function for reading in and creating a common table for all telocal count tables under one, provided as input argument, directory path
# Example of dir-input: "/my/dir/with/counts/"

# loads and merges all the cntTable tables, generating three separate data.tables, containing loci, subfamily and gene level counts
# To save all three results, use zeallot and following pattern of command:
# c(lociCnt,subfCnt, geneCnt)%<-%te_load(mypath)

telocal_load<-function(mypath){
  temp<-list.files(path=mypath, pattern="*.cntTable")
  temp<-paste0(mypath,temp)
  myfiles<-lapply(temp,read.delim)
  mytes<-lapply(myfiles,tesep)
  mygenes<-lapply(myfiles,genesep)
  tel_tab<-mytes %>% purrr::reduce(full_join, by="gene.TE")
  tel_tab<-as.data.table(tel_tab,key="gene.TE")
  colnames(tel_tab)[1]<-"ID"
  tel_tab[,ID:=str_split_fixed(ID,":" ,2)[,1]]
  tel_tab2<-separate(tel_tab, ID, sep="_dup", into=c("subid","ID"))
  tel_tab2[,ID:=NULL]
  te_sub<-colnames(tel_tab2)[!colnames(tel_tab2) %in% c("ID","subid")]
  te_sub2<-tel_tab2[,lapply(.SD,sum), by = .(subid),.SDcols=te_sub]
  
  gene_tab<-mygenes %>% purrr::reduce(full_join, by="gene.TE")
  gene_tab<-as.data.table(gene_tab,key="gene.TE")
  return(list(tel_tab,te_sub2,gene_tab))
}

# Example of run command
# c(TEloci_tab,TEsub_tab,Genes_tab)%<-%telocal_load("/absolute/path/to/your/TElocal_count_files/")

# Save tables using fwrite, write.csv or anything else
