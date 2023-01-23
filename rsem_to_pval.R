# 
# run deseq2 for pairwise statistical evaluation 
# of nf-core rna-seq pipeline output
# specifically for use on the ***merged gene counts tsv***
# 


#step 1
# make metadata table for deseq that looks like this 
# 


# #             id sampletype
# 1  cont_rnai_1    control   
# 2  cont_rnai_2    control
# 3  cont_rnai_3    control
# 4 sfmbt_rnai_1     mutant
# 5 sfmbt_rnai_2     mutant
# 6 sfmbt_rnai_3     mutant

# id needs to match rsem.merged.gene_counts.tsv column headers
# check your input like
# countData <- read.csv('http://germarium.ddns.net/nextflow_rui_sfmbt_fix/star_rsem/rsem.merged.gene_counts.tsv', header = TRUE, sep = "\t")
# names(countData)

# deseq2 will compare the control to mutant, based on ~sampletype call
# 
# the top group is the denominator of foldchange calc,
# so be sure to put controls at top of metaadata table

# for example, to make this table of your samples, 
# do something like this
# assuming control samples are the leftmost columns
# columns are alphanumerically arranged, caveat emptor
# note, control needs to called 'control' and mutant = 'mutant'
# for the deseq2 script to work



#countData <- read.csv('http://germarium.ddns.net/nextflow/rna/nextflow_rui_sfmbt_fix/star_rsem/rsem.merged.gene_counts.tsv', header = TRUE, sep = "\t")
# 
# an example where control comes before mut
metaData<-data.frame(id=colnames(countData[3:8]),
                     sampletype=c('control','control','control',
                                  'mutant','mutant','mutant'))

# an example where control comes after mut
countData <- read.csv('http://germarium.ddns.net/nextflow/rna/nextflow_rsem_suvar210_ninova_etAl/star_rsem/rsem.merged.gene_counts.tsv', header = TRUE, sep = "\t")

metaData<-data.frame(id=colnames(countData[,c(7:9,3:6)]), 
                     sampletype=c('control','control','control', 'mutant','mutant','mutant','mutant'
                                  
                                  ))





###########

###########

# step 2: run deseq2 
 # how to run:
#  outputFile<-rsem_genes_to_pval(countData)




rsem_genes_to_pval<-function (rsem_file_from_pipeline, #put path in quotes
                              metaData, 
                              mutantPrefix="mutant", #needs to be specific enough, a la grep
                              controlPrefix="control",   #needs to be specific enough
                              alpha = 0.05){
  #alpha is pval cutoff for sig test
  countData <- read.csv(rsem_file_from_pipeline, header = TRUE, sep = "\t")

  
  countData<-countData[,-2]
  
  # deseq2 wants integers not floats, and so we round the estimated count data with this 
  # expanding the number of columns included as needed
  rdcountData<-data.frame(gene_id=countData$gene_id, 
                          round(countData[,c(2:length(colnames(countData)))],0))
  
  # modify metadata table if more than 1 experimental (mutant)
  # or if want compare mutant1 v mutant2 (by putting mutant in place of controlPrefix)
  
 # metaData<-metaData[ grepl(mutantPrefix, metaData$id) | 
  #                      grepl(controlPrefix, metaData$id) , ]
  
  
  rdcountData<-as_tibble(rdcountData)
  
  # grab the columns you need, which will match metaData
  # control and mutant variable come down as quoted here, so no quotes needed  
  subsetCountData<-rdcountData %>% dplyr::select( c(any_of("gene_id")| 
                                                      starts_with( controlPrefix ) | 
                                                      starts_with( mutantPrefix )))   
  subsetCountData<-as.data.frame(subsetCountData)
    
  dds <- DESeqDataSetFromMatrix(countData=subsetCountData, 
                                colData=metaData, 
                                design=~sampletype, tidy = TRUE)
  
  
  dds = DESeq(dds)
  
  res <- results(dds, alpha = 0.05)
  
  
  summary(res)
  
  resLFC = lfcShrink(dds, coef='sampletype_mutant_vs_control')
  output<-as.data.frame(resLFC)
  head(resLFC)
  
  output<-data.frame(output, tx_id=rownames(res))
  
  return(output)
  
}
