---
title: "032720_RNAseq_analysis_PART3"
output: 
  html_document:
    keep_md: true
editor_options: 
  chunk_output_type: console
---



## Purpose

The goal of this section is to predict which transcription factors may regulate Card9-dependent and Card9-independent Dectin-1 regulated genes. To do this, we are using ImmGen ATAC-seq data to identify open chromatin regions in our cell type of interest (BM neutrophils) and the predicted transcription factor (TF) binding sites in these open regions. This information is then used with a pathway enrichment analysis approach to identify those predicted TF sites that are overrepresented near genes of one group or the other. 


```r
library(devtools)
library(Biobase)
library(goseq)
library(DESeq2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(geneLenDataBase)
library(tidyr)
library(dplyr)
library(readr)
library(VennDiagram)
library(ggplot2)
library(ggrepel)
library(plotly)
```

## 1. Compile ATACseq OCR data from ImmGen files

Processed data files were downloaded from ImmGen for BM neutrophils for each chromosome. Here we are compiling these separate chr files into a single dataframe.


```r
#loaded all csv files into a single file
files <- list.files(path = "Immgen_ATACseq_chromosomes", pattern = "*.csv", full.names = T)
tbl <- sapply(files, read_csv, simplify=FALSE) %>% 
bind_rows(.id = "id")

#493013 rows, 17 columns
dim(tbl)
```

```
## [1] 493013     17
```

```r
#modifying name, class, factors
tbl$Chr <- as.factor(tbl$Chr)
nlevels(tbl$Chr)
```

```
## [1] 20
```

```r
colnames(tbl)
```

```
##  [1] "id"                              "OCR ID"                         
##  [3] "Chr"                             "OCR Start"                      
##  [5] "OCR End"                         "OCR Summit"                     
##  [7] "mm10_60way_phastCons_scores"     "minus_log10_bestPvalue"         
##  [9] "Included in systematic analysis" "TSS"                            
## [11] "Genes within 100Kb"              "Cis associated genes"           
## [13] "signedLogPvalues"                "TF ID"                          
## [15] "TF Name"                         "Score"                          
## [17] "GN_BM"
```

```r
tbl$Genes.within.100kb <- tbl$`Genes within 100Kb`
tbl <- as.data.frame(tbl)
class(tbl)
```

```
## [1] "data.frame"
```

```r
colnames(tbl)
```

```
##  [1] "id"                              "OCR ID"                         
##  [3] "Chr"                             "OCR Start"                      
##  [5] "OCR End"                         "OCR Summit"                     
##  [7] "mm10_60way_phastCons_scores"     "minus_log10_bestPvalue"         
##  [9] "Included in systematic analysis" "TSS"                            
## [11] "Genes within 100Kb"              "Cis associated genes"           
## [13] "signedLogPvalues"                "TF ID"                          
## [15] "TF Name"                         "Score"                          
## [17] "GN_BM"                           "Genes.within.100kb"
```

## 2. Generate a tidy dataset, first subsetting on OCRs with values >10


```r
#CUTOFF - Remove rows where GN_BM is <10
OCR_threshhold <- 10
tbl_cutoff <- subset(tbl, tbl[,'GN_BM'] > OCR_threshhold)  

#38997 rows, 18 columns
dim(tbl_cutoff)
```

```
## [1] 38997    18
```

Gene names into rows


```r
#making tidy - splitting comma separated gene names into rows - took ~2min to run
start_time1 <- Sys.time()
tbl_1 <- tbl_cutoff %>% 
    mutate(Genes.within.100kb = strsplit(as.character(Genes.within.100kb), ",")) %>% 
    unnest(Genes.within.100kb)
end_time1 <- Sys.time()
end_time1 - start_time1
```

```
## Time difference of 1.18684 mins
```

```r
dim(tbl_1)
```

```
## [1] 184577     18
```

TF names into rows


```r
#making tidy - splitting transcription factor names into rows
#try replacing | with space
tbl_2 <- gsub("\\|", " ", tbl_1$`TF Name`)
class(tbl_2)
```

```
## [1] "character"
```

```r
#created character vector now add as new column to original dataframe
tbl_1$TF_name_space <- tbl_2


#separating each letter instead of each TF name.... took a very long time to run
start_time2 <- Sys.time()
tbl_3 <- tbl_1 %>% 
    mutate(TF_name_space = strsplit(as.character(TF_name_space), " ")) %>% 
    unnest(TF_name_space)
end_time2 <- Sys.time()
end_time2 - start_time2
```

```
## Time difference of 1.809073 hours
```

Looking at features of and saving the compiled file


```r
#now ~5 million rows long (with tidy gene + tf combinations)
dim(tbl_3)
```

```
## [1] 4795672      19
```

```r
head(tbl_3)
```

```
## # A tibble: 6 x 19
##   id    `OCR ID` Chr   `OCR Start` `OCR End` `OCR Summit` mm10_60way_phas…
##   <chr> <chr>    <fct>       <dbl>     <dbl>        <dbl>            <dbl>
## 1 Immg… OCR_379… chr10     3271113   3271363      3271238                0
## 2 Immg… OCR_379… chr10     3271113   3271363      3271238                0
## 3 Immg… OCR_379… chr10     3271113   3271363      3271238                0
## 4 Immg… OCR_379… chr10     3271113   3271363      3271238                0
## 5 Immg… OCR_379… chr10     3271113   3271363      3271238                0
## 6 Immg… OCR_379… chr10     3271113   3271363      3271238                0
## # … with 12 more variables: minus_log10_bestPvalue <dbl>, `Included in
## #   systematic analysis` <chr>, TSS <chr>, `Genes within 100Kb` <chr>,
## #   `Cis associated genes` <chr>, signedLogPvalues <dbl>, `TF ID` <chr>,
## #   `TF Name` <chr>, Score <chr>, GN_BM <dbl>, Genes.within.100kb <chr>,
## #   TF_name_space <chr>
```

```r
colnames(tbl_3)
```

```
##  [1] "id"                              "OCR ID"                         
##  [3] "Chr"                             "OCR Start"                      
##  [5] "OCR End"                         "OCR Summit"                     
##  [7] "mm10_60way_phastCons_scores"     "minus_log10_bestPvalue"         
##  [9] "Included in systematic analysis" "TSS"                            
## [11] "Genes within 100Kb"              "Cis associated genes"           
## [13] "signedLogPvalues"                "TF ID"                          
## [15] "TF Name"                         "Score"                          
## [17] "GN_BM"                           "Genes.within.100kb"             
## [19] "TF_name_space"
```

```r
#MAKE my columns of interest factors
tbl_3$TF_name_space <- as.factor(tbl_3$TF_name_space)
tbl_3$Genes.within.100kb <- as.factor(tbl_3$Genes.within.100kb)
#number of levels: 1205 TFs, 18071 genes
nlevels(tbl_3$TF_name_space)
```

```
## [1] 1205
```

```r
nlevels(tbl_3$Genes.within.100kb)
```

```
## [1] 18071
```

```r
#Save file
start_time3 <- Sys.time()
save(tbl_3, file = "Immgen_ATACSeq_GNBM_OCR_Compiled_GN10cutoff.Rdata")
end_time3 <- Sys.time()
end_time3 - start_time3
```

```
## Time difference of 29.90925 secs
```

## SessionInfo


```r
sessionInfo()
```

```
## R version 3.6.0 (2019-04-26)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Mojave 10.14.3
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
##  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] plotly_4.9.0                ggrepel_0.8.1              
##  [3] ggplot2_3.2.1               VennDiagram_1.6.20         
##  [5] futile.logger_1.4.3         readr_1.3.1                
##  [7] dplyr_0.8.3                 tidyr_1.0.0                
##  [9] org.Mm.eg.db_3.8.2          org.Hs.eg.db_3.8.2         
## [11] AnnotationDbi_1.47.1        DESeq2_1.25.14             
## [13] SummarizedExperiment_1.15.9 DelayedArray_0.11.8        
## [15] BiocParallel_1.19.3         matrixStats_0.55.0         
## [17] GenomicRanges_1.37.16       GenomeInfoDb_1.21.2        
## [19] IRanges_2.19.16             S4Vectors_0.23.25          
## [21] goseq_1.37.0                geneLenDataBase_1.21.0     
## [23] BiasedUrn_1.07              Biobase_2.45.1             
## [25] BiocGenerics_0.31.6         devtools_2.2.1             
## [27] usethis_1.5.1              
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.4-1         ellipsis_0.3.0          
##   [3] rprojroot_1.3-2          htmlTable_1.13.2        
##   [5] XVector_0.25.0           base64enc_0.1-3         
##   [7] fs_1.3.1                 rstudioapi_0.10         
##   [9] remotes_2.1.0            bit64_0.9-7             
##  [11] fansi_0.4.0              splines_3.6.0           
##  [13] geneplotter_1.63.0       knitr_1.25              
##  [15] pkgload_1.0.2            zeallot_0.1.0           
##  [17] jsonlite_1.6             Formula_1.2-3           
##  [19] Rsamtools_2.1.6          annotate_1.63.0         
##  [21] cluster_2.1.0            GO.db_3.8.2             
##  [23] dbplyr_1.4.2             compiler_3.6.0          
##  [25] httr_1.4.1               backports_1.1.5         
##  [27] assertthat_0.2.1         Matrix_1.2-17           
##  [29] lazyeval_0.2.2           cli_1.1.0               
##  [31] formatR_1.7              acepack_1.4.1           
##  [33] htmltools_0.4.0          prettyunits_1.0.2       
##  [35] tools_3.6.0              gtable_0.3.0            
##  [37] glue_1.3.1               GenomeInfoDbData_1.2.1  
##  [39] rappdirs_0.3.1           Rcpp_1.0.2              
##  [41] vctrs_0.2.0              Biostrings_2.53.2       
##  [43] nlme_3.1-141             rtracklayer_1.45.6      
##  [45] xfun_0.10                stringr_1.4.0           
##  [47] ps_1.3.0                 testthat_2.2.1          
##  [49] lifecycle_0.1.0          XML_3.98-1.20           
##  [51] zlibbioc_1.31.0          scales_1.0.0            
##  [53] hms_0.5.1                lambda.r_1.2.4          
##  [55] RColorBrewer_1.1-2       yaml_2.2.0              
##  [57] curl_4.2                 memoise_1.1.0           
##  [59] gridExtra_2.3            biomaRt_2.41.9          
##  [61] rpart_4.1-15             latticeExtra_0.6-28     
##  [63] stringi_1.4.3            RSQLite_2.1.2           
##  [65] genefilter_1.67.1        desc_1.2.0              
##  [67] checkmate_1.9.4          GenomicFeatures_1.37.4  
##  [69] pkgbuild_1.0.6           rlang_0.4.0             
##  [71] pkgconfig_2.0.3          bitops_1.0-6            
##  [73] evaluate_0.14            lattice_0.20-38         
##  [75] purrr_0.3.2              GenomicAlignments_1.21.7
##  [77] htmlwidgets_1.5.1        bit_1.1-14              
##  [79] processx_3.4.1           tidyselect_0.2.5        
##  [81] magrittr_1.5             R6_2.4.0                
##  [83] Hmisc_4.2-0              DBI_1.0.0               
##  [85] pillar_1.4.2             foreign_0.8-72          
##  [87] withr_2.1.2              mgcv_1.8-29             
##  [89] survival_2.44-1.1        RCurl_1.95-4.12         
##  [91] nnet_7.3-12              tibble_2.1.3            
##  [93] crayon_1.3.4             futile.options_1.0.1    
##  [95] utf8_1.1.4               BiocFileCache_1.9.1     
##  [97] rmarkdown_1.16           progress_1.2.2          
##  [99] locfit_1.5-9.1           data.table_1.12.2       
## [101] blob_1.2.0               callr_3.3.2             
## [103] digest_0.6.21            xtable_1.8-4            
## [105] openssl_1.4.1            munsell_0.5.0           
## [107] viridisLite_0.3.0        sessioninfo_1.1.1       
## [109] askpass_1.1
```

