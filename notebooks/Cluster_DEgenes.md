Cluster DE genes R Notebook
================

``` r
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggsci)
library(openxlsx)
options(future.globals.maxSize = 4000 * 1024^2)
```

# read data

``` r
seurat <- readRDS(file = "RDSfiles/hand2.bud.clustered.RDS")
Idents(seurat) <- "sub.cluster"
DimPlot(seurat) + scale_color_igv()
```

![](Cluster_DEgenes_files/figure-gfm/readdata-1.png)<!-- -->

``` r
levels(seurat)
```

    ##  [1] "endothelial and hematopoietic precursors"
    ##  [2] "cranial mesoderm"                        
    ##  [3] "ribosomal"                               
    ##  [4] "myocardium"                              
    ##  [5] "pericardium"                             
    ##  [6] "mesothelium"                             
    ##  [7] "PLPM"                                    
    ##  [8] "posterior hemangioblasts & kidney"       
    ##  [9] "paraxial mesoderm"                       
    ## [10] "uncommitted PLPM"                        
    ## [11] "pan-neuronal"                            
    ## [12] "neural plate"                            
    ## [13] "epidermal"                               
    ## [14] "endoderm 1"                              
    ## [15] "endoderm 2"                              
    ## [16] "notochord"                               
    ## [17] "prechordal plate & hatching gland 1"     
    ## [18] "prechordal plate & hatching gland 2"

``` r
abbreviate(gsub(" ",".",gsub("_",".",levels(seurat))), minlength = 20)
```

    ## endothelial.and.hematopoietic.precursors 
    ##                   "endthll.nd.hmtptc.pr" 
    ##                         cranial.mesoderm 
    ##                       "cranial.mesoderm" 
    ##                                ribosomal 
    ##                              "ribosomal" 
    ##                               myocardium 
    ##                             "myocardium" 
    ##                              pericardium 
    ##                            "pericardium" 
    ##                              mesothelium 
    ##                            "mesothelium" 
    ##                                     PLPM 
    ##                                   "PLPM" 
    ##        posterior.hemangioblasts.&.kidney 
    ##                   "pstrr.hmngblsts.&.kd" 
    ##                        paraxial.mesoderm 
    ##                      "paraxial.mesoderm" 
    ##                         uncommitted.PLPM 
    ##                       "uncommitted.PLPM" 
    ##                             pan-neuronal 
    ##                           "pan-neuronal" 
    ##                             neural.plate 
    ##                           "neural.plate" 
    ##                                epidermal 
    ##                              "epidermal" 
    ##                               endoderm.1 
    ##                             "endoderm.1" 
    ##                               endoderm.2 
    ##                             "endoderm.2" 
    ##                                notochord 
    ##                              "notochord" 
    ##      prechordal.plate.&.hatching.gland.1 
    ##                   "prchrdl.plt.&.htc..1" 
    ##      prechordal.plate.&.hatching.gland.2 
    ##                   "prchrdl.plt.&.htc..2"

# find all DE genes

``` r
markers <- FindAllMarkers(seurat, only.pos = T, verbose = F)
head(markers)
```

    ##                          p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## si:ch73-364h19.1 1.690412e-303   4.854271 0.970 0.084 6.144986e-299
    ## spi1b            1.212914e-268   6.601440 0.795 0.050 4.409183e-264
    ## hhex             2.236793e-266   4.651133 0.795 0.046 8.131189e-262
    ## sox7             2.632554e-256   4.942734 0.951 0.120 9.569859e-252
    ## tal1             9.000069e-231   4.293394 0.996 0.186 3.271705e-226
    ## egfl7            2.060682e-222   4.983454 0.750 0.062 7.490990e-218
    ##                                                   cluster             gene
    ## si:ch73-364h19.1 endothelial and hematopoietic precursors si:ch73-364h19.1
    ## spi1b            endothelial and hematopoietic precursors            spi1b
    ## hhex             endothelial and hematopoietic precursors             hhex
    ## sox7             endothelial and hematopoietic precursors             sox7
    ## tal1             endothelial and hematopoietic precursors             tal1
    ## egfl7            endothelial and hematopoietic precursors            egfl7

``` r
if (!dir.exists("results")){
  dir.create("results")
}
DEgenelist <- list()
for(cluster in levels(seurat)){
  clustername <- abbreviate(gsub(" ",".",gsub("_",".",cluster)), minlength = 20)
  DEgenelist[[clustername]] <- markers[markers$cluster == cluster,]
}
write.xlsx(DEgenelist, file = "results/hand2_bud_cluster_DEgenes.xlsx")
```

``` r
sessionInfo()
```

    ## R version 4.3.0 (2023-04-21)
    ## Platform: x86_64-apple-darwin20 (64-bit)
    ## Running under: macOS Monterey 12.6.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Denver
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] openxlsx_4.2.5.2   ggsci_3.0.3        dplyr_1.1.4        ggplot2_3.5.0     
    ## [5] Seurat_5.0.3       SeuratObject_5.0.1 sp_2.1-3          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3         
    ##   [4] rlang_1.1.3            magrittr_2.0.3         RcppAnnoy_0.0.22      
    ##   [7] spatstat.geom_3.2-9    matrixStats_1.3.0      ggridges_0.5.6        
    ##  [10] compiler_4.3.0         png_0.1-8              vctrs_0.6.5           
    ##  [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3       
    ##  [16] fastmap_1.2.0          labeling_0.4.3         utf8_1.2.4            
    ##  [19] promises_1.3.0         rmarkdown_2.26         purrr_1.0.2           
    ##  [22] xfun_0.43              jsonlite_1.8.8         goftest_1.2-3         
    ##  [25] highr_0.10             later_1.3.2            spatstat.utils_3.0-4  
    ##  [28] irlba_2.3.5.1          parallel_4.3.0         cluster_2.1.6         
    ##  [31] R6_2.5.1               ica_1.0-3              stringi_1.8.3         
    ##  [34] RColorBrewer_1.1-3     spatstat.data_3.0-4    limma_3.58.1          
    ##  [37] reticulate_1.35.0      parallelly_1.37.1      lmtest_0.9-40         
    ##  [40] scattermore_1.2        Rcpp_1.0.12            knitr_1.45            
    ##  [43] tensor_1.5             future.apply_1.11.2    zoo_1.8-12            
    ##  [46] sctransform_0.4.1      httpuv_1.6.15          Matrix_1.6-5          
    ##  [49] splines_4.3.0          igraph_2.0.3           tidyselect_1.2.1      
    ##  [52] abind_1.4-5            rstudioapi_0.16.0      yaml_2.3.8            
    ##  [55] spatstat.random_3.2-3  codetools_0.2-20       miniUI_0.1.1.1        
    ##  [58] spatstat.explore_3.2-7 listenv_0.9.1          lattice_0.22-6        
    ##  [61] tibble_3.2.1           plyr_1.8.9             withr_3.0.0           
    ##  [64] shiny_1.8.1.1          ROCR_1.0-11            evaluate_0.24.0       
    ##  [67] Rtsne_0.17             future_1.33.2          fastDummies_1.7.3     
    ##  [70] survival_3.5-8         polyclip_1.10-6        zip_2.3.1             
    ##  [73] fitdistrplus_1.1-11    pillar_1.9.0           KernSmooth_2.23-24    
    ##  [76] plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0        
    ##  [79] munsell_0.5.1          scales_1.3.0           globals_0.16.3        
    ##  [82] xtable_1.8-4           glue_1.7.0             lazyeval_0.2.2        
    ##  [85] tools_4.3.0            data.table_1.15.4      RSpectra_0.16-1       
    ##  [88] RANN_2.6.1             leiden_0.4.3.1         dotCall64_1.1-1       
    ##  [91] cowplot_1.1.3          grid_4.3.0             tidyr_1.3.1           
    ##  [94] colorspace_2.1-0       nlme_3.1-165           patchwork_1.2.0       
    ##  [97] presto_1.0.0           cli_3.6.3              spatstat.sparse_3.0-3 
    ## [100] spam_2.10-0            fansi_1.0.6            viridisLite_0.4.2     
    ## [103] uwot_0.1.16            gtable_0.3.4           digest_0.6.36         
    ## [106] progressr_0.14.0       ggrepel_0.9.5          farver_2.1.2          
    ## [109] htmlwidgets_1.6.4      htmltools_0.5.8.1      lifecycle_1.0.4       
    ## [112] httr_1.4.7             statmod_1.5.0          mime_0.12             
    ## [115] MASS_7.3-60.0.1
