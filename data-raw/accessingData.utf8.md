---
title: "ImmSig2"
author: "EH"
date: "May 9, 2018"
output: html_document
---



## Accessing IS2 data


```r
library(ImmuneSpaceR) # Current master branch 5/9 ok 
library(Biobase)
library(dplyr)
library(UpdateAnno) # devtools::install_github("rglab/UpdateAnno")

# updateAnno is an rglab pkg with an easy method for updating annotation to latest
# from org.Hs.eg.db: http://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html 
```

Data other than gene expression is easy to get with our current API methods.


```r
con <- CreateConnection("IS2") # virtual study specific
```

```
## No gene expression data
```

```r
con # print to see list of available data sets
```

```
## <ImmuneSpaceConnection>
##   Study: IS2
##   URL: https://www.immunespace.org/HIPC/IS2
##   User: unknown_user at not_a_domain.com
##   Available Datasets:
##     - demographics
##     - fcs_control_files
##     - mbaa
##     - elisa
##     - cohort_membership
##     - fcs_sample_files
##     - neut_ab_titer
##     - fcs_analyzed_result
##     - elispot
##     - hai
##     - gene_expression_files
```

```r
hai <- con$getDataset("hai") # example with hai
head(hai)
```

```
##    participant_id age_reported gender  race cohort study_time_collected
## 1:  SUB112829.269           26   Male White  Young                    0
## 2:  SUB112829.269           26   Male White  Young                    0
## 3:  SUB112829.269           26   Male White  Young                    0
## 4:  SUB112829.269           26   Male White  Young                   28
## 5:  SUB112829.269           26   Male White  Young                   28
## 6:  SUB112829.269           26   Male White  Young                   28
##    study_time_collected_unit                  virus value_reported
## 1:                      Days A/South Dakota/06/2007             40
## 2:                      Days     A/Uruguay/716/2007             40
## 3:                      Days       B/Florida/4/2006             20
## 4:                      Days A/South Dakota/06/2007             40
## 5:                      Days     A/Uruguay/716/2007             40
## 6:                      Days       B/Florida/4/2006             40
##                                                                 lsid
## 1:  urn:lsid:labkey.com:Study.Data-198:5005.SUB112829.269.0.0000.445
## 2:  urn:lsid:labkey.com:Study.Data-198:5005.SUB112829.269.0.0000.501
## 3:  urn:lsid:labkey.com:Study.Data-198:5005.SUB112829.269.0.0000.557
## 4: urn:lsid:labkey.com:Study.Data-198:5005.SUB112829.269.28.0000.473
## 5: urn:lsid:labkey.com:Study.Data-198:5005.SUB112829.269.28.0000.529
## 6: urn:lsid:labkey.com:Study.Data-198:5005.SUB112829.269.28.0000.585
```

Getting gene expression data is a bit more complicated and we are working on creating a cleaner way of doing this within the study-specific connection.  For now, here is a workaround.


```r
# Use global connection to all studies for accessing select gene expression matrices
# SDY215 and SDY270 not included b/c lack GE data per Patrick @ ImmPort
gl <- CreateConnection("") 
mats <- gl$cache$GE_matrices 
studies <- c("SDY212, SDY180, SDY80, SDY61, SDY269, SDY522, SDY404, SDY400, SDY224, SDY67")
studies <- strsplit(studies, ", ")[[1]]
mats <- mats[ mats$folder %in% studies, ]

# Pull probe level data with annotation from when it was deposited with ImmPort
# and update annotation here as well as summarize by maximum probe value per 
# gene symbol.
mx <- lapply(mats$name, gl$getGEMatrix, 
             outputType = "normalized",
             annotation = "latest", 
             reload = TRUE)
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```

```
## Downloading matrix..
```

```
## Downloading Features..
```

```
## Constructing ExpressionSet
```


```r
mxEdit <- lapply(mx, function(x){
    # calc probe average without log transform
    exprs <- data.frame(Biobase::exprs(x))
    exprs$prb_avg <- apply(exprs, 1 , function(x){ sum(2^x) / length(x) }) 
    
    # update gene symbols - uses org.Hs.eg.db pkg, which pulls Entrez GeneIds periodically
    fdat <- Biobase::fData(x)
    exprs$gs <- UpdateAnno::updateAnno(fdat$gene_symbol)
    
    # filter to max probe for each gene symbol
    maxPrb <- exprs %>%
        group_by(gs) %>%
        filter(prb_avg == max(prb_avg)) %>%
        ungroup()
    
    # clean up
    maxPrb <- maxPrb[ !(duplicated(maxPrb)), ]
    maxPrb <- data.frame(maxPrb[ !is.na(maxPrb$gs), ])
    maxPrb <- maxPrb[ , colnames(maxPrb) != "prb_avg"]
})

# combine all the matrices into one
allMx <- Reduce(f = function(x, y){ merge(x, y, by = "gs")}, mxEdit)
rownames(allMx) <- allMx$gs
allMx <- allMx[ , colnames(allMx) != "gs" ]

# to map colnames to subject ids we use the full view from gene_expression_files dataset
# and then generate a partcipant_id + study_time_collected unique id.
gef <- gl$getDataset("gene_expression_files", original_view = T)
ptid <- gef$participant_id[ match(colnames(allMx), gef$biosample_accession)]
day <- gef$study_time_collected[ match(colnames(allMx), gef$biosample_accession) ]
colnames(allMx) <- paste0(ptid, "_d", day)
```

Subsetting the hai and gene expression datasets for subjects that only have both


```r
comboSubs <- intersect(ptid, hai$participant_id)
hai <- hai[ hai$participant_id %in% comboSubs, ]
exprs <- allMx[ , ptid %in% comboSubs ]
```

