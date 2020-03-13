# TCGAworkflow
useful scripts for TCGA data analysis

## Installation
```devtools::github("alexyfyf/TCGAworkflow")```

## data preparation
```
data <- data_download("TCGA-CHOL", datatype = "RNA-seq")
```

## plot survival plot for single gene
```
plot_surv_gene("TET2", data$x$genes, data$x$counts, data$surv, path = "./")
```

