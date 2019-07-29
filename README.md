# TCGA_workflow
useful scripts for TCGA data analysis

## Installation
```devtools::github("alexyfyf/TCGA_workflow")```

## data preparation
```
data <- data_download("TCGA_CHOL", datatype = "RNA-seq)
```

## plot survival plot for single gene
```
plot_surv_gene("TET2", test$x$genes, test$x$counts, test$surv, path = "./")
```

