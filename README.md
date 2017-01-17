# ChipSeq

```{r eval=TRUE}
#In R console
# you need to have devtools
# it is ideal to have R installed in your directory:
# For example: 
#> .libPaths()
#[1] "/nethome/axy148/R/lib64/R/library"

library(devtools)
install_github("aiminy/ChipSeq",dependencies = T , force = T)

#If you use command line in pegasus terminal
R -e 'library(devtools);install_github("aiminy/ChipSeq")'

Rscript ~/ChipSeq/inst/bin/Run_Chip_Seq_interactive_model.r
```
