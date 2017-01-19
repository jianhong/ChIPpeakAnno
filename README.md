# ChipSeq

##Set up environment

Install the following softares firstly:

1. ngsplt(https://github.com/shenlab-sinai/ngsplot)

2. MACS(https://github.com/taoliu/MACS)

3. ChipSeq(https://github.com/aiminy/ChipSeq)

#Install ChipSeq

## install from github directly
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

# Download source, then install
```{r eval=FALSE}

#To download

git clone https://github.com/aiminy/ChipSeq

#To build 

R CMD build --resave-data ChipSeq

if knitr is not installed, please install it by performing the following step

R -e 'install.packages("knitr", dependencies = TRUE,repos = "http://cran.rstudio.org",lib="/nethome/axy148/R/lib64/R/library")'


#To Install
R CMD INSTALL ChipSeq_0.0.1.tar.gz -l ~/R/lib64/R/library/
```

##Peak Call and Annotation
```{r eval=FALSE}

input.file.dir="/projects/scratch/bbc/Project/Danny_chip/Alignment/BWA/"
output.file.dir="/scratch/projects/bbc/aiminy_project/
genome="Hs"
 
library(ChipSeq)
PeakCallAndAnnotation(input.file.dir,output.file.dir,genome)
```

##Bam files sorting, index and Visualization

+ Run ChipSeq to perform one-setp analysis for ChipSeq related data 

```{r eval=FALSE}
input.file.dir="/projects/scratch/bbc/Project/Danny_chip/Alignment/BWA/"
output.file.dir="/scratch/projects/bbc/aiminy_project/
genome="Hs"
 
BamFileSortIndexVisualization(input.file.dir,output.file.dir,genome)
 
```

## Perform the above steps in one streamlined procedure

You run this job on the linux cluster(pegasus)

```{bash eval=FALSE}


sh ~/R/lib64/R/library/ChipSeq/bin/BashRun.sh /projects/scratch/bbc/Project/Danny_chip/Alignment/BWA/ /scratch/projects/bbc/aiminy_project/ "Hs" 

#To check output 
ls -lrth /scratch/projects/bbc/aiminy_project/

You will find the follwing latest directories that was created with a timestamp

#For example,

#Peak call is located in  
/scratch/projects/bbc/aiminy_project/ReadBam_at_2016-12-01-18-58-57_PeakCall/
 
#Peak call is located in  
/scratch/projects/bbc/aiminy_project/ReadBam_at_2016-12-01-18-58-57_PeakAnnotation/

#Peak call is located in  
/scratch/projects/bbc/aiminy_project/ReadBam_at_2016-12-01-18-58-57_visualization/
```
