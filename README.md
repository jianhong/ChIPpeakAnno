# Intructions for using ChipSeq

##Set up environment

Install the following softares firstly:

1. ngsplt(https://github.com/shenlab-sinai/ngsplot)

2. MACS(https://github.com/taoliu/MACS)

# set Path

emacs .bashrc

export PYTHONPATH=/nethome/axy148/MACS-1.4.2-1/lib/python2.7/site-packages:$PYTHONPATH

export PATH=/nethome/axy148/MACS-1.4.2-1/bin:$PATH

export PATH=$HOME/NGS_tools/ngsplot/bin:$PATH

source .bashrc

#Install ChipSeq

## install from github directly
```{r or bash}
#In R console, you need to have devtools,it is ideal to have R installed in your directory:
#For example: 

> .libPaths()
[1] "/nethome/axy148/R/lib64/R/library"

#install test version
library(devtools)
install_github("aiminy/ChipSeq",ref='0.99.0',dependencies = T , force = T)

#install release version
library(devtools)
install_github("aiminy/ChipSeq",dependencies = T , force = T)



#If you use command line in pegasus terminal
R -e 'library(devtools);install_github("aiminy/ChipSeq",ref = '0.99.0',dependencies = T , force = T)'
```

## You can run ChipSeq using interactive model by follwing hints:

```{bash}
Rscript /nethome/axy148/R/lib64/R/library/ChipSeq/bin/Run_Chip_Seq_interactive_model.r
```
# Perform analysis in the streamlined batch model

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
