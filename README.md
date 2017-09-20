# Intructions for using ChipSeq

##Set up environment

Install the following softares firstly:

1. ngsplt(https://github.com/shenlab-sinai/ngsplot)

2. MACS(https://github.com/taoliu/MACS)

3. sudo apt-get install libboost-dev

4. R -e 'library(devtools);devtools::install_github("hms-dbmi/spp", build_vignettes = FALSE)'

# set Path

emacs .bashrc

export PYTHONPATH=/nethome/axy148/MACS-1.4.2-1/lib/python2.7/site-packages:$PYTHONPATH

export PATH=/nethome/axy148/MACS-1.4.2-1/bin:$PATH

export PATH=$HOME/NGS_tools/ngsplot/bin:$PATH

source .bashrc

#Install ChipSeq


# install
```{r}
R -e 'library(devtools);install_github("aiminy/ChipSeq")'

R -e 'library(ChipSeq);library(ChipSeq);ChipSeq:::parseToSampleInfo("/projects/scratch/bbc/Project/Danny_chip3/Filtered_bam","*.bam$","~/Danny_chip3","sample_infor_Danny_chip3.txt","Danny_chip3")'

# cluster
R -e 'library(PathwaySplice);library(DoGs);library(ChipSeq);ChipSeq:::submitJob4useRunSppR(""~/Danny_chip3/sample_infor_Danny_chip3.txt","/projects/scratch/bbc/Project/Danny_chip3/Filtered_bam","/scratch/projects/bbc/aiminy_project/Danny_chip3_Chipseq_QC")'

# local
R -e 'library(ChipSeq);re <- ChipSeq:::useRunSppR("~/Danny_chip3/sample_infor_Danny_chip3.txt","/media/aiminyan/DATA/Danny_chip3","/media/aiminyan/DATA/Danny_chip3_chipSeq_QC")'

```

