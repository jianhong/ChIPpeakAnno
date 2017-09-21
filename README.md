# Intructions for using ChipSeq

## Set up environment

Install the following softares firstly:

1. ngsplt(https://github.com/shenlab-sinai/ngsplot)

2. MACS(https://github.com/taoliu/MACS)

3. sudo apt-get install libboost-dev

4. R -e 'library(devtools);devtools::install_github("hms-dbmi/spp", build_vignettes = FALSE)'

## Set Path

emacs .bashrc

export PYTHONPATH=/nethome/axy148/MACS-1.4.2-1/lib/python2.7/site-packages:$PYTHONPATH

export PATH=/nethome/axy148/MACS-1.4.2-1/bin:$PATH

export PATH=$HOME/NGS_tools/ngsplot/bin:$PATH

source .bashrc

## Use ChipSeq

```{r}
# install ChipSeq
# 
R -e 'library(devtools);install_github("aiminy/ChipSeq")'

# Generate sample_infor_Danny_chip3.txt in ~/Danny_chip3
# 
R -e 'library(ChipSeq);library(ChipSeq);ChipSeq:::parseToSampleInfo("/projects/scratch/bbc/Project/Danny_chip3/Filtered_bam","*.bam$","~/Danny_chip3","sample_infor_Danny_chip3.txt","Danny_chip3")'

# Get CC plot on cluster
# 
R -e 'library(PathwaySplice);library(DoGs);library(ChipSeq);ChipSeq:::submitJob4useRunSppR(""~/Danny_chip3/sample_infor_Danny_chip3.txt","/projects/scratch/bbc/Project/Danny_chip3/Filtered_bam","/scratch/projects/bbc/aiminy_project/Danny_chip3_Chipseq_QC")'

# Get CC plot locally
# 
R -e 'library(ChipSeq);re <- ChipSeq:::useRunSppR("~/Danny_chip3/sample_infor_Danny_chip3.txt","/media/aiminyan/DATA/Danny_chip3","/media/aiminyan/DATA/Danny_chip3_chipSeq_QC")'

# Get heatmap locally

#install hg19 genome data if this data is not in your environment, 
# download ngsplotdb_hg19_75_3.00.tar.gz,ngsplotdb_hg19_75_3.00_enhancer.tar.gz and ngsplotdb_hg19_75_3.00_dhs.tar.gz to ~/Downloads
# install these data
ngsplotdb.py install ~/Downloads/ngsplotdb_hg19_75_3.00.tar.gz
ngsplotdb.py install ~/Downloads/ngsplotdb_hg19_75_3.00_enhancer.tar.gz
ngsplotdb.py install ~/Downloads/ngsplotdb_hg19_75_3.00_dhs.tar.gz

# if bam files are not sorted and indexed yet, run the foollowing 
R -e 'library(ChipSeq);re <- ChipSeq:::BamFileSortIndexVisualization("/media/aiminyan/DATA/Danny_chip3","*.bam$","/media/aiminyan/DATA/Danny_chip3_heatmap",5000,"Hs")'

# if bam files are already sorted and indexed, run the following 
R -e 'library(ChipSeq);re <- ChipSeq:::BamFileSortIndexVisualization("/media/aiminyan/DATA/Danny_chip3","*.bam$","/media/aiminyan/DATA/Danny_chip3_heatmap",5000,"Hs",bam.sort=TRUE)'

```
