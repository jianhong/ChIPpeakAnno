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

## Install ChipSeq
```{r}
# install ChipSeq
R -e 'library(devtools);install_github("aiminy/ChipSeq")'
```

## Generate sample_infor_Danny_chip3.txt in ~/Danny_chip3
```{r} 
R -e 'library(ChipSeq);library(ChipSeq);ChipSeq:::parseToSampleInfo("/projects/scratch/bbc/Project/Danny_chip3/Filtered_bam","*.bam$","~/Danny_chip3","sample_infor_Danny_chip3.txt","Danny_chip3")'
```  

## Get CC plot on cluster

### Get sample_infor_Danny_chip3.txt from the following link
https://www.dropbox.com/sh/z20u8ofiub8v7pk/AAAEYZwxHj9wUfZJxqf7QsNMa?dl=0
```{r}
R -e 'library(PathwaySplice);library(DoGs);library(ChipSeq);ChipSeq:::submitJob4useRunSppR(""~/Danny_chip3/sample_infor_Danny_chip3.txt","/projects/scratch/bbc/Project/Danny_chip3/Filtered_bam","/scratch/projects/bbc/aiminy_project/Danny_chip3_Chipseq_QC")'
```

## Get CC plot locally
```{r}
R -e 'library(ChipSeq);re <- ChipSeq:::useRunSppR("~/Danny_chip3/sample_infor_Danny_chip3.txt","/media/aiminyan/DATA/Danny_chip3","/media/aiminyan/DATA/Danny_chip3_chipSeq_QC")'
```

## Get heatmap locally

From https://drive.google.com/drive/folders/0B1PVLadG_dCKNEsybkh5TE9XZ1E, download the genome data file

```{r}
# For example, to install hg19 genome data, you need to download ngsplotdb_hg19_75_3.00.tar.gz,ngsplotdb_hg19_75_3.00_enhancer.tar.gz and ngsplotdb_hg19_75_3.00_dhs.tar.gz to ~/Download firstly

# Then install these data using the following commands:
ngsplotdb.py install ~/Downloads/ngsplotdb_hg19_75_3.00.tar.gz
ngsplotdb.py install ~/Downloads/ngsplotdb_hg19_75_3.00_enhancer.tar.gz
ngsplotdb.py install ~/Downloads/ngsplotdb_hg19_75_3.00_dhs.tar.gz
```

## If bam files are not sorted and indexed yet, run the foollowing with bam.sort=FALSE
```{r}
R -e 'library(ChipSeq);re <- ChipSeq:::BamFileSortIndexVisualization("/media/aiminyan/DATA/Danny_chip3","*.bam$","/media/aiminyan/DATA/Danny_chip3_heatmap",5000,"Hs")'
```  

## If bam files are already sorted and indexed, run the following with bam.sort=TRUE
```{r}
R -e 'library(ChipSeq);re <- ChipSeq:::BamFileSortIndexVisualization("/media/aiminyan/DATA/Danny_chip3","*.bam$","/media/aiminyan/DATA/Danny_chip3_heatmap",5000,"Hs",bam.sort=TRUE)'
```  

## Generate heatmap normalized by inputs, you need to make a Danny_chip3_4_ngs_plot.txt firstly
```{r}
R -e 'library(ChipSeq);ChipSeq:::ngs2("/media/H_driver/Aimin_project/Danny_chip3_4_ngs_plot.txt",5000,"/media/H_driver/Aimin_project/heapmapNormalizedByInput")'
```

## Generate Normalized-By-Subtract-Inputs-Rm-BlackList-Region BigWig files
```{r}
R -e 'library(ChipSeq);ChipSeq:::getBwUseBamCompare("/media/H_driver/Aimin_project/Danny_chip3.txt","/media/H_driver/Aimin_project/REF/consensusBlacklist.bed","/media/H_driver/Aimin_project/Norma_sub_Coverage_rm_bl")'
```

## Call peaks use macs2 without using input
```{r} 
 R -e 'library(ChipSeq);re <- ChipSeq:::peakCallAndAnnotationWithoutInput("/media/aiminyan/DATA/Danny_chip3","/media/H_driver/Danny_chip3_macs2_call","hs","macs2",0.0001)'
```

## Annotate peaks
```{r}
 R -e 'library(ChipSeq);re <- ChipSeq:::AnotationUsingChipSeeker3("/media/H_driver/Aimin_project/Danny_chip3_macs2_call/PeakCall","*.narrowPeak$","/media/H_driver/Aimin_project/Danny_chip3_macs2_annotation",txdb="hg19",DD=5000,distanceToTSS_cutoff=10000)'
```

## Call peaks use macs2 with input
```{r}
R -e 'library(ChipSeq);ChipSeq:::peakCall2("/media/H_driver/Aimin_project/Danny_chip3.txt","hs","/media/H_driver/Aimin_project/peak_call_with_input","macs2",0.00001)'
```

## Annotate peaks called with inputs
```{r}
R -e 'library(ChipSeq);re <- ChipSeq:::AnotationUsingChipSeeker3("/media/H_driver/Aimin_project/peak_call_with_input","*.narrowPeak$","/media/H_driver/Aimin_project/Danny_chip3_macs2_annotation_with_ input",txdb="hg19",DD=5000,distanceToTSS_cutoff=10000)'
```

## Get sequence around summit peak
```{r}
R -e 'library(ChipSeq);re <- ChipSeq:::getSummitSequence("/media/H_driver/Aimin_project/peak_call_with_input","*.xls$","hg19","/media/H_driver/Aimin_project/SummitSeq")'
```
