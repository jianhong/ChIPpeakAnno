#!/bin/bash

#Usage: sh ~/R/lib64/R/library/ChipSeq/bin/BashRun.sh input_file_dir out_file_dir

#Example: sh ~/R/lib64/R/library/ChipSeq/bin/BashRun.sh /scratch/projects/bbc/Project/Danny_chip/Alignment/BWA/ 9 ".bam" /scratch/projects/bbc/aiminy_project/Bam_marked_sorted/

#Example: sh ~/R/lib64/R/library/ChipSeq/bin/BashRun.sh /scratch/projects/bbc/aiminy_project/Bam_sorted/ 7 ".sorted.bam$" /scratch/projects/bbc/aiminy_project/Bam_marked_sorted/ test.txt  

#Example: sh ~/R/lib64/R/library/ChipSeq/bin/BashRun.sh /projects/scratch/bbc/Project/Danny_chip/Alignment/BWA/ 9 ".bam" /scratch/projects/bbc/aiminy_project/ "Hs"

#Example: sh ~/R/lib64/R/library/ChipSeq/bin/BashRun.sh /scratch/projects/bbc/aiminy_project/test_bam/ 7 ".bam" /scratch/projects/bbc/aiminy_project/ "Hs"

#Example: sh ~/R/lib64/R/library/ChipSeq/bin/BashRun.sh /scratch/projects/bbc/aiminy_project/test_bam/ 7 ".bam" /scratch/projects/bbc/aiminy_project/ "Hs" 

#Example  sh ~/R/lib64/R/library/ChipSeq/bin/BashRun.sh /scratch/projects/bbc/aiminy_project/test_bam/ 7 ".bam" /scratch/projects/bbc/aiminy_project/ "Hs" 

#Example  sh ~/R/lib64/R/library/ChipSeq/bin/BashRun.sh /scratch/projects/bbc/aiminy_project/test_bam/ /scratch/projects/bbc/aiminy_project/ "Hs" 

DIR="$1"
#INDEX="$2"
#pattern="$3"
results_dir="$2"
out="$3"

if [ -e $results_dir ]; then echo "already exists";
else
    mkdir $results_dir

fi

cat >  ~/R/lib64/R/library/ChipSeq/bin/Run_Rscript.sh <<EOF

#LSBATCH: User input                                                                                                                                                                 
#!/bin/bash                                                                                                                                                                          
#BSUB -P bbc                                                                                                                                                                         
#BSUB -J RunR                                                                                                                                                                
#BSUB -o %J.RunR.log                                                                                                                                                             
#BSUB -e %J.RunR.err                                                                                                                                                             
#BSUB -W 72:00                                                                                                                                                                       
#BSUB -n 8                                                                                                                                                    
#BSUB -q general
#BSUB -u aimin.yan@med.miami.edu              

Rscript ~/R/lib64/R/library/ChipSeq/bin/Set_up_R.r
wait

Rscript ~/R/lib64/R/library/ChipSeq/bin/Run_Chip_Seq.r $DIR $results_dir $out 

EOF

bsub -P bbc < ~/R/lib64/R/library/ChipSeq/bin/Run_Rscript.sh