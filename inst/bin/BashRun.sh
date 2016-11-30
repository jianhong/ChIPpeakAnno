#!/bin/bash

#Usage: sh ~/Code/BashRunMACS1-4-2_4_Danny_chip_seq.sh input_file_dir out_file_dir

#Example: sh ~/Code/BashRunMACS1-4-2_4_Danny_chip_seq3.sh /scratch/projects/bbc/Project/Danny_chip/Alignment/BWA/ 9 ".bam" /scratch/projects/bbc/aiminy_project/Bam_marked_sorted/

#Example: sh ~/Code/BashRunMACS1-4-2_4_Danny_chip_seq3.sh /scratch/projects/bbc/aiminy_project/Bam_sorted/ 7 ".sorted.bam$" /scratch/projects/bbc/aiminy_project/Bam_marked_sorted/ test.txt  

#Example: sh ~/Code/BashRunMACS1-4-2_4_Danny_chip_seq3.sh /projects/scratch/bbc/Project/Danny_chip/Alignment/BWA/ 9 ".bam" /scratch/projects/bbc/aiminy_project/ "Hs"

#Example: sh ~/Code/BashRunMACS1-4-2_4_Danny_chip_seq3.sh /scratch/projects/bbc/aiminy_project/test_bam/ 7 ".bam" /scratch/projects/bbc/aiminy_project/ "Hs"

#Example: sh ~R/lib64/R/library/ChipSeq/bin/BashRun.sh /scratch/projects/bbc/aiminy_project/test_bam/ 7 ".bam" /scratch/projects/bbc/aiminy_project/ "Hs" 

#Example  sh ~/R/lib64/R/library/ChipSeq/bin/BashRun.sh /scratch/projects/bbc/aiminy_project/test_bam/ 7 ".bam" /scratch/projects/bbc/aiminy_project/ "Hs" 


DIR="$1"
INDEX="$2"
pattern="$3"
results_dir="$4"
out="$5"


if [ -e $results_dir ]; then echo "already exists";
else
    mkdir $results_dir

fi


cat >  ~/Code/code_tmp/Run_Rscript.sh <<EOF

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

Rscript R_code/Set_up_R.R
wait

Rscript R_code/Run_Chip_Seq.R $DIR $INDEX $pattern $results_dir $out 

EOF

bsub -P bbc < ~/Code/code_tmp/Run_Rscript.sh
