#' GetSampleInfo
#'
#' GetSampleInfo 
#' 
#' @param input.sample.file
#' @param input.bam.file 
#'  
#' 
#' @return
#' @export
#'
#' @examples
#' 
#' input.sample.file <- '/Volumes/Bioinformatics$/2017/DannyNewData/SampleID_INFO_ChIP_new_Danny.csv'
#' 
#' input.bam.file <- '/Volumes/Bioinformatics$/2017/DannyNewData/DannyNewData.txt'
#'
#' re <- GetSampleInfo(input.sample.file,input.bam.file)
#' 
#' 
GetSampleInfo <- function(input.sample.file, input.bam.file) {
    
    re <- read.csv(input.sample.file)
    re.c <- colnames(re)
    
    bam.file <- read.table(input.bam.file, header = FALSE)
    
    # print(bam.file)
    
    yy <- apply(bam.file, 1, function(u) {
        y <- basename(u)
        pos <- regexpr("\\.", y)
        pos <- pos - 1
        # cat(pos,'\n')
        y <- substr(y, 1, pos)
        y
    })
    
    bam.file.sample.name <- cbind(bam.file, yy)
    
    # print(bam.file.sample.name)
    colnames(bam.file.sample.name) <- c("file.name", "ID")
    
    Ab.type <- unique(as.character(re[, colnames(re) == re.c[3]]))
    Cell.type <- unique(as.character(re[, colnames(re) == re.c[2]]))
    
    re1 <- cbind(re[, 1:3], paste0(re[, 2], "_", re[, 3]))
    
    re11 <- merge(re1, bam.file.sample.name, by = "ID")
    
    colnames(re11)[4] = "Cell_TF"
    
    chipseqCell <- function(n, es) {
        value <- list(name = n, es = es)
        attr(value, "class") <- "ChipSeqCell"
        value
    }
    
    cellType <- unique(as.character(re11[re11$Type_TF == "Input", ]$Type_Cell))
    
    
    yy <- lapply(cellType, function(u, re11) {
        n <- nchar(u)
        # cat(u, n,'\n')
        
        index1 <- which(substr(re11$Type_Cell, 1, n) == u)
        z <- re11[index1, ]
        index2 <- which(nchar(as.character(z$Type_Cell)) <= n + 2)
        
        z2 <- z[index2, ]
        
        cc <- chipseqCell(u, z2)
        
        cc
    }, re11)
    
    
    # print(yy)
    
    
    colnames(re1)[4] = "Cell_TF"
    
    Cell.Ab.type <- unique(as.character(re1[, 4]))
    
    re2 <- lapply(1:length(Cell.Ab.type), function(u, Cell.Ab.type, re1) {
        
        x = Cell.Ab.type[u]
        z = re1
        
        ZZ <- as.character(z[which(z[, 4] == x), 1])
        ZZ
        
    }, Cell.Ab.type, re1)
    
    names(re2) <- Cell.Ab.type
    
    # Get the combination of cell type and antibody with replicates
    
    re21 <- re2[lapply(re2, length) == 1]
    
    re22 <- re2[lapply(re2, length) == 2]
    
    re3 <- list(re1 = re1, re2 = re2, re21 = re21, re22 = re22, re11 = re11, y = yy)
    
    return(re3)
    
}

#' peakcallwithinput
#' 
#' @export
#' @example 
#' 
#' genome='Hs'
#' re <- peakcallwithinput(input.sample.file,input.bam.file,genome,output.dir,peakcaller)
#' 
#' 
peakcallwithinput <- function(input.sample.file, input.bam.file, genome = c("Hs", 
    "hs", "HS", "hS"), output.dir, peakcaller = c("macs14", "macs2"), peakPvalue) {
    
    re <- GetSampleInfo(input.sample.file, input.bam.file)
    
    cellInfo <- re$y
    
    output.dir.name = dirname(input.sample.file)
    
    temp3 = file.path(output.dir.name, output.dir)
    
    if (!dir.exists(temp3)) {
        dir.create(temp3)
    }
    
    peakcaller <- match.arg(peakcaller)
    
    genome <- match.arg(genome)
    cmd10 <- paste("-f BAM", "-g", genome, "-n", sep = " ")
    
    switch(peakcaller, macs2 = {
        
        PATH1=Sys.getenv("PATH")
        
        macs2_Lib = file.path("/nethome/axy148/NGS_tools/MACS/bin/")
        
        Sys.setenv(PATH = paste0(macs2_Lib, ":", PATH1))
        
        cmd1 <- Sys.which("macs2")[[1]]
        
        cat(cmd1, "\n")
        
        cmd9 = paste(cmd1, "callpeak -t", sep = " ")
        # cmd10 = '-f BAM -g hs -n ' cmd11 = ' -m 6 18 --bw=200 -p 0.00001'
        #cmd11 <- paste("-m 6 18 --bw=200", "-p", peakPvalue, sep = " ")
        cmd11 <- paste("-p", peakPvalue, sep = " ")
    }, {
        cmd9 = "macs14 -t "
        # cmd10 = '-f BAM -g hs -n ' cmd11 = ' -m 6,18 --bw=200 -p 0.00001'
        
        #cmd11 <- paste("-m 6,18 --bw=200", "-p", peakPvalue, sep = " ")
        
        cmd11 <- paste("-p", peakPvalue, sep = " ")
    })
    
    # #re.out<-file.name.2 '/nethome/axy148/MACS/bin/macs2 callpeak -t
    # /scratch/projects/bbc/BAM2BW/'2016-03-23-2-18_S7_'_bam_sorted_by_position.bam
    # -c
    # /scratch/projects/bbc/Project/Lluis/Alignment/2016-03-20-1-1_S5__BWA/2016-03-20-1-1_S5_.bam
    # -f BAM -g mm -n '2016-03-23-2-18_S7_'_bam_mm --outdir
    # '/scratch/projects/bbc/Peak_chip_seq/Peak_chip_broad' -B -q 0.01' cmd9='macs14
    # -t ' cmd10='-f BAM -g hs -n ' cmd11=' -m 6,18 --bw=200 -p 0.00001'
    
    cellInfo.run <- lapply(1:length(cellInfo), function(u, cellInfo, temp3) {
        
        # x=cellInfo[[u]]
        
        x.name = cellInfo[[u]]$name
        
        
        # cat(x_name,'\n')
        
        es <- cellInfo[[u]]$es
        
        x.input <- es[es$Type_TF == "Input", ]$file.name
        
        x.sample <- es[es$Type_TF != "Input", ]
        
        x.run <- apply(x.sample, 1, function(x) {
            
            y <- x
            
            ID <- y[1]
            Type_Cell <- y[2]
            Type_TF <- y[3]
            Cell_TF <- y[4]
            file.name <- y[5]
            
            
            
            xx <- file.name
            xx.name = paste(ID, gsub(" ", "-", Type_Cell), Type_TF, sep = "-")
            
            cmd12 = paste(cmd9, xx, "-c", x.input, cmd10, file.path(temp3, paste0(xx.name, 
                "_hs_1.00e-05_", peakcaller)), cmd11, sep = " ")
            # 
            cmd12
            # print(cmd12)
        })
        
        # cat(x.name,'\n')
        
        # print(x.run)
        
        x.run
        
        # cmd12=paste(cmd9,x,'–c',x.input,cmd10,paste0(temp3,'/',paste0(x_name,'_hs_1.00e-05_macs142')),cmd11,sep='
        # ')
        
        # print(x.input)
        
        # system(cmd12, intern = TRUE, ignore.stderr = TRUE)
        
        # re=read.table(u,header=FALSE) re<-as.character(re[,1])
        # #colnames(re)=c('Count','GeneName') re
    }, cellInfo, temp3)
    
    
    names(cellInfo.run) = unlist(lapply(cellInfo, function(u) {
        u$name
    }))
    
    zzz <- unlist(cellInfo.run)
    
    lapply(1:length(zzz), function(u, zzz) {
        
        cat(as.character(zzz[u][[1]]), "\n")
        cat("\n")
        
        # system(as.character(zzz[u][[1]]), intern = TRUE, ignore.stderr = TRUE)
        
        system(as.character(zzz[u][[1]]))
        
    }, zzz)
    
    
    re <- list(cellInforun = cellInfo.run, zzz = zzz)
    
    # AnntationUsingChipSeeker(temp3,)
    # dirname('/Volumes/Bioinformatics$/2017/DannyNewData/PeakCall/')
    
    AnntationUsingChipSeeker(temp3, "peaks.bed", temp3, DD = 5000)
    
    return(re)
    
    # AnnotatePeak2(paste0(temp3,'/'),'*macs142_peaks.bed',7,paste0(output.dir.name,'PeakAnnotation_at_',temp2),genome='Hs')
    
    # AnnotatePeak3(paste0(temp3,'/'),paste0(output.dir.name,'_PeakAnnotation'),
    # genome='Hs')
    
    # BamFileSortIndexVisualization2(re,genome)
    
}





