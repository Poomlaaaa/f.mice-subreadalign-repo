#mRNA Analysis of Innate-like T Cells in Mice Treated with LiNKTcellsofαGalCer/CD1d-NP

#Introduction: This repository contains the code and instructions for performing a full mRNA analysis of innate-like T cells in mice treated with LiNKTcellsofαGalCer/CD1d-NP 
#and untreated controls. The analysis includes various visualization techniques such as MA plots, PCA, and volcano plots, among others. The primary focus is on using Subread 
#for alignment and subsequent analysis of the data.

#Build an mouse index
#For mice I downloaded 19:04.2024: Release M34 (GRCm39)
#https://www.gencodegenes.org/
chmod +x GRCm39.primary_assembly.genome.fa.gz
chmod +x gencode.vM34.primary_assembly.annotation.gtf.gz
subread-buildindex -o /archive/pumla/mouse/indexM/GRCm39 /archive/pumla/mouse/indexM/GRCm39.primary_assembly.genome.fa

#Obtain a fastq file
fastq-dump --stdout -X 2 SRR13893223

#To obtain fast.gz file of a paired end read
module avail
module load sratoolkit/3.0.2
fastq-dump  --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR13893223

#To obtain a fastqc report
module load fastqc
fastqc SRR13893223_pass_1.fastq.gz
fastqc SRR13893223_pass_2.fastq.gz

#To obtain alignment
module load subread
chmod +x subread_align.sh
bash subread_align.sh

#I would already have mu subread_align.sh file in filezilla written as:
#!/bin/bash
subread-align -i /archive/pumla/mouse/indexM/GRCm39 -r /archive/pumla/mouse/fastqM/SRR13893223_pass_1.fastq.gz -R /archive/pumla/mouse/fastqM/SRR13893223_pass_2.fastq.gz -t 0 --multiMapping -B 1 -T 16 -o /archive/pumla/mouse/outputsM/SRR13893223.bam --sortReadsByCoordinates 2>&1 | tee -a subread_results_23.txt

#To obtain my counts
#First create a file like below to keep in pn filezilla for example a file called `inkt`;
 
 library(Rsubread)
 library(BiocParallel)
 library(limma)
 library(edgeR)
 library(stringr)


#featurecount####
bams <- c("SRR13893218.bam",
          "SRR13893219.bam",
          "SRR13893220.bam",
          "SRR13893221.bam",
          "SRR13893222.bam",
          "SRR13893223.bam")

counts <- featureCounts(files = bams, annot.ext = "/archive/pumla/mouse/annoM/gencode.vM34.primary_assembly.annotation.gtf.gz", 
                        isGTFAnnotationFile= TRUE,
                        GTF.featureType= c("exon"),
                        GTF.attrType="gene_id",
                        useMetaFeatures=T,
                        countMultiMappingReads=T,
                        isPairedEnd=TRUE,
                        nthreads=24
)

counts_df <- as.data.frame(counts$counts)

write.table(counts$counts, "/archive/pumla/mouse/annoM/countsM/counts.tsv", quote = F, col.names = T, sep = "\t")

#calculate fpkm values
z <- DGEList(counts=counts$counts, genes=counts$annotation[,c("GeneID","Length")])
z<- calcNormFactors(z)
RPKM <-rpkm(z)

write.table(RPKM, "/archive/pumla/mouse/annoM/countsM/FPKM.tsv", quote = F, col.names = T, sep = "\t")

done

#To obtain my counts
cat inkt.r
Rscript inkt.r


