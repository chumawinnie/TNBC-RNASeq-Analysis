                                                     #fastqsanger.gz_trimmed.fq.gz trimmed.fq.gz
import os

##using ensembl genome through ftp for indices
#SRR14585768.fastqsanger.gz_trimmed.fq

file = open("Acc.txt")
accessions = file.read()
accessions = accessions.split()

#hisat2 indices_location -U ~/trimmed/xyz.fastq -S ~/mapping/xyz.sam ~/unaligned/xyz.sam 

inst_dir = "hisat2" #using hisat2 for mapping 
indices = "/home/chuma-winner/Rna-seq-BioProject/genome/r64/genome" #indices location, genome is the file name that is required by hisat2
option_SE = "-U" #unpaired data, it is single end
input_dir = "/home/chuma-winner/Rna-seq-BioProject/trimmed/" #where trimmed fastq files are stored
output_sam_dir = "/home/chuma-winner/Rna-seq-BioProject/mapping/sam/" #hisat2 sam output is stored here
output_sorted_sam = "/home/chuma-winner/Rna-seq-BioProject/mapping/sam_sorted/" #samtools sorted sam files are stored here
output_bam = "/home/chuma-winner/Rna-seq-BioProject/mapping/bam/" #samtools bam files are stored here 
output_unaligned = "/home/chuma-winner/Rna-seq-BioProject/mapping/unaligned/" #hisat2 unaligned sam files are stored here
reports_dir = "/home/chuma-winner/Rna-seq-BioProject/mapping/reports/" #alignment stats are stored here
option_stringtie = " --dta " #for stringtie tags
cores = " -p 4"


for x in accessions:
   
    #command = inst_dir + " " + "-x " + indices + " " + option_SE + " " + input_dir + x + ".fastqsanger.gz_trimmed.fq.gz" + " -S " + output_sam_dir + x + ".sam" + " " + "--un " + output_unaligned + x + "_non_aligned.sam" + " --summary " + reports_dir + x + "_summary.txt" + cores + option_stringtie
    #print(command + " ....running....")
    #os.system(command)
    #command = "samtools sort" + " -o " + output_sorted_sam + x + "_sorted.sam " + output_sam_dir + x + ".sam" #sorting the original sam files
    #print(command)
    #os.system(command)
    #sorted sam files are being coverted to bam files 
    command = "samtools view" + " -bS " + output_sorted_sam + x + "_sorted.sam" " > " + output_bam + x + ".bam"
    print(command)
    os.system(command)


print("All done")
