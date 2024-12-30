import os


file = open("Acc.txt")  #reading sample accession to open files from directory for processing
accessions = file.read()
accessions = accessions.split() # a list of accessions is created by removing special characters !



command = "trim_galore in1.fastq in2.fastq -o out1.fastq -p out2.fastq"

inst_dir = "~/TrimGalore-0.6.7/trim_galore" #downlaoded from: https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.7.tar.gz
option_fastqc = "--fastqc" #to use fastqc after trimming
input_dir = "/home/chuma-winner/Rna-seq-BioProject/" #where the raw files are
output_dir =  "/home/chuma-winner/Rna-seq-BioProject/trimmed/" #where trimmed files will go

for x in accessions: #a list of accesions looping through the variables !
   
    command = inst_dir + " " + input_dir + x + ".fastqsanger.gz" + " " + " -o " + output_dir  + " " + option_fastqc + " -j 4" 
    print(command + " ....running....")
    os.system(command)

print("All done")
