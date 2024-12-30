
import os

#GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcript --merge
##using ensembl genome through ftp for indices
#transcript assembly and abundance estimation

file = open("Acc.txt")
accessions = file.read()
accessions = accessions.split()


command = "stringtie -eB ~/_sorted.bam -G merged-transcript.gtf -o"

inst_dir = "stringtie"
reference = "/home/chuma-winner/Rna-seq-BioProject/stringtie/assemblies/merged-transcript.gtf" #"S.cerevease genome annotation" 
options = " -eB "
input_dir = "/home/chuma-winner/Rna-seq-BioProject/mapping/bam/"
output_dir = "/home/chuma-winner/Rna-seq-BioProject/ballgown/R-cvrg-tables/"
cores = " -p 4 "


for x in accessions:
   
    command = inst_dir + cores + " -G " + reference + " " + input_dir + x + ".bam" + " " + " -o " + output_dir + x + "/" + x + ".gtf" + options
    print(command + " ....running....")
    os.system(command)


print("All done")
