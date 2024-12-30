import os

import glob

# All files ending with .gz/fastq
files = glob.glob("/home/chuma-winner/Rna-seq-BioProject/*.gz")


inst_dir = "fastqc"
input_dir = "/home/chuma-winner/Rna-seq-BioProject/"
output_dir = "/home/chuma-winner/Rna-seq-BioProject/fastqc_output/"

for x in files:

    command = inst_dir + " " + x + " -o " + output_dir
    print(x + " ....running....")
    os.system(command)

print("All done")
