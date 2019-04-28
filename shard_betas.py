'''
Reads hg38_training_examples_* files (the sharded hg38_training_examples table
downloaded from the TCGA dataset in BigQuery) with the following columns:
disease_code, aliquot_barcode, probe_id, beta_value.
Also reads the list of aliquot barcodes for each disease code, and assigns
each aliquot to one of the sharded methylation betas files, so as to have 10
aliquots per shard.
Generates sharded methylation betas files, writing betas for each aliquot into
the assigned shard file.

Copyright (c) 2019 Maxim Vaysburd

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''
import pandas as pd
import glob
import sys
import math
import resource

rsrc = resource.RLIMIT_DATA
soft, hard = resource.getrlimit(rsrc)
print 'Soft limit starts as :', soft

resource.setrlimit(rsrc, (5*1024*1024*1024, hard))

soft, hard = resource.getrlimit(rsrc)
print 'Soft limit changed to :', soft

NUM_ALIQUOT_PER_BETA_FILE = 10

# Read aliquot barcode indices, put into aliquot_barcode_indices dictionary
max_beta_file_idx = 0

#list_of_disease_codes = ["BLCA", "BRCA", "CESC", "COAD", "DLBC", "ESCA", "GBM", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PRAD", "READ", "STAD", "TGCT", "THCA", "UCEC"]
# list_of_disease_codes = [
#     "BLCA", "CESC", "COAD", "DLBC", "ESCA", "GBM", "KICH",
#     "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO",
#     "PAAD", "PRAD", "READ", "STAD", "TGCT", "THCA", "UCEC"
# ]

list_of_disease_codes = sys.argv
list_of_disease_codes.pop(0)
print "list of disease codes: " + ",".join(list_of_disease_codes)

for disease_code in list_of_disease_codes:
    aliquot_barcode_indices = {} #maps aliquots to file indices
    print "processing: " + disease_code
    df = pd.read_csv("aliquot_barcode_" + disease_code)
    for index, row in df.iterrows():
        file_idx = row["idx"] / NUM_ALIQUOT_PER_BETA_FILE
        aliquot = row["aliquot_barcode"]
        aliquot_barcode_indices[aliquot] = file_idx
        max_beta_file_idx = max(max_beta_file_idx, file_idx)
        print aliquot + " --> file " + str(file_idx)

    # Create a file for each file index
    vertical_beta_files = []
    for i in range(max_beta_file_idx + 1):
        filepath = "sharded_betas_" + disease_code + "_%i" % i
        print "creating file " + filepath
        f = open(filepath, 'w')
        f.write("disease_code,aliquot_barcode,probe_id,beta_value\n")
        vertical_beta_files.append(f)


    # Get the list of files (shards) with betas
    listofbetafiles = sorted(glob.glob("hg38_training_examples_*" + disease_code + "*"))

    for filename in listofbetafiles:
        print "Processing:" + filename
        df = pd.read_csv(filename)
        for df_index, row in df.iterrows():
            dcode = row["disease_code"]
            if dcode != disease_code:
                continue
            aliquot_barcode = row["aliquot_barcode"]
            probe_id = row["probe_id"]
            beta_value = row["beta_value"]
            file_idx = aliquot_barcode_indices[aliquot_barcode]
            f = vertical_beta_files[file_idx]
            f.write(dcode + "," + aliquot_barcode + "," + probe_id + "," + str(beta_value) + "\n")
        for f in vertical_beta_files:
            f.flush()

    for f in vertical_beta_files:
        f.close()
