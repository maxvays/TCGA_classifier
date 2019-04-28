'''
Generates files with methylation betas for each aliquot in one row,
as comma-separated values.

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
import os

rsrc = resource.RLIMIT_DATA
soft, hard = resource.getrlimit(rsrc)
print 'Soft limit starts as :', soft

resource.setrlimit(rsrc, (5*1024*1024*1024, hard))

soft, hard = resource.getrlimit(rsrc)
print 'Soft limit changed to :', soft

# Read probe_id indices, put into probe_indices dictionary
probe_indices = {} #maps probe_ids to indices
df = pd.read_csv("probe_id_indices")
for index, row in df.iterrows():
    probe_indices[row["probe_id"]] = row["index"] - 1

numprobes = len(probe_indices)

#list_of_disease_codes = ["BLCA", "BRCA", "CESC", "COAD", "DLBC", "ESCA", "GBM", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PRAD", "READ", "STAD", "TGCT", "THCA", "UCEC"]
list_of_disease_codes = sys.argv
list_of_disease_codes.pop(0)

for disease_code in list_of_disease_codes:
    # Get the list of files (shards) with betas
    listofbetafiles = sorted(glob.glob("sharded_betas_" + disease_code + "_*"))
    for filename in listofbetafiles:
        statinfo = os.stat(filename)
        if statinfo.st_size < 1024:
            continue
        out_filename = filename.replace("sharded", "horizontal")
        print "Processing:" + filename + ", output: " + out_filename
        out_file = open(out_filename, "w")
        aliquot_probe_dic = {} #maps aliquots to list of betas
        aliquot_beta_count = {} #maps aliquots to the number of betas inserted into aliquot_probe_dic
        num_done = 0
        df = pd.read_csv(filename)
        for df_index, row in df.iterrows():
            aliquot_barcode = row["aliquot_barcode"]
            probe_id = row["probe_id"]
            beta_value = row["beta_value"]
            probe_index = probe_indices[probe_id]
            if aliquot_barcode not in aliquot_probe_dic:
                listofbetas = [0 for i in range(numprobes)]
                aliquot_probe_dic[aliquot_barcode] = listofbetas
                aliquot_beta_count[aliquot_barcode] = 0
            listofbetas = aliquot_probe_dic[aliquot_barcode]
            listofbetas[probe_index] = beta_value
            aliquot_beta_count[aliquot_barcode] = aliquot_beta_count[aliquot_barcode] + 1
            if aliquot_beta_count[aliquot_barcode] == numprobes:
                num_done += 1
                print num_done, aliquot_barcode
                out_file.write(disease_code + "," + aliquot_barcode + "," + ','.join(map(str,listofbetas)) + "\n")
                del aliquot_probe_dic[aliquot_barcode]
        out_file.close()
