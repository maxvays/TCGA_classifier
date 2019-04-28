'''
Generates training examples files, with labels set to 1 for training examples
representing methylation betas for the target study (BRCA) and labels set to 0
for training examples representing methylation betas for all other studies.
Methylation betas for training examples are retrieved from horizontal-betas
files.

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

list_of_disease_codes = [
    "BLCA", "BRCA", "CESC", "COAD", "DLBC", "ESCA", "GBM", "KICH", "KIRC",
    "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PRAD", "READ",
    "STAD", "TGCT", "THCA", "UCEC"
]

EXAMPLES_PER_DISEASE_CODE = 750
NUM_ROWS_IN_HORIZONTAL_BETAS_FILE = 10
NUM_DISEASE_CODES = len(list_of_disease_codes)
DISEASE_CODE_FOR_TARGET_STUDY = "BRCA"
NUM_EXAMPLES_REQUIRED = EXAMPLES_PER_DISEASE_CODE / NUM_ROWS_IN_HORIZONTAL_BETAS_FILE
NUM_EXAMPLES_PER_FILE_PER_DISEASE_CODE = NUM_ROWS_IN_HORIZONTAL_BETAS_FILE
NUM_EXAMPLES_PER_FILE = NUM_DISEASE_CODES * NUM_EXAMPLES_PER_FILE_PER_DISEASE_CODE
NUM_TRAINING_FILES = EXAMPLES_PER_DISEASE_CODE / NUM_EXAMPLES_PER_FILE_PER_DISEASE_CODE

num_codes = len(list_of_disease_codes)

training_file_idx = 0

# map: disease code --> list of beta files
beta_files = {}

# map: disease code --> number of beta files
beta_files_counts = {}

# disease code --> index of beta file from which we are currently reading
disease_code_file_idx = {}

# initialize maps for each disease code
for disease_code in list_of_disease_codes:
    list_of_beta_files = sorted(glob.glob("horizontal_betas_" + disease_code + "_*"))
    beta_files[disease_code] = list_of_beta_files
    beta_files_counts[disease_code] = len(list_of_beta_files)
    disease_code_file_idx[disease_code] = 0

# generate each training file: use equal number of training examples from each
# disease code.
for training_file_idx in range(NUM_TRAINING_FILES):
    print "procssing training file " + str(training_file_idx)
    training_file = open("tcga_training_examples_" + str(training_file_idx), 'w')
    for disease_code in list_of_disease_codes:
        label = "1" if disease_code == DISEASE_CODE_FOR_TARGET_STUDY else "0"
        print "processing " + disease_code
        if beta_files_counts[disease_code] == 0:
            print "skipping " + disease_code + " (there are no horizonatal beta files)"
            continue
        filename = beta_files[disease_code][disease_code_file_idx[disease_code]]
        print "processing " + filename
        f = open(filename, 'r')
        # put entire contents of this horizontal betas file into the training examples file
        num_times = num_codes - 1 if disease_code == DISEASE_CODE_FOR_TARGET_STUDY else 1
        for line in f:
            for ntimes in range(num_times):
                training_file.write(label + "," + line)
        f.close()
        disease_code_file_idx[disease_code] = (disease_code_file_idx[disease_code] + 1) % beta_files_counts[disease_code]
    training_file.close()
