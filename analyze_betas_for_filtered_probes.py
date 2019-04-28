'''
Calculate average, standard deviation, t-stat, and p-value for filtered
(top 100) probe IDs in BRCA vs. non-BRCA training examples.

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

import tensorflow as tf
import pandas as pd
import math
import glob
import sys
import os
import numpy as np
#import scipy
from scipy import stats

################################################################

filtered_probes_file = "top_100_probes.csv"

################################################################

study_with_label_1 = "BRCA"

list_of_disease_codes = [
        "BLCA", "BRCA", "CESC", "COAD", "DLBC", "ESCA", "GBM", "KICH", "KIRC",
        "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PRAD", "READ",
        "STAD", "TGCT", "THCA", "UCEC"
    ]

# Read probe_id indices, put into probe_indices dictionary
probe_indices = {} #maps probe_ids to indices
df = pd.read_csv("probe_id_indices")
for index, row in df.iterrows():
    probe_indices[row["probe_id"]] = row["index"] - 1

# Read indices of probes to include in the analysis list.
included_probes_indices = []
included_probes_ids = []
included_probes_betas_label0 = {}
included_probes_betas_label1 = {}
df = pd.read_csv(filtered_probes_file)
for index, row in df.iterrows():
    probe_id = row["probe_id"]
    probe_index = probe_indices[probe_id]
    included_probes_ids.append(probe_id)
    included_probes_indices.append(probe_index)
    included_probes_betas_label0[probe_index] = []
    included_probes_betas_label1[probe_index] = []

for disease_code in list_of_disease_codes:
    beta_files = sorted(glob.glob("horizontal_betas_" + disease_code + "_*"))
    for fname in beta_files:
        print "processing ", fname
        f = open(fname, 'r')
        for line in f:
            betas = line.split(",")
            study = betas.pop(0)
            aliquot = betas.pop(0)
            label = 1 if study == study_with_label_1 else 0
            for probe_idx in included_probes_indices:
                beta = float(betas[probe_idx])
                betas_list = included_probes_betas_label0[probe_idx] if label == 0 \
                             else included_probes_betas_label1[probe_idx]
                betas_list.append(beta)
        f.close()

for i in range(len(included_probes_ids)):
    probe_id = included_probes_ids[i]
    probe_index = included_probes_indices[i]
    betas_label0 = np.array(included_probes_betas_label0[probe_index])
    betas_label1 = np.array(included_probes_betas_label1[probe_index])
    t, p = stats.ttest_ind(betas_label0, betas_label1, equal_var=False)
    print "%s,%d,%d,%0.6f,%0.6f,%d,%0.6f,%0.6f,%0.2f,%0.2E" % (probe_id,
                                                               probe_index,
                                                               betas_label0.size,
                                                               np.mean(betas_label0),
                                                               np.std(betas_label0),
                                                               betas_label1.size,
                                                               np.mean(betas_label1),
                                                               np.std(betas_label1),
                                                               t,
                                                               p)
