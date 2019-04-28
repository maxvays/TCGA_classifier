'''
Sorts probe IDs by the sums of the absolute values of weights in the
classifier's hidden layer (dnn/hiddenlayer_0/kernel), and outputs
(probe_index, probe_id, sum_weights) for the top 1000 probe IDs.

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
import getopt

print "reading probe indices"
probe_indices = {} #maps probe indices to probe_ids
df = pd.read_csv("probe_id_indices")
for index, row in df.iterrows():
    probe_indices[row["index"]] = row["probe_id"]

rows_list = []

print "reading feature weights"
weights = pd.read_csv("weights_dnn_hiddenlayer_0_kernel.csv")
for index, row in weights.iterrows():
    probe_id = probe_indices[index+1]
    values = row.tolist()
    #print values, "\n"
    sum = 0
    for i in range(len(values)-1):
        sum += abs(float(values[i]))
    #print str(index+1) + "," + probe_id + "," + str(sum)
    rows_list.append([index+1, probe_id, sum])
    #sys.exit(0)

summary_df = pd.DataFrame(rows_list, columns=["probe_index", "probe_id", "sum_weights"])
summary_df = summary_df.sort_values(by=['sum_weights'], ascending=False)

top_1000 = summary_df[:1000]

for index, row in top_1000.iterrows():
    values = row.tolist()
    print ",".join([str(x) for x in values])
