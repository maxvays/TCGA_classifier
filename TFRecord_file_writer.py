'''
Reads training examples files (methylation betas and labels as comma-separated
values, one row per aliquot) and generates TensorFlow TFRecord files, containing
training examples (methylation betas and corresponding 0/1 label values) in the
format that can be fed directly into a TensorFlow classifier.

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

destination_dir = "BRCA"

try:
    os.stat(destination_dir)
except:
    os.mkdir(destination_dir)

list_of_training_examples_files = sorted(glob.glob("tcga_training_examples_*"))

for filename in list_of_training_examples_files:
    f = open(filename, 'r')
    tfr_file_name = destination_dir + "/" + filename + ".tfr"
    print "generating TFRecord file " + tfr_file_name + " from " + filename
    writer = tf.python_io.TFRecordWriter(tfr_file_name)
    for line in f:
        betas = line.split(",")
        label = int(betas.pop(0))
        disease_code = betas.pop(0)
        aliquot = betas.pop(0)
        betas = [float(x) for x in betas]
        features = {
          "betas": tf.train.Feature(float_list = tf.train.FloatList(value=betas)),
          "label": tf.train.Feature(int64_list = tf.train.Int64List(value=[label])),
        }
        example = tf.train.Example(
            features=tf.train.Features(feature=features))
        writer.write(example.SerializeToString())
    writer.close()
