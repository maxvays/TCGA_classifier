'''
Trains TensorFlow DNNClassifier on training examples in TFRecord format, with
each training example comprising methylation betas for one aliquot and the 0/1
label.
Outputs statistics on the number of input features with non-zero weights in the
hidden layer (dnn/hiddenlayer_0/kernel).
Upon completion of training, writes the weights for all input features in the
hidden layer (dnn/hiddenlayer_0/kernel) into a separate file for analysis.

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
import glob
import sys
import math
import resource
import getopt
from tensorflow.python.ops import variables
from tensorflow.python.framework import dtypes

opts, args = getopt.getopt(
	sys.argv[1:], "",
	["regularization=", "eval_start_idx=", "study=", "tag=", "learning_rate="])

regularization = None
eval_start_idx = None
study = None
tag = None
learning_rate = None

NUM_TRAINING_FILES = 75
NUM_EVAL_FILES = 15

for o, a in opts:
	if o == "--regularization":
		print "setting regularization = " + a
		regularization = float(a)
	elif o == "--eval_start_idx":
		print "setting eval_start_idx = " + a
		eval_start_idx = int(a)
	elif o == "--study":
		print "setting study = " + a
		study = a
	elif o == "--tag":
		print "setting tag = " + a
		tag = a
        elif o == "--learning_rate":
                print "setting learning rate = " + a
                learning_rate = float(a)

if regularization == None or eval_start_idx == None or study == None or tag == None or learning_rate == None:
	print("specify --regularization, --eval_start_idx, --study, --tag, --learning_rate")
	sys.exit(1)

modeldir = "model_{tag}_{study}_e{eval_start_idx}_lr{learning_rate}_r{regularization}".format(
	tag=tag,
	study=study,
	eval_start_idx=eval_start_idx,
  	regularization=regularization,
        learning_rate=learning_rate
)

print "model dir: " + modeldir

TRAINING_FILE_NAME = "{study}/tcga_training_examples_{idx}.tfr"

training_files = []
eval_files = []

for i in range(NUM_TRAINING_FILES):
	filename = TRAINING_FILE_NAME.format(study=study, idx=i)
	if i >= eval_start_idx and i < eval_start_idx + NUM_EVAL_FILES:
		eval_files.append(filename)
	else:
		training_files.append(filename)

rsrc = resource.RLIMIT_DATA
soft, hard = resource.getrlimit(rsrc)
print 'Soft limit starts as :', soft

resource.setrlimit(rsrc, (4*1024*1024*1024, hard))

soft, hard = resource.getrlimit(rsrc)
print 'Soft limit changed to :', soft

probe_indices = {} #maps probe_ids to indices
df = pd.read_csv("probe_id_indices")
print "reading probe indices"
for index, row in df.iterrows():
	probe_indices[row["probe_id"]] = row["index"] - 1

# numprobes = len(probe_indices)
numprobes = 25
print "done reading, number of probes: " + str(numprobes)

feature_columns = [
	tf.feature_column.numeric_column('betas', shape=(numprobes,), dtype=tf.float32)
]

label_column = tf.feature_column.numeric_column('label', shape=(1,), dtype=tf.int64)

features_spec = tf.feature_column.make_parse_example_spec(
	feature_columns + [label_column]
)

def input_fn(file_list):
	dataset = tf.contrib.data.make_batched_features_dataset(
		file_pattern=file_list,
		batch_size=16,
		features=features_spec,
		num_epochs=None,
		shuffle=False,
		)
	it = dataset.make_one_shot_iterator()
	features = it.get_next()
	labels = features.pop('label')
	#labels = tf.stack([label, 1-label])
	return features, labels

classifier = tf.estimator.DNNClassifier(
	feature_columns=feature_columns,
	hidden_units=[128],
	model_dir=modeldir,
	n_classes=2,
	optimizer=tf.train.ProximalAdagradOptimizer(
		learning_rate=learning_rate,
		l1_regularization_strength=regularization))

#training_files = sorted(glob.glob("BRCA/tcga_training_examples_*"))
#total_num_files = len(training_files)
#num_eval_files = total_num_files / 10
#eval_files = [training_files.pop() for i in range(num_eval_files)]

print "training files: " + ",".join(training_files)
print "evaluation files: " + ",".join(eval_files)

train_spec = tf.estimator.TrainSpec(
	input_fn=lambda: input_fn(training_files),
	max_steps=200000*1000
)

eval_spec = tf.estimator.EvalSpec(
	input_fn=lambda: input_fn(eval_files)
)

tf.estimator.train_and_evaluate(classifier, train_spec, eval_spec)

variable_names = classifier.get_variable_names()

for var_name in variable_names:
	if var_name.startswith("dnn/hiddenlayer_0/kernel"):
		values = classifier.get_variable_value(var_name)
		n_nonzero9 = 0
                n_nonzero6 = 0
                n_nonzero3 = 0
                n_nonzero0 = 0
                n_nonzero1 = 0
                n_nonzero2 = 0
		n_total_weights = 0
		n_nonzero_weights9 = 0
                n_nonzero_weights6 = 0
                n_nonzero_weights3 = 0
                n_nonzero_weights0 = 0
                n_nonzero_weights1 = 0
                n_nonzero_weights2 = 0

                wf_name = "weights_" + var_name.replace("/", "_") + ".csv"
                wf = open(wf_name, "w")

		for values_for_one_feature in values:
			has_nonzeros9 = False
                        has_nonzeros6 = False
                        has_nonzeros3 = False
                        has_nonzeros0 = False
                        has_nonzeros1 = False
                        has_nonzeros2 = False
			for x in values_for_one_feature:
                                wf.write("%f," % (x*1.0e6))
				n_total_weights += 1
                                if abs(x) > 1.0:
					has_nonzeros0 = True
					n_nonzero_weights0 += 1
                                elif abs(x) > 1.0e-1:
                                        has_nonzeros1 = True
                                        n_nonzero_weights1 += 1
                                elif abs(x) > 1.0e-2:
                                        has_nonzeros2 = True
                                        n_nonzero_weights2 += 1
                                elif abs(x) > 1.0e-3:
					has_nonzeros3 = True
					n_nonzero_weights3 += 1
                                elif abs(x) > 1.0e-6:
					has_nonzeros6 = True
                                        n_nonzero_weights6 += 1
                                elif abs(x) > 1.0e-9:
					has_nonzeros9 = True
					n_nonzero_weights9 += 1
                        wf.write("\n")
			if has_nonzeros9:
				n_nonzero9 += 1
                        if has_nonzeros6:
				n_nonzero6 += 1
			if has_nonzeros3:
				n_nonzero3 += 1
			if has_nonzeros0:
				n_nonzero0 += 1
			if has_nonzeros1:
				n_nonzero1 += 1
			if has_nonzeros2:
				n_nonzero2 += 1
		print var_name, "# inputs with nonzero weights", n_nonzero0, n_nonzero1, n_nonzero2, n_nonzero3, n_nonzero6, n_nonzero9
                print var_name, "# nonzero weights", n_total_weights, n_nonzero_weights0, n_nonzero_weights3, n_nonzero_weights6, n_nonzero_weights9
                wf.close()

########################################################################

weights_file = open("weights_" + modeldir + ".csv", "w")
var_name = "dnn/hiddenlayer_0/kernel"
values_for_all_features = classifier.get_variable_value(var_name)

weights_file.write(values_for_all_features)

weights_file.close()

########################################################################
