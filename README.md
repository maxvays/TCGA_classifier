# TCGA_classifier
Identifying epigenetic signature of breast cancer using TCGA dataset and TensorFlow DNNClassifier

## shard_betas.py
Reads hg38_training_examples_* files (the sharded hg38_training_examples table
downloaded from the TCGA dataset in BigQuery) with the following columns:
disease_code, aliquot_barcode, probe_id, beta_value.
Also reads the list of aliquot barcodes for each disease code, and assigns
each aliquot to one of the sharded methylation betas files, so as to have 10
aliquots per shard.
Generates sharded methylation betas files, writing betas for each aliquot into
the assigned shard file.

## horizontal_betas.py
Generates files with methylation betas for each aliquot in one row,
as comma-separated values.

## training_examples.py
Generates training examples files, with labels set to 1 for training examples
representing methylation betas for the target study (BRCA) and labels set to 0
for training examples representing methylation betas for all other studies.
Methylation betas for training examples are retrieved from horizontal-betas
files.

## TFRecord_file_writer.py
Reads training examples files (methylation betas and labels as comma-separated
values, one row per aliquot) and generates TensorFlow TFRecord files, containing
training examples (methylation betas and corresponding 0/1 label values) in the
format that can be fed directly into a TensorFlow classifier.

## TFRecord_file_writer_filtered.py
Reads training examples files (methylation betas and labels as comma-separated
values, one row per aliquot) and generates TensorFlow TFRecord files, containing
training examples (methylation betas and corresponding 0/1 label values) using
only methylation betas for the 25 probe IDs with the highest weights in the
neural network of the trained classifier, as identified in the previous step.

## classifier.py
Trains TensorFlow DNNClassifier on training examples in TFRecord format, with
each training example comprising methylation betas for one aliquot and the 0/1
label.
Outputs statistics on the number of input features with non-zero weights in the
hidden layer (dnn/hiddenlayer_0/kernel).
Upon completion of training, writes the weights for all input features in the
hidden layer (dnn/hiddenlayer_0/kernel) into a separate file for analysis.

## analyze_weights.py
Sorts probe IDs by the sums of the absolute values of weights in the
classifier's hidden layer (dnn/hiddenlayer_0/kernel), and outputs
(probe_index, probe_id, sum_weights) for the top 1000 probe IDs.

## analyze_betas_for_filtered_probes.py
Calculate average, standard deviation, t-stat, and p-value for filtered
(top 100) probe IDs in BRCA vs. non-BRCA training examples.

