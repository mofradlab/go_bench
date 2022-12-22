# GO Bench
go_bench is a python package for generating Gene Ontology datasets from the SwissProt and GOA databases. It contains functions for downloading and processing protein annotations to generate training, validation, and testing sets from predetermined random splits. 

General usage, with download links for required datasets, is shown in gen_datasets.ipynb. To replicate results of [gobench.org](https://www.gobench.org/), it is important to use the same train/validation/test split used the internal server, instead of the script included in this package. These are included in the data_split directory. 

# Install
go_bench can be installed from git with `pip install git+https://github.com/mofradlab/go_bench.git`. It requires python >= 3.6, and lists numpy, pandas, scipy, and goatools as dependencies. 

# Data Requirements
Most go_bench functions parse or manipulate gene ontology data. Raw data can be downloaded from SwissProt, the GOA database, and the Gene Ontology, with links for all requirements included in the go_bench_usage.ipynb notebook. Reproducing the datasets from [gobench.org](https://www.gobench.org/) requires an additional set of train/validation/test splits included in the data_splits directory. 

# Usage
The full process of generating a go_bench dataset from scratch is included in go_bench_usage.ipynb. Examples of parsing go_bench datasets and processing them for machine learning are included in https://github.com/mofradlab/GO_Bench_Sample, or alternatively, our colab notebook: https://colab.research.google.com/drive/1af4WHuNXsn4O6L9H1UvVNq8-2qdBCIGP?usp=sharing. 

We divide the gobench codebase into several sections, focused on loading and manipulating gene ontology annotation data, evaluating model predictions, and outputting datasets. Important functions, with context, are outlined for each section below. **For details, see function doc-strings for extensive notes on usage.**

## Data Loading
### load_protein_annotations
Loading function which takes in a GOA tab formatted file in which each row contains a protein annotation and
outputs a dictionary mapping each protein ID in the file to each annotation for that protein included in the 
file (after filtering). 
### load_protein_sequences
Loads proteins sequences from a tab or fasta file and filters by a whitelist set. Returns sequences and 
prot_ids representing the protein id associated with each sequence, both in the same order. 
### load_GO_tsv_file
Loading function for the GOBench specific GO tsv file. Each row from the file must give information
for a unique protein, in the format of a protein id followed by a list of comma separated gene ontology ids that
the protein has been annotated with. E.g. Protein_ID  GO:00002316,GO:00005547,GO:00008904
The output is a dictionary mapping protein ids to lists of associated gene ontology ids, in the same format as that
given by load_tools.load_protein_annotations. 
### convert_to_sparse_matrix
A helper function for converting a dictionary mapping proteins to gene ontology annotations
to a sparse binary matrix, with a row for each protein and a column for each gene ontology term. 
A true value in the matrix at row i and col j indicates that protein i is annotated with GO term j.
### read_sparse
Reads a CAFA formatted file making GO annotation predictions for a set of proteins. 
Each row of the file should have the format ProtID GOID Score. Output is a sparse 
floating point matrix. A non-zero value s in the matrix at row i and col j indicates 
that protein i is annotated with GO term j with score s.
### write_sparse
Convert a sparse matrix of GO annotations to a CAFA formatted tsv file.

## Gene Annotation Prediction Metrics
### threshold_stats
Takes in model predictions and true classifications for a multiclass, multilabel GO annotation 
problem. Calculates F-score, S-score, and F-est metrics for a range of 100 thresholds between 0 and 1. 
Returns the range of values for each so that F-max and S-min can be calculated, and results can be plotted. 

## GO Bench Pipeline Processing
### pipeline
Comprehensive function to parse a GOA file, filter the result by evidence codes, term frequency, namespace, and date,
propogate terms through the GO tree, and output GO Bench tsv files. Used in gobench.org to generate outputs. 

# Contact
Please send questions or suggestions to amdickson (at) berkeley.edu. Also see our lab website at https://biomechanics.berkeley.edu/. 

