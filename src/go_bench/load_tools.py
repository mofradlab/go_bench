import os
from collections import defaultdict
import gzip

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, dok_matrix, lil_matrix
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#from Bio.UniProt.GOA import GAF20FIELDS
GAF20FIELDS = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 
            'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence', 
            'With', 'Aspect', 'DB_Object_Name', 'Synonym', 
            'DB_Object_Type', 'Taxon_ID', 'Date', 'Assigned_By', 
            'Annotation_Extension', 'Gene_Product_Form_ID']

def to_data_num(year, month, day):
    return int(year)*10000+int(month)*100+int(day)

#The annotations dict maps protein ids to a list of (GO annotation, evidence_code) tuples. 

neg_qualifiers = {'NOT|colocalizes_with', 'NOT|acts_upstream_of_or_within', 
                  'NOT|acts_upstream_of', 'NOT|acts_upstream_of_or_within_negative_effect', 
                  'NOT|enables', 'NOT|part_of', 'NOT|is_active_in', 'NOT|located_in', 
                  'NOT|contributes_to', 'NOT|involved_in', 'acts_upstream_of_or_within_negative_effect', 
                  'acts_upstream_of_negative_effect'}

def load_protein_annotations(goa_path, annotation_codes, min_date=None, max_date=None):
    """
    Loading function which takes in a GOA tab formatted file in which each row contains a protein annotation and
    outputs a dictionary mapping each protein ID in the file to each annotation for that protein included in the 
    file (after filtering). 
  
    Parameters:
    goa_path (string): File path leading to a GOA formatted tabular file. 
    Example at ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/old//UNIPROT/goa_uniprot_all.gaf.203.gz

    annotations_codes (Set[String]): List of annotations codes to count as valid. If the annotation code 
    for a row in the GOA file isn't included in annotation_codes, the row is ignored. Valid annotation codes,
    and groups of similar annotation codes, are included as variables in gobench.utils, with the full list being
    ('ISO', 'IGI', 'ISA', 'IGC', 'IEP', 'RCA', 'HDA', 'HGI', 'IKR', 'TAS', 'HEP', 'ND', 'IBA', 'IMP', 'EXP', 
    'IDA', 'IC', 'ISM', 'ISS', 'NAS', 'IRD', 'IEA', 'IPI', 'HMP')

    min_date/max_date (Int|None): Minimum and maximum dates for GOA rows included in result. Must be formatted according to 
    load_tools.to_data_num(year, month, day). 
  
    Returns:
    Dict[String, List[String]]: Dictionary mapping protein ids to a list of gene ontology ids. 
    """
    annotation_codes = set(annotation_codes)
    annot_dict = defaultdict(set)
    df_iter = pd.read_csv(goa_path, dtype=str,
                          sep='\t',
                          comment='!',
                          names=GAF20FIELDS,
                          chunksize=int(1e6)) 
    
    if(min_date is None and max_date is None):
        for zdf in df_iter:
            # For now, remove all with a qualifier
            zdf = zdf[zdf.Evidence.isin(annotation_codes) & ~zdf.Qualifier.isin(neg_qualifiers)]
            for tup in zdf.itertuples():
                annot_dict[tup[2]].add(tup[5])
    else:
        if(min_date is None):
            min_date = 0
        if(max_date is None):
            max_date = 1e10
            
        for zdf in df_iter:
            # For now, remove all with a qualifier
            dates = zdf.Date.astype(int)
            zdf = zdf[zdf.Evidence.isin(annotation_codes) & ~zdf.Qualifier.isin(neg_qualifiers) & dates.ge(min_date) & dates.le(max_date)]
            for tup in zdf.itertuples():
                annot_dict[tup[2]].add(tup[5])
    return annot_dict


def load_protein_sequences(path, prot_whitelist=None):
    """
    Loads proteins sequences from a tab or fasta file and filters by a whitelist set. Returns sequences and 
    prot_ids representing the protein id associated with each sequence, both in the same order. 

    Parameters:
    path (string): Path to .fasta or .tab file, either gzipped with an optional .gz extension
    prot_whitelist (Iterable[String]): An optional set of protein ids to load from the file, for cases where
    only a known subset of proteins are needed. 
    """
    if(prot_whitelist):
        prot_whitelist = set(prot_whitelist)

    sequences = []
    prot_ids = []
    if(path[-3:] == "tab" or path[-6:] == "tab.gz"):
        
        df_iter = pd.read_csv(path, dtype=str,
                                sep='\t',
                                header=0,
                                chunksize=int(1e5)) 
        print("Loading sequences")
        count = 0
        for df in df_iter:
            print(f"{count*1e5} rows read")
            if(prot_whitelist):
                df = df[df["Entry"].isin(prot_whitelist)]
            sequences.extend([s.lower() for s in df["Sequence"]])
            prot_ids.extend(df["Entry"])
            count += 1
    elif(path[-8:] == "fasta.gz"):
        with gzip.open(path, 'rt') as fp:
            seqs = SeqIO.parse(fp, 'fasta')
            count = 0
            for rec in seqs:
                seq_id = rec.id
                if('|' in seq_id):
                    seq_id = rec.id.split('|')[1]
                seq = rec.seq
                if(seq_id in prot_whitelist):
                    sequences.append(str(seq.lower()))
                    prot_ids.append(seq_id)
                if(count % 100000 == 0):
                    print(f"{count} proteins processed")
                count += 1
    elif(path[-5:] == "fasta"):
        with open(path, 'r') as fp:
            seqs = SeqIO.parse(fp, 'fasta')
            count = 0
            for rec in seqs:
                seq_id = rec.id
                if('|' in seq_id):
                    seq_id = rec.id.split('|')[1]
                seq = rec.seq
                if(prot_whitelist is None or seq_id in prot_whitelist):
                    sequences.append(str(seq.lower()))
                    prot_ids.append(seq_id)
                if(count % 10000 == 0):
                    print(f"{count} proteins processed")
                count += 1
    return sequences, prot_ids

#Data management methods and classes
def load_GO_tsv_file(path):
    """
    Loading function for the GOBench specific GO tsv file. Each row from the file must give information
    for a unique protein, in the format of a protein id followed by a list of comma separated gene ontology ids that
    the protein has been annotated with. E.g. Protein_ID  GO:00002316,GO:00005547,GO:00008904
    The output is a dictionary mapping protein ids to lists of associated gene ontology ids, in the same format as that
    given by load_tools.load_protein_annotations. 
    
    """
    prot_dict = {}
    df_iter = pd.read_csv(path, dtype=str,
                          sep='\t',
                          header=0,
                          chunksize=int(1e4))
    for zdf in df_iter:
        for tup in zdf.itertuples():
            prot_id = tup[1]
            annotations = tup[2]
            prot_dict[prot_id] = list(annotations.split(", ")) if isinstance(annotations, str) else []
    return prot_dict

#Convert a protein annotation dict to a sparse matrix. 
def convert_to_sparse_matrix(protein_annotation_dict, term_list, prot_id_list):
    """
    A helper function for converting a dictionary mapping proteins to gene ontology annotations
    to a sparse binary matrix, with a row for each protein and a column for each gene ontology term. 
    A true value in the matrix at row i and col j indicates that protein i is annotated with GO term j.

    Rows and columns are ordered by their respective orders in prot_id_list and term_list. If a protein or term
    not included in their list, it will be excluded from the resulting matrix. 

    Parameters:
    protein_annotation_dict (Dict[String, List[String]]): Protein annotation dictionary
    term_list (List[String]): An ordered list of gene ontology terms.
    prot_id_list (List[String]): An ordered list of protein ids. 

    """
    term_col_mappings = {term:i for i, term in enumerate(term_list)}
    prot_row_mappings = {prot:i for i, prot in enumerate(prot_id_list)}

    labels = lil_matrix((len(prot_id_list), len(term_list)), dtype=np.int8)

    for place, prot_id in enumerate(prot_id_list):
        if(prot_id in protein_annotation_dict):
            for go_id in protein_annotation_dict[prot_id]:
                if(go_id in term_col_mappings):
                    labels[place, term_col_mappings[go_id]] = 1
    labels = labels.tocsr()
    return labels

def read_sparse(fn, prot_rows, GO_cols):
    """
    Reads a CAFA formatted file making GO annotation predictions for a set of proteins. 
    Each row of the file should have the format ProtID GOID Score. 

    Output is a sparse floating point matrix. A non-zero value s in the matrix at row i and col j 
    indicates that protein i is annotated with GO term j with score s.

    Rows and columns are ordered by their respective orders in prot_rows and GO_cols. If a protein or term
    not included in their list, it will be excluded from the resulting matrix. 
    """
    prm = {prot:i for i, prot in enumerate(prot_rows)}
    tcm = {term:i for i, term in enumerate(GO_cols)}
    sparse_probs = dok_matrix((len(prot_rows), len(GO_cols)))
    df = pd.read_csv(fn)
    for (i, prot, go_id, prob) in df.itertuples():
        if(prot in prm and go_id in tcm):
            sparse_probs[prm[prot], tcm[go_id]] = prob
    return csr_matrix(sparse_probs)

def write_sparse(fn, preds, prot_rows, GO_cols, go, min_certainty):
    """
    Convert a sparse matrix of GO annotations to a CAFA formatted tsv file.

    Parameters:
    fn (String): Path to output file. 
    preds (csr_matrix): Sparse matrix with non-zero entries at row i, col j indicating an annotation for protein i with
    GO term j. 
    prot_rows (List[String]): An ordered list of protein ids. 
    GO_cols (List[String]): An ordered list of GO terms. 
    go (GODag): Gene Ontology object from the GODag library. 
    E.g. from goatools.obo_parser import GODag
    godag = GODag('data/go.obo')
    min_certainty (Float): Cutoff at which low scores should be rounded down to zero to increase sparsity. 
    """
    with open(fn, mode='w') as f:
        f.write("g,t,s\n")
        for row, col in zip(*preds.nonzero()):
            prot_id = prot_rows[row]
            go_id = GO_cols[col]
            val = preds[row, col]
            if(val > min_certainty and go_id in go):
                f.write(f"{prot_id},{go_id},{val}\n")
