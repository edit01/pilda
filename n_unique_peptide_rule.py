#!/usr/bin/python

import argparse
import re
from collections import defaultdict

mol_weights = {
    'I':131.1736, 'L':131.1736, 'K':146.1882, 'M':149.2124,
    'F':165.1900, 'T':119.1197, 'W':204.2262, 'V':117.1469,
    'R':174.2017, 'H':155.1552, 'A':89.0935,  'N':132.1184,
    'D':133.1032, 'C':121.1590, 'E':147.1299, 'Q':146.1451,
    'G':75.0669,  'P':115.1310, 'S':105.0930, 'Y':181.1894}

#Given a string of amino acid, compute the molecular weight
def compute_weight(aa_string):
    weight = 0
    for aa in aa_string:
        weight += mol_weights[aa]
    return weight

def parse_sample_file(in_file, threshold):
    with open(getattr(args, 'in'), 'r') as sf:
        col_name = sf.readline().strip().split('\t')
        col_index = {col:i for i, col in enumerate(col_name)}                      

        peptides = defaultdict(list)#peptide => [run_id, run_id, run_id]
        for line in sf:
            peptide = line.strip().split('\t')
            sequence = peptide[col_index['peptide']]
            q_value = peptide[col_index['qvalue']]
            run_id = peptide[col_index['run_identifier']]
            mod_string = peptide[col_index['modification_string']]
            #print q_value
            if float(q_value) <= float(threshold):
                #print sequence, q_value, mod_string
                peptides[sequence].append(run_id)                
        
    return peptides

def parse_protein_library(protein_database, sample_file):
    with open(protein_database, 'r') as pd:

        peptides = defaultdict(list)#peptide => [protein1, protein2, protein3] AKA multimapping        
        for line in pd:
            line = line.strip()

            #Line 1.Protein identifier
            #Line 2.Peptide sequence
            if re.match('^>', line):
                protein_id = line.lstrip('>')
            else:
                pep_seq = line
                peptides[pep_seq].append(protein_id)
        
        #Filter out multimapping peptides from within the library
        filtered_dict = {sequence:proteins[0] for sequence, proteins in peptides.items() if len(proteins) == 1}        
        
        #Filter out peptides we don't care about. Keep peptides that we see in sample file
        #This also leaves only uniquely mapping peptides in the sample file
        filtered_dict2 = { sequence:protein for sequence, protein in filtered_dict.items() if sequence in sample_file}
        

    return filtered_dict2

def process(args):

    #Parse/Filter by qval, sample file
    sample_file = parse_sample_file(getattr(args, 'in'), args.threshold)

    #Parse/Filter/Filter trypsin protein library    
    protein_lib = parse_protein_library(args.protein_database, sample_file)

    '''There is a mismatch between sample_file & protein_lib. Meaning, there are
    peptide sequences in the sample file that are not found in the protein library.'''
    #Peptides in the sample file are thrown out if they are not found within our library
    sample_file_2 = {sequence:runs for sequence, runs in sample_file.items() if sequence in protein_lib}

    #Sort peptides into runs + Do Spectral counting
    by_runs = defaultdict(list)
    for sequence, runs in sample_file_2.items():
        spectral_counts = defaultdict(int) #of times this sequence appeared in this run AKA spectral count
        for run in runs:
            spectral_counts[run] += 1
        
        for run, count in spectral_counts.items():
            by_runs[run].append((sequence, count))

    
    #Predict the proteins
    with open(args.out, 'w') as fh:
        fh.write("%s\t%s\t%s\n" % ("run_identifier", "protein", "total_spectral_count"))

        for run, peptides in by_runs.items():

            #Predict the proteins for this run, from this set of peptides
            predicted = defaultdict(list)
            for (sequence, count) in peptides:
                #Protein lookup
                protein_id = protein_lib[sequence]
                predicted[protein_id].append((sequence, count))
            
            #Write out to file
            for protein, peptides in predicted.items():
                if len(peptides) >= int(args.min):
                    #Sum the spectral counts from all peptides that contribute to this protein
                    total_spectral_count = 0
                    for (sequence, count) in peptides:                    
                        total_spectral_count += count

                    fh.write("%s\t%s\t%s\n" % (run, protein, total_spectral_count))

    return None

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-in')
    parser.add_argument('-out')
    parser.add_argument('-min')
    parser.add_argument('-threshold')
    parser.add_argument('-protein_database')
    parser.add_argument('-all_params')

    args = parser.parse_args()

    KMIN_UNIQUE_PARAMS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    THRESHOLD_PARAMS = [.01, .05, .10, .15, .20, .25, .30, .35, .40, .45, .50]
    
    if args.all_params and int(args.all_params) == 1:
        for k in KMIN_UNIQUE_PARAMS:
            args.min = k            
            for t in THRESHOLD_PARAMS:
                args.threshold = t
                args.out = "min.%s.threshold.%s.ISB18.txt" % (args.min, int(args.threshold*100))
                print "Processing..:", args.out
                process(args)
    else:
        process(args)











