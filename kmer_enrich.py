#!/usr/bin/python

from __future__ import division
import argparse
import scipy.stats

__author__ = "Edward Kong"
__copyright__ = "Copyright 2016"
__credits__ = ["ELK"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Edward Kong"
__email__ = "edward.kong@yale.edu"

# Usage:      python kmer_enrich.py -i <input FASTA file> -k <kmer as a string>
# Example:    python kmer_enrich.py -i catsection.FASTA -k "CATTAG"
# Note:       Input: FASTA file, can be multiple sequences long.
#                    kmer: directly give a kmer as a string
#             Output: Prints the following information (for each sequence):
#                     Name: name of the fasta sequence from its identifier
#                     -- "expected" number of hits given independent bernoulli trials (strong assumption)
#                     -- actual number of hits (kmer matches to the base sequence)
#                     -- p value based on testing on the binomial distribution
#             catsection.FASTA: some example sequences


# Set up argument parser
parser = argparse.ArgumentParser(description='GCTtoNet')
parser.add_argument('-i', '--input', help='input file, FASTA nucleotide file', required=True)
parser.add_argument('-k', '--kmer', help='A string kmer', type=str, required=True)
args = parser.parse_args()


def run_kmer(inputfile, kmer):
    # list of sequence names and sequences contained in the FASTA file
    names = []
    sequences = []
    # index for the current sequence (set to -1 so the first loop starts at 0)
    curr_index = -1

    # Open the input file
    with open(inputfile, 'r') as f:
        # extract content
        rawdata = f.readlines()
        # parse lines in the rawdata
        for line in rawdata:
            # Read the > character as denoting a new nuceotide sequence
            if line[0] == '>':
                # store sequence name to later write to csv
                names.append(line.rstrip('\n'))
                curr_index += 1
                sequences.append('')
            else:
                # extend the current nucleotide sequence
                sequences[curr_index] += line.rstrip('\n')

    # Loop over all sequences in the list
    for s in range(0, len(sequences)):
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Calculate probability that a random sequence of length k matches our given kmer
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # (this is the parameter p in our binomial distribution model)
        # (A basic calculation implicitly assumes each nucleotide appears 25% of the time.
        # p = 0.25 ** k
        # But here, I take into account the overall base frequencies

        # First, establish base pair frequencies for the sequence:
        # bases := number of non N bases
        bases = len(sequences[s].rstrip('N'))
        A_frac = sequences[s].count('A') / bases
        T_frac = sequences[s].count('T') / bases
        C_frac = sequences[s].count('C') / bases
        G_frac = 1 - A_frac - T_frac - C_frac

        k = len(kmer)
        # initialize the prior probability as 1 (proba of finding a 0-kmer)
        p = 1
        # walk through the k-mer, using base frequencies to calculate the probability of
        # randomly encountering the particular k-mer in the reference sequence
        for bp in kmer:
            if bp == 'A':
                p *= A_frac
            elif bp == 'T':
                p *= T_frac
            elif bp == 'C':
                p *= C_frac
            elif bp == 'G':
                p *= G_frac
            else:
                p *= 0
                print('Error: k-mer contains a non-nucleotide!')
        print('Probability of a random match to the k-mer: ' + str(p))

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Parse the string to count actual / expected hits
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # n is the length of the sequence (so the total # of trials is n-k+1)
        # when calculating expected and actual hits, we omit any window with N's (poor quality /unknown bases)
        # Let w denote the number of 'valid' windows by the above criteria
        # so the (adjusted) expected number of hits for independent bernoulli trials is (w-k+1)*p
        n = len(sequences[s])

        # note: we could have just used sequences[s].count(kmer) to get the hits
        # however, we cannot just use len(sequences[s].rstrip('N')) to get w
        hits = 0
        w = 0
        # since we can't get the number of windows with a simple .count() operation, we use the following loop
        for i in range(0, n-k+1):
            # count the window if it doesn't contain any unknown bases
            if sequences[s][i:i+k].find('N') == -1:
                w += 1
                # conditional on the window not containing N's, see if the window matches
                if sequences[s][i:i+k] == kmer:
                    hits += 1

        # calculate expected hits:
        expected = (w-k+1)*p

        # Test for significance using a binomial test
        # arguments: actual hits, total trials, true probability, form of alt hypothesis
        p_value = scipy.stats.binom_test(hits, w-k+1, p, alternative='two-sided')

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Print output
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        print('Name: ' + names[s].replace('>', ''))
        print('Expected matches: ' + str(expected) +
              ' | Actual matches: ' + str(hits) +
              ' | P-value of two-sided test: ' + str(p_value))


# Run the k-mer enrichment algorithm
run_kmer(args.input, args.kmer)
