# 2.6_kmer_enrichment
Calculate the enrichment of a particular k-mer given a reference sequence. Also gives the p-value based on the native base-pair frequencies of the query sequence and a binomial distribution 

Usage:      python kmer_enrich.py -i <input FASTA file> -k <kmer as a string>

Example:    python kmer_enrich.py -i catsection.FASTA -k "CATTAG"

Note:       Input: FASTA file, can be multiple sequences, as long as each sequence begins with '>' as is standard. 

  --kmer: directly give a kmer as a string

Output: Prints the following information (for each sequence):

  Name: name of the fasta sequence from its identifier

  -- "expected" number of hits given independent bernoulli trials 

  -- actual number of hits (kmer matches to the base sequence)

  -- p value based on testing on the binomial distribution

  catsection.FASTA: some example sequences

  catgenome.FASTA: cat genome as a reference
