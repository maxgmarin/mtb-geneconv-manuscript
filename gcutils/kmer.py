# pgqc/module1.py


import pandas as pd
import numpy as np

#import time
from tqdm import tqdm

import screed
import mmh3


##### k-mer functions #####

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1
    
    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)
        
    return kmers


def read_kmers_from_file(filename, ksize):
    all_kmers = []
    for record in screed.open(filename):
        sequence = record.sequence
        
        kmers = build_kmers(sequence, ksize)
        all_kmers += kmers

    return all_kmers


def hash_kmer(kmer):
    # calculate the reverse complement
    rc_kmer = screed.rc(kmer)
    
    # determine whether original k-mer or reverse complement is lesser
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer
        
    # calculate murmurhash using a hash seed of 42
    hash = mmh3.hash64(canonical_kmer, 42)[0]
    if hash < 0: hash += 2**64
        
    # done
    return hash

## Note that hashing collections of k-mers doesn't change Jaccard calculations:
def hash_kmers(kmers):
    hashes = []
    for kmer in kmers:
        hashes.append(hash_kmer(kmer))
    return hashes

def hash_kmers_ToSet(kmers):
    hashes = set()
    for kmer in kmers:
        hashes.add(hash_kmer(kmer))
    return hashes

def hash_kmers_ToUnqNP(kmers):
    hashes = []
    for kmer in kmers:
        hashes.append(hash_kmer(kmer))
        
    return np.unique(np.array(hashes))

    

def read_kmers_from_file_ToHashesDict(filename, ksize):

    all_hashes_Set_Dict = {}
    seqLen_Dict = {}
    
    NumParsedRecords = 0
    
    for record in screed.open(filename):
        
        ShortName = record.name.split(" ")[-1]

        NumParsedRecords += 1
        sequence = record.sequence

        kmers = build_kmers(sequence, ksize)
        hashes_Set = hash_kmers_ToSet(kmers)
        
        all_hashes_Set_Dict[ShortName] = hashes_Set
        seqLen_Dict[ShortName] = len(sequence)

    print(NumParsedRecords, " total records were parsed")
    
    return all_hashes_Set_Dict, seqLen_Dict



def jaccard_containment_FromSets(a, b):
    '''
    This function returns the Jaccard Containment between sets a and b.
    '''
    
    intersection = len(a.intersection(b))
    
    return intersection / len(a)

def jaccard_similarity_FromSets(a, b):
    '''
    This function returns the Jaccard Similarity between sets a and b.
    '''
    intersection = len(a.intersection(b))
    union = len(a.union(b))
    
    return intersection / union

def jaccard_containment_MaxVal_FromSets(a, b):
    '''
    This function returns the maximum possible Jaccard Containment between sets a and b.
    '''
    
    intersection = len(a.intersection(b))

    min_Len = min(len(a), len(b) )

    return intersection / min_Len

########################################################################################################################




































