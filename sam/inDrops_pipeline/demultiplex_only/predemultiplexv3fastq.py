#import statements
import os, subprocess
import sys
import itertools
import operator
from collections import defaultdict, OrderedDict
import errno

# cPickle is a faster version of pickle that isn't installed in python3
# inserted try statement just in case
try:
   import cPickle as pickle
except:
   import pickle

from io import BytesIO

import numpy as np
import re
import shutil
import gzip

# product: product(A, B) returns the same as ((x,y) for x in A for y in B).
# combination: Return r length subsequences of elements from the input iterable.
from itertools import product, combinations
import time
import matplotlib
import yaml
import pysam

import tempfile
import string
from contextlib import contextmanager



# -----------------------
#
# Helper functions
#
# -----------------------



def string_hamming_distance(str1, str2):
    """
    Fast hamming distance over 2 strings known to be of same length.
    In information theory, the Hamming distance between two strings of equal 
    length is the number of positions at which the corresponding symbols 
    are different.

    eg "karolin" and "kathrin" is 3.
    """
    return sum(itertools.imap(operator.ne, str1, str2))

def rev_comp(seq):
    tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join(tbl[s] for s in seq[::-1])


def to_fastq(name, seq, qual):
    """
    Return string that can be written to fastQ file
    """
    return name+'\n'+seq+'\n+\n'+qual+'\n'
        
def seq_neighborhood(seq, n_subs=1):
    """
    Given a sequence, yield all sequences within n_subs substitutions of 
    that sequence by looping through each combination of base pairs within
    each combination of positions.
    """
    for positions in combinations(range(len(seq)), n_subs):
    # yields all unique combinations of indices for n_subs mutations
        for subs in product(*("ATGCN",)*n_subs):
        # yields all combinations of possible nucleotides for strings of length
        # n_subs
            seq_copy = list(seq)
            for p, s in zip(positions, subs):
                seq_copy[p] = s
            yield ''.join(seq_copy)
            
            
def weave_fastqs(fastqs):
    last_extension = [fn.split('.')[-1] for fn in fastqs]
    if all(ext == 'gz' for ext in last_extension):
        processes = [subprocess.Popen("gzip --stdout -d %s" % (fn), shell=True, stdout=subprocess.PIPE) for fn in fastqs]
        streams = [r.stdout for r in processes]
    elif all(ext == 'bz2' for ext in last_extension):
        processes = [subprocess.Popen("bzcat %s" % (fn), shell=True, stdout=subprocess.PIPE) for fn in fastqs]
        streams = [r.stdout for r in processes]
    elif all(ext == 'fastq' for ext in last_extension):
        streams = [open(fn, 'r') for fn in fastqs]
    else:
        raise("ERROR: Different files are compressed differently. Check input.")

    while True:
        names = [next(s)[:-1].split()[0] for s in streams]
        seqs = [next(s)[:-1] for s in streams]
        blanks = [next(s)[:-1]  for s in streams]
        quals = [next(s)[:-1]  for s in streams]
        assert all(name==names[0] for name in names)
        yield names, seqs, quals

    for s in streams:
        s.close()
            
            
def get_index_neighbors(libs):
    """
    input - list of DNA indexes
    output - dictionary with key = index neighbor, value = the only index the neighbor maps to
    """
    # Prepare error corrected index sets
    ix_neigh_2_ix = {}
    ix_neigh_2_ix = dict(zip(libs, libs))
    index_neighborhoods = [set(seq_neighborhood(lib, 1)) for lib in libs]
    for lib, clibs in zip(libs, index_neighborhoods):
        # Quick check that error-correction maps to a single index
        for clib in clibs:
            if sum(clib in hood for hood in index_neighborhoods)==1:
                ix_neigh_2_ix[clib] = lib
    return ix_neigh_2_ix


# -----------------------
#
# run
#
# -----------------------
start = time.time()

#load yaml
with open(sys.argv[1],'r') as stream:
    yamo = yaml.load(stream)

libraries = OrderedDict()
runs = OrderedDict()


for run in yamo['sequencing_runs']:
    if run['version']!='v3':
        continue
        
    savedir = os.path.join(run['dir'], 'demultiplexed')
    try:
        os.rmdir(savedir) #useful while developing the code, remove eventually
    except:
        pass
    os.mkdir(savedir) #for saving demultiplexed fastqs
    
    if 'split_affixes' in run:
        split_affixes = run['split_affixes']
    else:
        split_affixes = ['']
    run_libraries = [adict['library_index'] for adict in run['libraries']]
    run_library_names = [adict['library_name'] for adict in run['libraries']]
    ix_neigh_2_ix = get_index_neighbors(run_libraries)
    expected_indexes = ix_neigh_2_ix.keys()
    
    #prepare gzipper for saving. The resulting files will not be split by lane
    library_buffers = {}
    library_gzippers_R1 = {}
    library_gzippers_R2 = {}
    library_gzippers_R3 = {}
    library_gzippers_R4 = {}
    for lib,name in zip(run_libraries,run_library_names):
        library_gzippers_R1[lib] = {}
        library_gzippers_R2[lib] = {}
        library_gzippers_R3[lib] = {}
        library_gzippers_R4[lib] = {}
        libdir = savedir+'/%s_%s'%(lib,name)
        os.mkdir(libdir)
        for affix in split_affixes:
            library_gzippers_R1[lib][affix] = gzip.open(libdir+'/%s_%s_R1.fastq.gz'%(lib,affix),mode='a')
            library_gzippers_R2[lib][affix] = gzip.open(libdir+'/%s_%s_R2.fastq.gz'%(lib,affix),mode='a')
            library_gzippers_R3[lib][affix] = gzip.open(libdir+'/%s_%s_R3.fastq.gz'%(lib,affix),mode='a')
            library_gzippers_R4[lib][affix] = gzip.open(libdir+'/%s_%s_R4.fastq.gz'%(lib,affix),mode='a')
        
    #for saving read with ambiguous library indexes
    ambdir = savedir+'/ambiguous_library_index'
    os.mkdir(ambdir)
    ambig_library_gzippers_R1 = {}
    ambig_library_gzippers_R2 = {}
    ambig_library_gzippers_R3 = {}
    ambig_library_gzippers_R4 = {}
    for affix in split_affixes:
        ambig_library_gzippers_R1[affix] = gzip.open(ambdir+'/ambig_lib_index_%s_R1.fastq.gz'%affix,mode='a')
        ambig_library_gzippers_R2[affix] = gzip.open(ambdir+'/ambig_lib_index_%s_R2.fastq.gz'%affix,mode='a')
        ambig_library_gzippers_R3[affix] = gzip.open(ambdir+'/ambig_lib_index_%s_R3.fastq.gz'%affix,mode='a')
        ambig_library_gzippers_R4[affix] = gzip.open(ambdir+'/ambig_lib_index_%s_R4.fastq.gz'%affix,mode='a')
    
    for affix in split_affixes:
        input_filename = os.path.join(run['dir'], run['fastq_path'].format(split_affix=affix, read='{read}'))
        
        # Paths for the 4 expected FastQs
        input_fastqs = []
        for r in ['R1', 'R2', 'R3', 'R4']:
            input_fastqs.append(input_filename.format(read=r))

        for names, seqs, quals in weave_fastqs(input_fastqs):
            
            # Python 3 compatibility in mind!
            seqs = [s.decode('utf-8') for s in seqs]

            if seqs[2] in expected_indexes: #funny that R3 contains the index sequence, this read is called R4 on Basespace.
                lib = ix_neigh_2_ix[seqs[2]]
                library_gzippers_R1[lib][affix].write(to_fastq(names[0], seqs[0], quals[0]))
                library_gzippers_R2[lib][affix].write(to_fastq(names[1], seqs[1], quals[1]))
                library_gzippers_R3[lib][affix].write(to_fastq(names[2], seqs[2], quals[2]))
                library_gzippers_R4[lib][affix].write(to_fastq(names[3], seqs[3], quals[3]))
            else:
                ambig_library_gzippers_R1[affix].write(to_fastq(names[0], seqs[0], quals[0]))
                ambig_library_gzippers_R2[affix].write(to_fastq(names[1], seqs[1], quals[1]))
                ambig_library_gzippers_R3[affix].write(to_fastq(names[2], seqs[2], quals[2]))
                ambig_library_gzippers_R4[affix].write(to_fastq(names[3], seqs[3], quals[3]))
                
    for lib in run_libraries:
        for affix in split_affixes:
            library_gzippers_R1[lib][affix].close()
            library_gzippers_R2[lib][affix].close()
            library_gzippers_R3[lib][affix].close()
            library_gzippers_R4[lib][affix].close()
        
    for affix in split_affixes:
        ambig_library_gzippers_R1[affix].close()
        ambig_library_gzippers_R2[affix].close()
        ambig_library_gzippers_R3[affix].close()
        ambig_library_gzippers_R4[affix].close()

print time.time()-start
