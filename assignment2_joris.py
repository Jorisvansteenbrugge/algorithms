#!/usr/bin/env python
"""
Author: Joris van Steenbrugge
Student number: 950416798110
Implementation of the SSAHA algorithm

Hints:
- write a function to generate a hash table of k-mers from a set of sequences
- write a function to return a sorted list of hits in the hash table for a query sequence
- write a function to find the best hit
- write a fasta parser to read in the Arabidopsis data

"""
 
""" 
Definitions:
    Q               --  query sequence
    D               --  Data base of sequences D = {S1, s2, .., Sd}
                            each sequence in D is labeled with an index i
    ktuple          -- contiguous sequences of k bases long. 
                        each sequence of length n has (n-k+1) overlapping 
                        tuples
    offset wj(S)    --  the k-tuple of S that has offset j
""" 

# import statements
from sys import argv
from Bio import SeqIO

def parse_fasta(file_name):
    # sequences = []
    # seq = ""
    # with open(file_name) as in_file:
    #     for line in in_file:
    #         if not line.strip():
    #             continue
    #         if line.startswith('>'):
    #             if seq != "":
    #                 sequences.append(seq)
    #                 seq = ""
    #         else:
    #             seq += line.strip()

    # sequences.append(seq)
    # return sequences
    return [entry.seq for entry in SeqIO.parse(file_name,'fasta')]
    #for entry in SeqIO.parse(file_name,'fasta'):
     #   sequences.append(entry.seq)

def create_hashtable(seqs, k = 2):
    """Returns a hashtable of non-overlaping k-mers for a database of sequences.

        The hashtable is implemented in a python dictionary where k-mers are 
        keys and offsets are values. The offsets are stored in a tuple where 
        the first index corresponds to the sequence in the database and the 
        second index to the starting index of the k-mer in that sequence.

        Keyword Arguments:
            seqs -- list, database of DNA sequences
            k -- integer, wordsize 

        Returns:
            hashtable -- dictionary, {kmer:offset}

    """
    hashtable = {}
    for i, seq in enumerate(seqs): 
        for x in range(0, len(seq), k):
            kmer = seq[x: x + k]
            
            """
            "It's easier to ask forgiveness than it is to get permission."
             ~ Grace Hopper
            """
            if len(kmer) != k:
                continue
            try:
                hashtable[kmer].append((i, x))
            except KeyError:
                hashtable[kmer] = [(i, x)]
        
    return hashtable

def sequence_search(query, hashtable, k = 2):
    M = []
    for t in range(0, len(query) - k + 1):
        kmer = query[t:t+k]

        #If the key is not present skip this iteration step
        if not kmer in hashtable:
            continue

        r_positions = hashtable[kmer]

        for hit in r_positions:
            index = hit[0]
            offset = hit[1]
            shift = offset - t

            M.append((index, shift, offset))

    return sorted(M, key = lambda x: (x[0], x[1]))

def get_longest_match(M_list, l_seqs):
    all_matches = []
    for i in range(l_seqs):
        cur_index_matches = [x for x in M_list if x[0] == i]
        uniq_offsets = set([x[1] for x in cur_index_matches])

        for offset in uniq_offsets:
            matches = [match for match in cur_index_matches if \
                match[1]  == offset]
            all_matches.append(matches)

    return max(all_matches, key=len)

def pretty_print_table(hashtable):
    for key,val in hashtable.items():
        values = [",".join(map(str, v)) for v in val]
        print key + "\t" + "\t".join(values)

def print_alignment(longest_match, seqs, k = 2):
    subject = seqs[longest_match[0][0]]

    matching_query = [subject[match[2]:match[2]+k] for \
        match in longest_match]

    matching_query = "".join(matching_query)
    pos = subject.index(matching_query)


    print subject
    print " " * pos + "|" *len(matching_query)
    print " " * pos + matching_query

def SSAHA(query, seqs, k = 2, print_info = True):
    hashtable = create_hashtable(seqs, k)
    M_list = sequence_search(query, hashtable, k)
    longest_match = get_longest_match(M_list, len(seqs))


    if print_info:
        ##print Hashtable
        print("HashTable: ")
        pretty_print_table(hashtable)

        
        print("Number of hits: {}".format(str(len(M_list))))
        print("First hit: {}".format(M_list[0]))
        print("Last hit: {}".format(M_list[-1]))
    
    #Print the actual alignment
    print_alignment(longest_match, seqs, k)

if __name__ == "__main__":
    # the code below should produce the results necessary to answer the questions
    # in other words, if we run your code, we should see the data that you used
    # to answer the questions
    
    ara_genome_file = '/home/joris/algorithms/TAIR10.fasta'
    query = 'TGCAACAT'
    
    s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
    seqs = [s1,s2,s3]

    #With test database
    SSAHA(query, seqs, 2, print_info = True)
    
    genome = parse_fasta(ara_genome_file)
    #Tair database
    SSAHA(query, genome, 15, print_info = False)

    
