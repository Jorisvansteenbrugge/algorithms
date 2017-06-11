#!/usr/bin/env python
"""
Author: Joris van Steenbrugge
Student number: 950416798110
Implementation of the SSAHA algorithm
"""

from sys import argv, exit
from Bio import SeqIO

def parse_fasta(file_name):
    """Returns a list with sequences from a fasta file

        The fasta file is parsed using Biopython as fasta parsing is
        trivial for this assignment.

        Keyword Arguments:
            file_name -- string, file path of the fasta file
    """
    return [str(entry.seq) for entry in SeqIO.parse(file_name, 'fasta')]

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
        for x in range(0, len(seq), k): #create all chained k-mers for seq
            kmer = seq[x: x + k]
             
            # Sometimes the last kmer is truncated and not the actual length
            if len(kmer) != k:
                continue

            """
            "It's easier to ask forgiveness than it is to get permission."
             ~ Grace Hopper
             (also faster in python)
            """
            try:
                hashtable[kmer].append((i, x))
            except KeyError: #happens if the key does not yet exist
                hashtable[kmer] = [(i, x)]
        
    return hashtable

def sequence_search(query, hashtable, k = 2):
    """Return a list containing tuples of hits

        The tuples contain (index, shift, offset) each.

        Keyword Arguments:
            query -- string, DNA sequence to compare with the database
            hashtable -- dictionary, {kmer: offset(s)}
            k -- int, size the kmers should be.
    """
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
    """Returns the longest match of the query and the database.

        Keyword Arguments:
            M_list -- list, matches between the query and the database
            l_seqs -- int, The number of sequences in the database

    """
    all_matches = []
    for i in range(l_seqs):
        cur_index_matches = [x for x in M_list if x[0] == i]
        uniq_offsets = set([x[1] for x in cur_index_matches])

        for offset in uniq_offsets:
            matches = [match for match in cur_index_matches if \
                match[1]  == offset]
            all_matches.append(matches)
    max_val = None
    try:
        max_val = max(all_matches, key=len)
    except ValueError:
        print("No match was found")
    return max_val

def pretty_print_table(hashtable):
    """Converts the hashtable to a printable format.
    
    The output will be comparable to table 1 in N. Zemin et al.

        Keyword Arguments:
            hashtable -- dict, {kmer:offset}
    """
    for key,val in hashtable.items():
        values = [",".join(map(str, v)) for v in val]
        print(key + "\t" + "\t".join(values))

def print_alignment(query, longest_match, seqs, k = 2):
    """Prints the alignment in a human interpretable way.

        Not suitable for large sequences as the console will be flooded 
        with output.

        Keyword Arguments:
            query --  string, DNA sequence
            longest_match -- list, positions of the longest match
            seqs -- list, DNA sequences in the database

    """
    subject = seqs[longest_match[0][0]]

    matching_query = [subject[match[2]:match[2]+k] for \
        match in longest_match]

    matching_query = "".join(matching_query)
    pos = subject.index(matching_query)

    pos_match_in_query = query.index(matching_query)

    print(subject)
    print(" " * (pos - pos_match_in_query) + "." * pos_match_in_query \
                + "|" * len(matching_query))
    print(" " * (pos - pos_match_in_query) + query)


def rev_comp(seq):
    """Returns the reverse complement of a DNA sequence

        Keyword Arguments:
            seq -- string, DNA sequence
    """
    comp_dic = {"A": "T", "C": "G",
                "T": "A", "G": "C"}
    return "".join([comp_dic[x] for x in seq[::-1]])
    
def SSAHA(query, seqs, k = 2, print_info = True):
    """Performs the SSAHA on a DNA sequence query and database

        Keyword Arguments:
            query -- string or list, containing the DNA sequence(s)
            seqs -- list, containing the DNA sequences to serve as database
            k -- int, size of the k-mers
            print_info -- bool, if yes detailed information such as the hash
                            table is printed to the console. Not suitable for
                            large sequences.
    """
    hashtable = create_hashtable(seqs, k)
    if type(query) == list: #This corresponds with the TAIR database search
        for i, Q in enumerate(query):
            M_list = sequence_search(Q, hashtable, k)
            longest_match = get_longest_match(M_list, len(seqs))
            if longest_match:
                print("Longest match for query {}:".format(
                        str(i)))
                print(longest_match)

        print("Now for the reverse complement")
        for i, Q in enumerate(query):
            Q = rev_comp(Q)
            M_list = sequence_search(Q, hashtable, k)
            longest_match = get_longest_match(M_list, len(seqs))
            if longest_match:
                print("Longest match for query's {} rev_comp :".format(
                        str(i)))
                print(longest_match)

    else: #Happens if the query is a string (for the test query)
        M_list = sequence_search(query, hashtable, k)
        longest_match = get_longest_match(M_list, len(seqs))

        print("The hashtable contains {} k-mers of size {}".format(
                            str(len(hashtable.keys())), str(k)))
        if print_info:
            print("HashTable: ")
            pretty_print_table(hashtable)
        
            print("Number of hits: {}".format(str(len(M_list))))
            print("First hit: {}".format(M_list[0]))
            print("Last hit: {}".format(M_list[-1]))
    
        #Visualize the alignment
        if longest_match:
            print_alignment(query,longest_match, seqs, k)

if __name__ == "__main__": 
    if len(argv) != 3:
        print("Usage: python assignment2_joris.py <genome.fa> <query.fa>")
        exit(2)

    ara_genome_file = argv[1]
    query_file = argv[2] 

    query = "TGCAACAT"
    s1 = "GTGACGTCACTCTGAGGATCCCCTGGGTGTGG"
    s2 = "GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT"
    s3 = "GGATCCCCTGTCCTCTCTGTCACATA"
    seqs = [s1,s2,s3]

    print("Using the test database:")
    SSAHA(query, seqs, 3, print_info = True)

    
    print("\n\nUsing the Tair database:\n")
    queries = parse_fasta(query_file)
    genome = parse_fasta(ara_genome_file)
    l_seq = sum([len(x) for x in genome])


    print("The genome contains: {} sequences".format(
            str(len(genome))))
    print("The total sequence length is: {} bp".format(
            str(l_seq)))
    SSAHA(queries, genome, 15, print_info = False)

