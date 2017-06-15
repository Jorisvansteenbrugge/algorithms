#!/usr/bin/env python

"""
Author: Joris van Steenbrugge

Description: Train a Hidden Markov Model and sample possible sequences 
 based on the probabilities.
"""
from __future__ import division
from sys import argv
from random import random
from Bio import SeqIO
import numpy

# Background amino acid probabilities
pa = { 'A': 0.074, 'C': 0.025, 'D': 0.054, 'E': 0.054, 'F': 0.047, 'G': 0.074,
       'H': 0.026, 'I': 0.068, 'L': 0.099, 'K': 0.058, 'M': 0.025, 'N': 0.045,
       'P': 0.039, 'Q': 0.034, 'R': 0.052, 'S': 0.057, 'T': 0.051, 'V': 0.073,
       'W': 0.013, 'Y': 0.034 } 


class HMM():
    """HMM object to store an HMM model

    This object is designed to keep track of all HMM states, emissions, and 
    transitions. It may be used in your implementation, but may also be ignored,
    and replaced by a data structure of choice
    """
    # Emission probabilities for the match and insert states
    e_m   = []; e_i   = pa; 
    
    # Transition probabilities from/to matches, inserts and deletions
    t_mm  = []; t_mi  = []; t_md = [];
    t_im  = []; t_ii  = []
    t_dm  = []; t_dd  = []; 
    
    def __init__(self, match_states, insert_states):
        """Initialize HMM object with number of match states
        nmatches: int, number of match states
        """

        self.match_states_idx  = match_states
        self.insert_states_idx = insert_states
        self.nmatches = len(match_states)

        print("There are {} match states".format(self.nmatches))

        self.e_m   = [dict(pa) for _ in range(0, self.nmatches)]
        for i in range(0, self.nmatches):
            for j in pa.keys():
                self.e_m[i][j] = 0.0
        self.e_i   = pa;

        self.t_mm  = [0.0] * (self.nmatches + 1)
        self.t_mi  = [0.0] * (self.nmatches + 1)
        self.t_im  = [1.0] * (self.nmatches + 1)
        self.t_ii  = [0.0] * (self.nmatches + 1)
        self.t_md  = [0.0] * (self.nmatches + 1)
        self.t_dm  = [1.0] * (self.nmatches + 1)
        self.t_dd  = [0.0] * (self.nmatches + 1)

    def get_probabilities(self, sequences):
        """Wrapper to calculate the transmissionand emission probabilities

            Keyword Arguments:
                sequences -- list of strings, aligned amino acid sequences
        """
        for pos, match_idx in enumerate(self.match_states_idx):
            for amino in pa.keys(): 
                self.e_m[pos][amino] = self._emissions_match_state(amino, 
                                                sequences, match_idx)
        self._transmission_probabilities(sequences)

    def print_probabilities(self):
        """Prints the different transmission probabilities to the console.
        """
        print("Transmission:")
        print("mm {}".format(tuple(self.t_mm)))
        print("mi {}".format(tuple(self.t_mi)))
        print("md {}".format(tuple(self.t_md)))
        print("im {}".format(tuple(self.t_im)))
        print("ii {}".format(tuple(self.t_ii)))
        print("dm {}".format(tuple(self.t_dm)))
        print("dd {}".format(tuple(self.t_dd)))

    def get_emission_matrix(self):
        """Yields an emission probability matrix for a position specific profile.

            Probabilities are rounded to 3 decimals.
        """

        #Header line with all amino acids
        aminos = self.e_m[0].keys()
        print("\t{}".format("\t".join(aminos)))

        #for each match state take the emission state
        for idx, match_state in enumerate(self.e_m):
            #Join the emission state dictionary to a tab separated string
            # containing rounded probability scores.
            values = "\t".join(["%.3f" % round(x[1],3) \
                            for x in match_state.items()])
            yield("{}\t{}".format(idx+1, values))


    def _emissions_match_state(self, amino, sequences, match_idx):
        """Returns the emission probabilities in each match state.

            This function is intended to be called from within the
            class object privately.-=p

            .

            Keyword Arguments:
                amino -- char, the one letter amino acid to count the 
                            emission probability for.
                sequences -- list of strings, containint the aligned aa 
                                sequences
                match_idx -- int, th position of the current match state
        """

        total = [seq[match_idx] for seq in sequences if seq[match_idx] != "-"]
        total = len(total)
        count = 0
        for seq in sequences:
            if seq[match_idx] == amino:
                count += 1
        return count / total 

    def _pos_prob(self, sequences, pos):
        """Returns four transmission probabilities at the start or end position

            Keyword Arguments:
                sequences -- list of strings, aligned amino acid sequences
                pos -- int, should be either 0 for start or -1 for end position

            As the start and end positions have no predecessors and successors,
            respectively , the transmission probabilities are calculated here 
            by calculating the chance for a match/deletion/insertion at the 
            position itself.
        """

        #Get all sequences at the specific position
        pos_seq = [seq[pos] for seq in sequences]
        length = len(pos_seq)
        
        #Either a match or insert state
        pos_state = self._get_state(0)

        mm = 0.0
        mi = 0.0
        md = 0.0
        im = 0.0

        c = pos_seq.count("-")
        # If the position is a match state the number of - indicates the number
        # of deletions.
        if pos_state == 'match':
            md = (c / length)
            mm = ((length - c) / length)

        #If the position is a insert state the number of "-" indicates the 
        # number of m -> positions.
        else: # state == "insert"
            mm = (c / length)
            mi = ((length - c) / length)


        return mm, mi, md, im

    def _get_state(self, idx):
        """Returns what state a position is in
    
            Keyword Arguments:
                idx -- int, the position
        """
        if idx in self.match_states_idx:
            return 'match'
        else:
            return 'insert'

    def _transmission_probabilities(self, sequences):
        """Calculates the transmission probabilities.

            Keyword Arguments:
                sequences -- list of strings, aligned amino acid sequences

            Calculates probabilities for various transmissions (mm,md, mi) and
            stores them in class contained lists.
            
            Transitions that are not calculated: im, dd, ii due to time 
            restrictions.

        """
        self.t_mm[0], self.t_mi[0], \
        self.t_md[0], self.t_im[0] = self._pos_prob(sequences, 0)

        self.t_mm[-1], self.t_mi[-1], \
        self.t_md[-1], self.t_im[-1] = self._pos_prob(sequences, -1)

        for i, M in enumerate(self.match_states_idx):
            next_state  = self._get_state(M + 1)
            current_seq = [seq[M] for seq in sequences]

            try: # The last match state is possibly the last state
                next_seq = [seq[M + 1] for seq in sequences]
                length   = len(next_seq)

                if next_state == "match":
                    c = next_seq.count("-")
                    mm = ((length - c) / length)
                    md = (c / length)

                    self.t_mm[i + 1] = mm
                    self.t_md[i + 1] = md
                    self.t_mi[i + 1] = 0.0
                else: #next_state == "insert"
                    c  = next_seq.count("-")
                    mi = ((length - c) / length)
                    mm = (c / length)

                    self.t_mi[i + 1] = mi
                    self.t_md[i + 1] = 0.0
                    self.t_mm[i + 1] = mm
            except IndexError:
                pass

def sample(events):
    """Return a key from dict based on the probabilities 

    events: dict of {key: probability}, sum of probabilities
    should be 1.0. 
    """
    k = events.keys()
    cum = [0 for i in k]

    cum[0] = events[k[0]]
    for i in range(1,len(events)):
        cum[i] = cum[i-1] + events[k[i]]
    # Should not be necessary, but for safety
    cum[len(cum)-1] = 1.0

    r = random()
    pick = ''
    i = 0
    while (pick == '') and (i < len(cum)):
        if r < cum[i]:
            pick = k[i]
        i = i + 1
    return pick
    
def get_states(sequences):
    """Returns the indices of match states.

        Keyword Arguments:
            sequences -- list of strings, aligned amino acid sequences.

        This is implemented in a way that the length of the indices list is 
        the number of match states.
    
    """
    match_states = []
    insert_states = []
    half_length = len(sequences) / 2.0

    for i in range(len(sequences[0])):
        has_amino = 0
        for seq in sequences:
            if seq[i] in pa:
                has_amino += 1
        if has_amino >= half_length:
            match_states.append(i)
        else: 
            insert_states.append(i)

    return match_states, insert_states
    
def calc_traverse_path(a, p):
    """Returns a randomly picked value in a based on probabilities in p

        Keyword Arguments:
            a -- list, containing values to be considered
            p -- list, containing probabilites for each considered value
    """
    return numpy.random.choice(a = a, p = p)

def traverse(hmm):
    """Returns a randomly selected sequence based on probabilities in the HMM

        Keyword Arguments:
            hmm -- hmm class object

        iteratively randomly selects a transition based on probabilities.
        In match states an amino acid is selected based on the match state
        emission probability list (contained in the hmm class).
        In insert states an amino acid is selected based on the background
        distribution emissions.
    """
    sequence = ""
    state = "match"
    for idx in range(len(hmm.t_mm)):
        

        if state == "match":
            if idx != 0:
                sequence += sample(hmm.e_m[idx - 1])

            mm = hmm.t_mm[ idx ]
            md = hmm.t_md[ idx ]
            mi = hmm.t_mi[ idx ]

            a = ("mm", "md", "mi")
            p = ( mm,   md,   mi )
            
        elif state == "insert":
            if idx != 0:
                sequence += sample(hmm.e_i)

            im = hmm.t_im[ idx ]
            ii = hmm.t_ii[ idx ]

            a = ("im", "ii")
            p = ( im,   ii )

        elif state == "delete":
            dd = hmm.t_dd[ idx ]
            dm = hmm.t_dm[ idx ]

            a = ("dd", "dm")
            p = ( dd,   dm )

        choice = calc_traverse_path(a, p)

        if choice == "mm":
            state = "match"
        elif choice == "md":
            state = "delete"
        elif choice == "mi":
            state = "insert"
        elif choice == "im":
            state = "match"
        elif choice == "ii":
            state = "insert"
        elif choice == "dd":
            state = "delete"
        elif choice == "dm":
            state = "match"
    return sequence

def parse_fasta(file_name):
    """Returns a list with sequences from a fasta file

        It still feels a bit redundant to 

        Keyword Arguments:
            file_name -- string, file path of the fasta file
    """
    sequences = []
    seq = []
    with open(file_name) as in_file:
        for line in in_file:
            if line.startswith(">"): #entry start
                if len(seq) == 0:
                    pass
                else:#entry is present
                    sequences.append("".join(seq))
                    seq = []
            else:
                seq.append(line.strip())
    sequences.append("".join(seq))
    return sequences
    #return [str(entry.seq) for entry in SeqIO.parse(file_name, "fasta")]
           

def find_HMM(infile):
    """Wrapper function to train the HMM.

        Keyword Arguments:
            infile -- filepath to a fasta file containin aligned amino acid
                        sequences.
    """
    sequences = parse_fasta(infile)
    match_states, insert_states = get_states(sequences)
    hmm = HMM(match_states, insert_states)
    hmm.get_probabilities(sequences)

    return hmm

if __name__ == "__main__":
    #Q1
    hmm = find_HMM(infile = "/home/joris/hmm_data/test.fasta")
    #Q2
    hmm.print_probabilities()
    #Q3
    print("\n".join(list(hmm.get_emission_matrix())))

    #Q4
    for i in range(10):
        print(traverse(hmm))

    # #Q6
    hmm = find_HMM(infile = "/home/joris/hmm_data/test_large.fasta")
    hmm.print_probabilities()
    print("\n".join(hmm.get_emission_matrix()))
    print(traverse(hmm))