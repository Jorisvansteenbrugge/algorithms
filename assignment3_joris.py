#!/usr/bin/env python

"""
Author: 

Description: this is a script to ...
"""
from __future__ import division
from sys import argv
from random import random
from Bio import SeqIO

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

        self.e_m   = [dict(pa) for _ in range(0, self.nmatches)]
        for i in range(0, self.nmatches):
            for j in pa.keys():
                self.e_m[i][j] = 0.0
        self.e_i   = pa;

        self.t_mm  = [0.0] * (self.nmatches + 1)
        self.t_mi  = [0.0] * (self.nmatches + 1)
        self.t_im  = [0.0] * (self.nmatches + 1)
        self.t_ii  = [0.0] * (self.nmatches + 1)
        self.t_md  = [0.0] * (self.nmatches + 1)
        self.t_dm  = [0.0] * (self.nmatches + 1)
        self.t_dd  = [0.0] * (self.nmatches + 1)

    def get_probabilities(self, sequences):
        for pos, match_idx in enumerate(self.match_states_idx):
            for amino in pa.keys(): 
                self.e_m[pos][amino] = self._emissions_match_state(amino, 
                                                sequences, match_idx)
        self._transmission_probabilities(sequences)

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
        pos_seq = [seq[pos] for seq in sequences]
        length = len(pos_seq)
        
        pos_state = self._get_state(0)

        mm = 0.0
        mi = 0.0
        md = 0.0
        im = 0.0

        c = pos_seq.count("-")
        if pos_state == 'match':
            md = (c / length)
            mm = ((length - c) / length)
        else: # state == "insert"
            mm = (c / length)
            mi = ((length - c) / length)


        return mm, mi, md, im

    def _get_state(self, idx):
        if idx in self.match_states_idx:
            return 'match'
        else:
            return 'insert'


    def _transmission_probabilities(self, sequences):
        """
        Missing: im, dd, ii 

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



def parse_fasta(file_name):
    """Returns a list with sequences from a fasta file

        The fasta file is parsed using Biopython as fasta parsing is
        trivial for this assignment.

        Keyword Arguments:
            file_name -- string, file path of the fasta file
    """
    return [str(entry.seq) for entry in SeqIO.parse(file_name, "fasta")]
           
if __name__ == "__main__":

    # implement main code here
    infile = "/home/joris/algorithms/test.fasta"
    sequences = parse_fasta(infile)

    match_states, insert_states = get_states(sequences)

    hmm = HMM(match_states, insert_states)
    hmm.get_probabilities(sequences)
