
from __future__ import division
import argparse 

blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-seqA", dest = "seqA", required = True,
        help = "The first protein sequence")
    parser.add_argument("-seqB", dest = "seqB", required = True,
        help = "The second protein sequence")
    parser.add_argument("-gap_penalty", dest = "gap_penalty", required = False,
        default = 4, type = int,
        help = "The scoring penalty for introducing gaps in the alignment")
    parser.add_argument("-end_gap_penalty", dest = "end_gap_penalty", required = False,
        default = 2, type = int,    
        help = "The scoring penalty for introducing gaps at the start or beginning of the alignment")

    return parser.parse_args()

def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int,parts[1:])))
    return order, blosum_matrix

BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()

def score(res1, res2):
    """Return similarity score from BLOSUM62 matrix for two residues
    
    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]

def score_at_position(i, j, sequences):
    res_a = sequences[0][i - 1]
    res_b = sequences[1][j - 1]
    # print("{} -> {}".format(res_a, res_b))
    scr = score(res_a, res_b)
    return scr

def gap_penalty(n, E):
    return ((-n) * E)

def create_empty_matrix(sequences, Nones= True):      
    m = len(sequences[0]) + 1
    n = len(sequences[1]) + 1
    if Nones:
        return [[None for i in range(n)] for j in range(m)]
    else:
        return [ [ list()  for i in range(n)] for j in range(m)]

def matrix_template(sequences, A = 2):
    """Fill the matrix's first row and first collumn with gap penalties.

        Keyword Arguments:
            sequences -- A tuple reversely sorted on length
                containing two protein sequences 
            A -- End gap penalty score

        Returns:
            A matrix implemented as two dimensional list, 
                to be queried as matrix[row][col]
    """
    #Create an empty matrix with the corresponding dimensions
    matrix = create_empty_matrix(sequences)
    matrix[0][0] = 0

    #End/start gap penalties in each column
    for i in range(1, len(matrix[0])):
        matrix[0][i] = gap_penalty(i, A)

    #End/start gap penalties in each row
    for j in range(1, len(matrix)):
        matrix[j][0] = gap_penalty(j, A)
    return matrix

def parent_matrix_template(sequences):
    """Fill the parent matrix's first row and first collumn.
    The positions are tracing back to position (0,0) so that traceback 
    function will retrain the whole sequences in the alingment.

        Keyword Arguments:
            sequences -- A tuple reversely sorted on length
                         containing two protein sequences 
        Returns:
            A matrix implemented as two dimensional list,
                to be queried as matrix[row][col]
    """
    parent_matrix = create_empty_matrix(sequences, False)
    parent_matrix[0][0] = (0,0)
    for i in range(1, len(parent_matrix[0])):
        parent_matrix[0][i] = (0,i -1)
    for j in range(1, len(parent_matrix)):
        parent_matrix[j][0] = (j - 1, 0)

    return parent_matrix



def pretty_print(matrix):
    for line in matrix:
        print "\t".join(map(str, line))



def fill_matrix(matrix, parent_matrix, sequences, E = 4, A = 2):

    for i in range(1, len(matrix)):
        left_gap = True
        for j in range(1, len(matrix[0])):
            up_gap = True

            if left_gap:
                left = matrix[i][j-1] + gap_penalty(1, A)
            else:
                left = matrix[i][j-1] + gap_penalty(1, E)

            if up_gap:
                up = matrix[i-1][j] + gap_penalty(1, A)
            else:
                up = matrix[i-1][j] + gap_penalty(1, E)
            diag = matrix[i-1][j-1] + score_at_position(i, j, sequences)

            idx, value = max(enumerate((diag,left,up)), key = lambda x: x[1])

            if idx == 0:
                parent = (i - 1, j - 1)
                left_gap = False
                up_gap = False
            elif idx == 1:
                parent = (i, j - 1)
                up_gap = False
            elif idx == 2:
                parent = (i - 1, j)
                left_gap = False
            matrix[i][j] = value
            parent_matrix[i][j] = parent
    return matrix, parent_matrix




def traceback(matrix, parent_matrix, sequences):
    i = (len(matrix) - 1)
    j = (len(matrix[0]) - 1)


    alignA = ""
    alignB = ""

    seqA = sequences[1]
    seqB = sequences[0]
    while i != 0 and j != 0:
        next_j = parent_matrix[i][j][1]
        next_i = parent_matrix[i][j][0]

        if j == next_j:
            alignA += "-"
        else:
            alignA += seqA[j-1]

        if i == next_i:
            alignB += "-"
        else:
            alignB += seqB[i-1]

        i = next_i
        j = next_j

    return alignA[::-1], alignB[::-1]

        
def get_identity(seq1, seq2):
    num_ident = len([1 for i in range(len(seq1)) if seq1[i] == seq2[i]])
    
    return (num_ident * 100) / len(seq1)


def align(seq1, seq2, E = 4, A = 2):
    sequences = sorted((seq1, seq2), key = len, reverse = False)
    
        matrix = matrix_template(sequences, A)
        parent_matrix = parent_matrix_template(sequences)
        matrix, parent_matrix = fill_matrix(matrix, parent_matrix, sequences, E, A)
        

        aligned_seq1, aligned_seq2 = traceback(matrix, parent_matrix, sequences)
        alignment_indicator = []
        for i in range(len(aligned_seq2)):
            if aligned_seq1[i] == aligned_seq2[i]:
                alignment_indicator.append('|')
            else:
                alignment_indicator.append(' ')

        print(aligned_seq1)
        print("".join(alignment_indicator))
        print(aligned_seq2)
        print("")
        print("Identity: {:.4f}".format(get_identity(
                                            aligned_seq1, aligned_seq2)
                                          )
             )

    if __name__ == "__main__":
        args = parse_arguments()
        seq1 = "THISLINE"
        seq2 = "ISALIGNED"
        
        align(args.seqA, args.seqB, E = args.gap_penalty,
         A = args.end_gap_penalty)


