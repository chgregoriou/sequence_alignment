#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the Dynamic Programming assignment in the Algorithms in Sequence Analysis course at the VU.
    
INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    Gregoriou Charalambos 2740346
"""


import argparse
import pickle


def parse_args():
    """Parses inputs from commandline and returns them as a Namespace object."""

    parser = argparse.ArgumentParser(prog='python3 align.py', formatter_class=argparse.RawTextHelpFormatter,
                                     description='  Aligns the first two sequences in a specified FASTA\n'
                                                 '  file with a chosen strategy and parameters.\n'
                                                 '\ndefaults:\n  strategy = global\n  substitution matrix = pam250\n'
                                                 '  gap penalty = 2')
        
    parser.add_argument('fasta', help='path to a FASTA formatted input file')
    parser.add_argument('output', nargs='*', help='path to an output file where the alignment is saved\n'
                                                  '  (if a second output file is given,\n'
                                                  '   save the score matrix in there)')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='print the score matrix and alignment on screen', default=False)
    parser.add_argument('-s', '--strategy', dest='strategy',
                        choices=['global', 'semiglobal', 'local'], default="global")
    parser.add_argument('-m', '--matrix', dest='substitution_matrix',
                        choices=['pam250', 'blosum62', 'identity'], default='pam250')
    parser.add_argument('-g', '--gap_penalty', dest='gap_penalty', type=int,
                        help='must be a positive integer', default=2)

    args = parser.parse_args()

    args.align_out = args.output[0] if args.output else False
    args.matrix_out = args.output[1] if len(args.output) >= 2 else False
    # Fancy inline if-else statements. Use cautiously!
                      
    if args.gap_penalty <= 0:
        parser.error('gap penalty must be a positive integer')

    return args


def load_substitution_matrix(name):
    """Loads and returns the specified substitution matrix from a pickle (.pkl) file."""
    # Substitution matrices have been prepared as nested dictionaries:
    # the score of substituting A for Z can be found with subst['A']['Z']
    # NOTE: Only works if working directory contains the correct folder and file!
    
    with open('substitution_matrices/%s.pkl' % name, 'rb') as f:
        subst = pickle.load(f)
    return subst
    

def load_sequences(filepath):
    """Reads a FASTA file and returns the first two sequences it contains."""
    
    seq1 = []
    seq2 = []
    with open(filepath,'r') as f:
        for line in f:
            if line.startswith('>'):
                if not seq1:
                    current_seq = seq1
                elif not seq2:
                    current_seq = seq2
                else:
                    break # Stop if a 3rd sequence is encountered
            else:
                current_seq.append(line.strip())
    
    if not seq2:
        raise Exception('Error: Not enough sequences in specified FASTA file.')
    
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)
    return seq1, seq2


def align(seq1, seq2, strategy, substitution_matrix, gap_penalty):
    """Do pairwise alignment using the specified strategy and parameters."""
    # This function consists of 3 parts:
    #
    #   1) Initialize a score matrix as a "list of lists" of the appropriate length.
    #      Fill in the correct values for the first row and column given the strategy.
    #        (local / semiglobal = 0  --  global = stacking gap penalties)
    #   2) Fill in the rest of the score matrix using Dynamic Programming, accounting
    #      for the selected alignment strategy, substitution matrix and gap penalty.
    #   3) Perform the correct traceback routine on your filled in score matrix.
    #
    # Both the resulting alignment (sequences with gaps and the corresponding score)
    # and the filled in score matrix are returned as outputs.
    #
    # NOTE: You are strongly encouraged to think about how you can reuse (parts of)
    #       your code between steps 2 and 3 for the different strategies!

    ### 1: Initialize
    M = len(seq1)+1
    N = len(seq2)+1
    score_matrix = []
    for i in range(M):
        row = []
        score_matrix.append(row)
        for j in range(N):
            row.append(0)
    
    if strategy == 'global':
        # for global alignment the first row and column are initialized with the sum of the gap penalty
        for i in range(1, N):
            score_matrix[0][i] = score_matrix[0][i-1] - gap_penalty
        for i in range(1, M):
            score_matrix[i][0] = score_matrix[i-1][0] - gap_penalty
        #####################
        #  END CODING HERE  #
        #####################

    ### 2: Fill in Score Matrix
 

    def dp_function(diagonal, top, left, aa1, aa2):
        """
        This function takes three values and the amino acid of each sequence in a specific position.
        Computes the three values and takes the max and assigns it to the cell
        """
        # map list for keeping track from where the max value was taken
        map_list = [(1, 1), (1, 0), (0, 1)]
        maximum = top-gap_penalty
        # in case coming from top
        trace_value = map_list[1]
        # find the maximum computed value. in case of equality we always prioritize coming from the top cell
        # so that leads us to the high road
        if maximum < diagonal+substitution_matrix[aa1][aa2]:
            maximum = diagonal+substitution_matrix[aa1][aa2]
            # in case of coming from the diagonal
            trace_value = map_list[0]

        if maximum < left-gap_penalty:
            maximum = left-gap_penalty
            # in case coming from left
            trace_value = map_list[2]

        # in case of local alignment if a negative value is computed we reassign the value to 0
        if strategy == "local" and maximum < 0:
            maximum = 0
            # if we get the value from this case then we stop the alignment
            trace_value = "stop"
        # this list has the steps taken to each cell in one row
        row_trace.append(trace_value)
        return maximum

    # initialize the trace list
    # the first row always comes form the left
    # thi list contains the steps taken to each cell for all the matrix
    trace_list = [[(0, 1)]*N]
    for i in range(1, M):
        # in the first column always comes from the top
        row_trace = [(1, 0)]
        for j in range(1, N):
            # for each cell compute its value
            score_matrix[i][j] = dp_function(score_matrix[i-1][j-1], score_matrix[i-1][j], score_matrix[i][j-1],
                                             seq1[i-1], seq2[j-1])
        # append the steps of each row to the list
        trace_list.append(row_trace)
            

    ### 3: Traceback
    
    aligned_seq1 = ''  # These are dummy values! Change the code so that
    aligned_seq2 = ''  # aligned_seq1 and _seq2 contain the input sequences
    # with gaps inserted at the appropriate positions.

    if strategy == "global":
        # for global alignment the alignment score is the last cell of the matrix
        align_score = score_matrix[-1][-1]
        # we start the traceback from the cell that contains the alignment score
        row = len(trace_list)-1
        column = len(trace_list[0])-1
        # as long as we do not reach the start of the matrix we have to keep going back
        while row > 0 or column > 0:
            # take the cell that we came here
            step = trace_list[row][column]
            # if one of the steps is equal to 0 then we have to append "-" gap to one sequence
            # otherwise we append in the the aligned sequences the corresponding amino acids
            if step[1] == 0:
                aligned_seq2 += "-"
                aligned_seq1 += seq1[row-1]
            elif step[0] == 0:
                aligned_seq1 += "-"
                aligned_seq2 += seq2[column-1]
            else:
                aligned_seq1 += seq1[row-1]
                aligned_seq2 += seq2[column-1]
            # subtract from the indices the corresponding value row, column
            row -= step[0]
            column -= step[1]
    elif strategy == "semiglobal":
        # if the strategy is semi-global the alignment score is the maximum value among the cells
        # of the last column and last row. In case of equality we take the cell with the lowest row index
        # and highest column index
        align_score = score_matrix[0][-1]
        row = 0
        column = N - 1
        # locating the alignment score
        for i in range(1, M):
            if score_matrix[i][-1] > align_score:
                align_score = score_matrix[i][-1]
                row = i
        for i in range(N-1, 0, -1):
            if score_matrix[-1][i] > align_score:
                align_score = score_matrix[-1][i]
                column = i
                row = M - 1
        # Initializing the alignment sequences with the the values that have been left behind
        # since we don't start the traceback from the last cell of the matrix
        for i in range(N-1, column, -1):
            aligned_seq1 += "-"
            aligned_seq2 += seq2[i-1]
        for i in range(M-1, row, -1):
            aligned_seq2 += "-"
            aligned_seq1 += seq1[i-1]
        # loop through the steps taken according to the traceback list until the first cell is reached
        while row > 0 or column > 0:
            # the step taken to end up on this cell
            # one of (1,1): from diagonal, (0,1): from left, (1,0): from top
            step = trace_list[row][column]
            # if the second item of the tuple is 0 the we stay in the same column
            # therefore appending a gap into sequence 2 (horizontal)
            if step[1] == 0:
                aligned_seq2 += "-"
                aligned_seq1 += seq1[row-1]
            # if the first item of the tuple is 0 then we stay in the same row
            # therefore appending a gap into sequence 1 (vertical)
            elif step[0] == 0:
                aligned_seq1 += "-"
                aligned_seq2 += seq2[column-1]
            # If both items in the tuple are one then we move one column left and one row up
            # therefore we append the letters accordingly
            else:
                aligned_seq1 += seq1[row-1]
                aligned_seq2 += seq2[column-1]
            # subtracting from row the value of item 1, and form the column the value of item 2 in the tuple
            row -= step[0]
            column -= step[1]
    elif strategy == "local":
        # In the local alignment the alignment score is the maximum value of the matrix
        # in case of equality we take the rightmost cell with the lowest row index
        align_score = score_matrix[1][1]
        row = 1
        column = 1
        for i in range(1, M):
            for j in range(1, N):
                if score_matrix[i][j] > align_score or (j > column and score_matrix[i][j] >= align_score):
                    align_score = score_matrix[i][j]
                    row = i
                    column = j
        # the local alignment stops when a zero is reached that came from the re-assignment
        # or when the start of the matrix is reached
        while row > 0 and column > 0:
            step = trace_list[row][column]
            if step == "stop":
                break
            if step[1] == 0:
                aligned_seq2 += "-"
                aligned_seq1 += seq1[row-1]
            elif step[0] == 0:
                aligned_seq1 += "-"
                aligned_seq2 += seq2[column-1]
            else:
                aligned_seq1 += seq1[row-1]
                aligned_seq2 += seq2[column-1]
            row -= step[0]
            column -= step[1]

    # reverse the alignments
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    alignment = (aligned_seq1, aligned_seq2, align_score)
    return (alignment, score_matrix)


def print_score_matrix(s1,s2,mat):
    "Pretty print function for a score matrix."
    
    # Prepend filler characters to seq1 and seq2
    s1 = '-' + s1
    s2 = ' -' + s2
    
    # Print them around the score matrix, in columns of 5 characters
    print(''.join(['%5s' % aa for aa in s2])) # Convert s2 to a list of length 5 strings, then join it back into a string
    for i, row in enumerate(mat):               # Iterate through the rows of your score matrix (and keep count with 'i').
        vals = ['%5i' % val for val in row]    # Convert this row's scores to a list of strings.
        vals.insert(0,'%5s' % s1[i])           # Add this row's character from s2 to the front of the list
        print(''.join(vals))                   # Join the list elements into a single string, and print the line.


def print_alignment(a):
    "Pretty print function for an alignment (and alignment score)."
    
    # Unpack the alignment tuple
    seq1 = a[0]
    seq2 = a[1]
    score = a[2]
    
    # Check which positions are identical
    match = ''
    for i in range(len(seq1)): # Remember: Aligned sequences have the same length!
        match += '|' if seq1[i] == seq2[i] else ' ' # Fancy inline if-else statement. Use cautiously!
            
    # Concatenate lines into a list, and join them together with newline characters.
    print('\n'.join([seq1,match,seq2,'','Score = %i' % score]))


def save_alignment(a,f):
    "Saves two aligned sequences and their alignment score to a file."
    with open(f,'w') as out:
        out.write(a[0] + '\n') # Aligned sequence 1
        out.write(a[1] + '\n') # Aligned sequence 2
        out.write('Score: %i' % a[2]) # Alignment score


def save_score_matrix(m,f):
    "Saves a score matrix to a file in tab-separated format."
    with open(f,'w') as out:
        for row in m:
            vals = [str(val) for val in row]
            out.write('\t'.join(vals)+'\n')
    

def main(args=False):
    # Process arguments and load required data
    if not args:
        args = parse_args()
    
    sub_mat = load_substitution_matrix(args.substitution_matrix)
    seq1, seq2 = load_sequences(args.fasta)

    # Perform specified alignment
    strat = args.strategy
    gp = args.gap_penalty
    alignment, score_matrix = align(seq1, seq2, strat, sub_mat, gp)

    # If running in "verbose" mode, print additional output
    if args.verbose:
        print_score_matrix(seq1,seq2,score_matrix)
        print('') # Insert a blank line in between
        print_alignment(alignment)
    
    # Save results
    if args.align_out: save_alignment(alignment, args.align_out)
    if args.matrix_out: save_score_matrix(score_matrix, args.matrix_out)


if __name__ == '__main__':
    main()
