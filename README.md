# sequence_alignment
Performs pairwise alignment of protein sequences using dynamic programming based on substitution matrices:
1. pam250
2.  blosum62
3.  identity

Arguments:
1. path to fasta file that contains two sequenceses
2. path to an output text file
3. -v, --verbose (whether to print output on the terminal)
4. -s, --strategy (global, semiglobal, local) - alignemnt method - default: global
5. -m, --matrix (pam250, blosum62, identity) - substitution matrix to use - default: pam250
6. -g, --gap_penalty (positive integer penalty to use) - default: 2

Run through the comand line:

`python align.py PATH.fasta output.txt`
