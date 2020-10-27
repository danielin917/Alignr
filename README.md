# Alignr

Alignr performs both global and local alignment on amino acid sequence files in FASTA format based on the Needleman-Wunsch and Smith-Waterman algorithms. It will also calculate the empirical P-Value score of these alignments by performing random shuffles and aligning against these shuffled sequences.

The following command will build the executable. C++11 features are used so c++0x is a required flag for the compiler.                                       
                                                                                 
g++ -std=c++0x main.cc src/alignment.cc                                              
                                                                                                                                                                    
To run the code after building execute using the following args:                 
./a.out <local | global> file1_path file2_path <optional: num_permutations>   
                                                                                 
By default after doing intial alignment permutation runs of the second sequence will kick off in order to calculate the P-Value. If you would like to skip the permutation step simply add 0 to the end and the permutation step will be skipped. Code loosely requires valid amino acid chars and a FASTA like format to pick up an identifier.   

By default this code is hardcoded to use the BLOSUM62 score matrix.

TODO:
- Allow pass in of score matrix file
- Incorporate cmake
- Allow for DNA alignment/arbitrary text using passed in score matrx
