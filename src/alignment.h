/*
 * This file contains the interface for an Alignment object that aligns two
 * provided sequences.
 *
 * Author: Dan Lin, danielin@uw.edu
 *
 */

#include <string>
#include <vector>

class Alignment {
 public:
  const std::string& sequence_str1() { return sequence_1_.sequence_str; }
  const std::string& sequence_str2() { return sequence_2_.sequence_str; }

  void set_sequence_str1(const std::string& str) {
    sequence_1_.sequence_str = str;
  }
  void set_sequence_str2(const std::string& str) {
    sequence_2_.sequence_str = str;
  }
  // Constructor. Initializes using the contents at 'file_path_1' and
  // 'file_path_2' which should each contain a sequence in FASTA format.
  Alignment(const std::string& file_path_1,
            const std::string& file_path_2);

  // Get the best local alignment using the Smith-Waterman algorithm. If
  // 'print' is true we will print the local alignment to stdout.
  int GetLocalAlignment(bool print = true) const;

  // Get the best global alignment using the Needleman-Wunsch algorithm. If
  // 'print' is true we will print the global alignment to stdout.
  int GetGlobalAlignment(bool print = true) const;


  // Prints both seqences of this alignment object.
  void PrintSequences() const;

 private:
  // Struct that encapsulates all of the state associated with a given
  // sequence.
  struct Sequence {
    // Indentifier string for this particular sequence.
    std::string identifier;

    // The sequence itself in string form.
    std::string sequence_str;
  };

  // Open the FASTA file at the 'sequence_file' path and place the amino acid
  // sequence into 'sequence_string'.
  void ReadSequenceFromFile(const std::string& sequence_file,
                            Sequence *sequence_string);

  // Print 'str' with columns to stdout.
  void PrintWithColumns(const std::string& str) const;

  // Print 'matrix' with columns to stdout.
  void PrintMatrix(const std::vector<std::vector<int>>& matrix) const;

  // Prints the traceback found in 'align_matrix' to stdout. Traverses from
  // 'start_row' and 'start_col' using 'trace_matrix' for backtracking.
  void PrintTraceback(const std::vector<std::vector<int>>& align_matrix,
                      const std::vector<std::vector<char>>& trace_matrix,
                      int start_row,
                      int start_col,
                      bool is_global) const;

  // Print the prefix that comes before a printed sequence using 'identifier'
  // for the sequence and the 'index' of the next aligned char in the sequence.
  void PrintLinePrefix(const std::string& identifier, int index) const;

  // Helper function that retrieves the score for the file in 'align_matrix'
  // at 'row' and 'col'. 'allow_reset' is true if we allow alignments to be
  // restarted from from fresh midway while comparing strings. 'score' will be
  // filled with the besc score found and 'trace' will he filled with a
  // character indicating from which tile we traversed from to the returned
  // score.
  // For trace:
  // x: Denotes that we are traversing by matching the two characters.
  // l: Denotes that we are traversing by adding a gap to the second sequence.
  // u: Denotes that we are traversing by adding a gap to the first sequence.
  // r: Denotes that we are resetting the string as the score is below zero.
  void GetNextScore(const std::vector<std::vector<int>>& align_matrix,
                    int row, int col,
                    bool allow_reset,
                    int *score,
                    char *trace) const;

  // Creates and returns an alignment matrix based on the two sequences.
  std::vector<std::vector<int>> CreateAlignmentMatrix() const;

  // Creates and returns a trace matrix based on the two sequences.
  std::vector<std::vector<char>> CreateTracebackMatrix() const;

 private:
  // First sequence that will be used for alignment.
  Sequence sequence_1_;

  // Second sequence that will be used for alignment.
  Sequence sequence_2_;
};
