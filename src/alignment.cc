/*
 * This file contains the implementation of the Alignment class which uses
 * different algorithms to align two sequences.
 *
 * Author: Dan Lin, danielin@uw.edu
 *
 */

#include <cctype>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "alignment.h"

// Column width that will be used for printing sequences.
static const int kColumnWidth = 60;

// Max dimension for a matrix that will be printed.
static const int kMaxPrintedWidth = 14;

// Score that will be substracted when adding a gap to the alignment.
static const int kGapScore = 4;

using namespace std;

//-----------------------------------------------------------------------------

// Enums representing different amino acids.
enum class AminoAcid {
  // Alanine.
  kA = 0,
  // Arginine.
  kR = 1,
  // Asparagine.
  kN = 2,
  // Aspartic acid.
  kD = 3,
  // Cysteine.
  kC = 4,
  // Glutamine.
  kQ = 5,
  // Glutamic acid.
  kE = 6,
  // Glycine.
  kG = 7,
  // Histidine.
  kH = 8,
  // Isoleucine.
  kI = 9,
  // Leucine.
  kL = 10,
  // Lysine.
  kK = 11,
  // Methionine.
  kM = 12,
  // Phenylalanine.
  kF = 13,
  // Proline.
  kP = 14,
  // Serine.
  kS = 15,
  // Threonine.
  kT = 16,
  // Tryptophan.
  kW = 17,
  // Tyrosine.
  kY = 18,
  // Valine.
  kV = 19,
  // Number of amino acids represented by this enum class.
  kNumAminoAcids = 20
};

//-----------------------------------------------------------------------------

// Blosum62 scoring table that will be used to calculate alignment scores.
static const vector<vector<int>> kBlosum62 =
  {
  /*         A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V */
  /* A */ {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0 },
  /* R */ { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3 },
  /* N */ { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3 },
  /* D */ { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3 },
  /* C */ {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1 },
  /* Q */ { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2 },
  /* E */ { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2 },
  /* G */ {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3 },
  /* H */ { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3 },
  /* I */ { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3 },
  /* L */ { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1 },
  /* K */ { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2 },
  /* M */ { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1 },
  /* F */ { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1 },
  /* P */ { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2 },
  /* S */ {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2 },
  /* T */ {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0 },
  /* W */ { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3 },
  /* Y */ { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1 },
  /* V */ {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4 }
  };

//-----------------------------------------------------------------------------

static AminoAcid CharToAmino(const char c) {
  const char upper_c = toupper(c);
  switch(upper_c) {
    case 'A':
      return AminoAcid::kA;
    case 'R':
      return AminoAcid::kR;
    case 'N':
      return AminoAcid::kN;
    case 'D':
      return AminoAcid::kD;
    case 'C':
      return AminoAcid::kC;
    case 'Q':
      return AminoAcid::kQ;
    case 'E':
      return AminoAcid::kE;
    case 'G':
      return AminoAcid::kG;
    case 'H':
      return AminoAcid::kH;
    case 'I':
      return AminoAcid::kI;
    case 'L':
      return AminoAcid::kL;
    case 'K':
      return AminoAcid::kK;
    case 'M':
      return AminoAcid::kM;
    case 'F':
      return AminoAcid::kF;
    case 'P':
      return AminoAcid::kP;
    case 'S':
      return AminoAcid::kS;
    case 'T':
      return AminoAcid::kT;
    case 'W':
      return AminoAcid::kW;
    case 'Y':
      return AminoAcid::kY;
    case 'V':
      return AminoAcid::kV;
    default:
      cout << upper_c << " is not a recognized amino acid" << endl;
      assert(false);
      break;
  }
  __builtin_unreachable();
}

//-----------------------------------------------------------------------------

Alignment::Alignment(const string& sequence_file_1,
                     const string& sequence_file_2) {
  assert(!sequence_file_1.empty());
  assert(!sequence_file_2.empty());

  ReadSequenceFromFile(sequence_file_1, &sequence_1_);
  ReadSequenceFromFile(sequence_file_2, &sequence_2_);
}

//-----------------------------------------------------------------------------

void Alignment::ReadSequenceFromFile(const string& sequence_file,
                                     Sequence *const sequence) {
  assert(!sequence_file.empty());
  assert(sequence);

  ifstream ifs;
  ifs.open(sequence_file.c_str(), std::ifstream::in);
  if (!ifs.good()) {
    cout << "Unable to open file: " << sequence_file << endl;
    assert(false);
  }

  while (true) {
    char peek_char = ifs.peek();
    if (peek_char == EOF) {
      // All input has been read.
      break;
    }

    string line_str;
    getline(ifs, line_str);
    if (peek_char == '>') {
      if (!sequence->identifier.empty()) {
        continue;
      }

      // Now get the next space or endline.
      int start = line_str.find("|");
      int end = line_str.find("|", start + 1);
      sequence->identifier = line_str.substr(start + 1, end - start - 1);
      continue;
    }
    sequence->sequence_str +=line_str;
  }
}

//-----------------------------------------------------------------------------

int Alignment::GetLocalAlignment(const bool print) const {
  if (print) {
    cout << sequence_1_.identifier << " vs " << sequence_2_.identifier << endl;
  }
  // Initialize an empty alignment matrix.
  vector<vector<int>> align_matrix = CreateAlignmentMatrix();
  vector<vector<char>> trace_matrix = CreateTracebackMatrix();
  const string& sequence_str_1 = sequence_1_.sequence_str;
  const string& sequence_str_2 = sequence_2_.sequence_str;

  int total_max = 0;
  int max_row = 0;
  int max_col = 0;
  for (int ii = 1; ii < sequence_str_1.size() + 1; ++ii) {
    for (int kk = 1; kk < sequence_str_2.size() + 1; ++kk) {
      GetNextScore(align_matrix, ii, kk,
                   true /* allow_reset */,
                   &align_matrix[ii][kk],
                   &trace_matrix[ii][kk]);

      if (align_matrix[ii][kk] > total_max) {
        // Update the global max score for all local alignments seen so far.
        total_max = align_matrix[ii][kk];
        max_row = ii;
        max_col = kk;
      }
    }
  }

  if (!print) {
    return total_max;
  }
  cout << "Max Alignment Score: " << total_max << endl;
  PrintTraceback(align_matrix, trace_matrix, max_row, max_col,
                 false /* is_global */);

  if (align_matrix.size() <= kMaxPrintedWidth &&
      align_matrix.back().size() <= kMaxPrintedWidth) {
    PrintMatrix(align_matrix);
  }
  return total_max;
}

//-----------------------------------------------------------------------------

int Alignment::GetGlobalAlignment(const bool print) const {
  if (print) {
    cout << sequence_1_.identifier << " vs " << sequence_2_.identifier << endl;
  }
  // Initialize an empty alignment matrix.
  vector<vector<int>> align_matrix = CreateAlignmentMatrix();
  vector<vector<char>> trace_matrix = CreateTracebackMatrix();
  const string& sequence_str_1 = sequence_1_.sequence_str;
  const string& sequence_str_2 = sequence_2_.sequence_str;

  for (int ii = 0; ii < sequence_str_1.size() + 1; ++ii) {
    for (int kk = 0; kk < sequence_str_2.size() + 1; ++kk) {
      if (ii == 0 && kk == 0) {
        continue;
      }
      GetNextScore(align_matrix, ii, kk,
                   false /* allow_reset */,
                   &align_matrix[ii][kk],
                   &trace_matrix[ii][kk]);
    }
  }

  if (!print) {
    return align_matrix.back().back();
  }
  cout << "Max Alignment Score: " << align_matrix.back().back() << endl;
  PrintTraceback(align_matrix, trace_matrix,
                 align_matrix.size() - 1,
                 align_matrix.back().size() - 1,
                 true  /* is_global */);
  if (align_matrix.size() <= kMaxPrintedWidth &&
      align_matrix.back().size() <= kMaxPrintedWidth) {
    PrintMatrix(align_matrix);
  }
  return align_matrix.back().back();
}

//-----------------------------------------------------------------------------

void Alignment::GetNextScore(const vector<vector<int>>& align_matrix,
                             const int row,
                             const int col,
                             const bool allow_reset,
                             int *const score,
                             char *const trace) const {

  const string& sequence_str_1 = sequence_1_.sequence_str;
  const string& sequence_str_2 = sequence_2_.sequence_str;

  int max_score;
  if (row > 0 && col > 0) {
    // Character indexes in string that we are attempting to align. We must
    // subtract 1 from row and col because the matrix contains one extra row
    // and column for the empty string comparison.
    const int blosum_index_1 =
      static_cast<int>(CharToAmino(sequence_str_1[row - 1]));
    const int blosum_index_2 =
      static_cast<int>(CharToAmino(sequence_str_2[col - 1]));
    max_score = align_matrix[row - 1][col - 1] +
                  kBlosum62[blosum_index_1][blosum_index_2];
  } else {
    // If we are doing global alignment we may be starting with gaps. In this
    // case we  will set this initial score to negative infinity so as to take
    // on the score traversing the edge of the matrix.
    max_score = numeric_limits<int>::min();
  }
  // We will fill the traceback as follows:
  // x: Denotes that we are traversing by matching the two characters.
  // l: Denotes that we are traversing by adding a gap to the second sequence.
  // u: Denotes that we are traversing by adding a gap to the first sequence.
  // r: Denotes that we are resetting the string as the score is below zero.
  *trace = 'x';

  // Get best score depending on operation performed using previous best
  // scores.
  if (row > 0 && align_matrix[row - 1][col] - kGapScore > max_score) {
    // We have moved from the 'left' tile by adding a gap to the second
    // sequence.
    *trace = 'l';
    max_score = align_matrix[row - 1][col] - kGapScore;
  }

  if (col > 0 && align_matrix[row][col - 1] - kGapScore > max_score) {
    // We have moved from the 'up' tile by adding a gap to the first
    // sequence.
    *trace = 'u';
    max_score = align_matrix[row][col - 1] - kGapScore;
  }

  if (allow_reset && 0 > max_score) {
    *trace = 'r';
    max_score = 0;
  }
  *score = max_score;
}

//-----------------------------------------------------------------------------

void Alignment::PrintTraceback(const vector<vector<int>>& align_matrix,
                               const vector<vector<char>>& trace_matrix,
                               const int start_row,
                               const int start_col,
                               const bool is_global) const {

  vector<pair<char, int>> aligned_sequence_1;
  vector<pair<char, int>> aligned_sequence_2;
  vector<char> aligned_relation;
  int ii = start_row;
  int kk = start_col;

  while (true) {
    assert(ii >= 0);
    assert(kk >= 0);
    // Traceback defined as follows.
    // x: Denotes that we are traversing by matching the two characters.
    // l: Denotes that we are traversing by adding a gap to the second
    //    sequence.
    // u: Denotes that we are traversing by adding a gap to the first
    //    sequence.
    // r: Denotes that we are resetting the string as the score is below zero.
    if (trace_matrix[ii][kk] == 'x') {
      // Trace back to top left corner and add both letters.
      aligned_sequence_1.push_back(
        pair<char, int>(sequence_1_.sequence_str[ii - 1], ii - 1));
      aligned_sequence_2.push_back(
        pair<char, int>(sequence_2_.sequence_str[kk - 1], kk - 1));
      char char1 = aligned_sequence_1.back().first;
      char char2 = aligned_sequence_2.back().first;
      const int blosum_index_1 = static_cast<int>(CharToAmino(char1));
      const int blosum_index_2 = static_cast<int>(CharToAmino(char2));
      if (char1 == char2) {
        aligned_relation.push_back(char1);
      } else if (kBlosum62[blosum_index_1][blosum_index_2] > 0){
        aligned_relation.push_back('+');
      } else {
        aligned_relation.push_back(' ');
      }
      --ii;
      --kk;
      continue;
    }

    if (trace_matrix[ii][kk] == 'l') {
      // Trace back to left tile and add gap to first sequence.
      aligned_sequence_1.push_back(
        pair<char, int>(sequence_1_.sequence_str[ii - 1], ii - 1));
      aligned_sequence_2.push_back(pair<char, int>('-', kk));
      aligned_relation.push_back(' ');
      --ii;
      continue;
    }

    if (trace_matrix[ii][kk] == 'u') {
      // Trace back to above tile and add gap to second sequence.
      aligned_sequence_1.push_back(pair<char, int>('-', ii));
      aligned_sequence_2.push_back(pair<char, int>(
        sequence_2_.sequence_str[kk - 1], kk - 1));
      aligned_relation.push_back(' ');
      --kk;
      continue;
    }

    if (trace_matrix[ii][kk] == 'r') {
      // We have reached a reset meaning the current string was started after
      // this tile. This is equivalent to an initial zero also the trace matrix
      // is prefilled with resets. Global alignment can fluctuate below and
      // above zero so this symbol is useful in traceback. If this is a global
      // alignment then we must make sure all chars have been read otherwise
      // we can return once we have run out of chars in one sequence.
      if (!is_global) {
        break;
      }

      while (ii > 0) {
        aligned_sequence_1.push_back(
          pair<char, int>(sequence_1_.sequence_str[ii - 1], ii - 1));
        aligned_sequence_2.push_back(pair<char, int>('-', kk));
        aligned_relation.push_back(' ');
        --ii;
      }
      assert(ii == 0);

      while (kk > 0) {
        aligned_sequence_1.push_back(pair<char, int>('-', ii));
        aligned_sequence_2.push_back(pair<char, int>(
          sequence_2_.sequence_str[kk - 1], kk - 1));
        aligned_relation.push_back(' ');
        --kk;
      }
      assert(kk == 0);

      // We have reached the top left tile No more letters to process.
      break;
    }
    cout << "Unrecognized tracback: " << trace_matrix[ii][kk] << endl;
    assert(false);
  }

  assert(aligned_sequence_1.size() == aligned_sequence_2.size());
  assert(aligned_relation.size() == aligned_sequence_2.size());
  while (!aligned_sequence_1.empty() &&
         !aligned_sequence_2.empty()) {

    PrintLinePrefix(sequence_1_.identifier, aligned_sequence_1.back().second);
    for (int ii = 0; ii < kColumnWidth && !aligned_sequence_1.empty(); ++ii) {
      cout << aligned_sequence_1.back().first;
      aligned_sequence_1.pop_back();
    }
    cout << endl;

    PrintLinePrefix("", -1);
    for (int ii = 0; ii < kColumnWidth && !aligned_relation.empty(); ++ii) {
      cout << aligned_relation.back();
      aligned_relation.pop_back();
    }
    cout << endl;

    PrintLinePrefix(sequence_2_.identifier, aligned_sequence_2.back().second);
    for (int ii = 0; ii < kColumnWidth && !aligned_sequence_2.empty(); ++ii) {
      cout << aligned_sequence_2.back().first;
      aligned_sequence_2.pop_back();
    }
    cout << endl;
    cout << endl;
  }
}

//-----------------------------------------------------------------------------

void Alignment::PrintLinePrefix(const string& identifier,
                                const int index) const {

    // We will pad the prefix so as to align the sequence starts.
    static const int kLinePrefixLength = 15;
    const string prefix =
      index < 0 ? "" : identifier + ": " + to_string(index) + " ";
    cout << prefix;
    for (int ii = prefix.size(); ii < kLinePrefixLength; ++ii) {
      cout << " ";
    }
}

//-----------------------------------------------------------------------------

vector<vector<int>> Alignment::CreateAlignmentMatrix() const {
  // Alignment matrix that will be used to fill in our running scores. Note
  // that an extra row and column are added denoting the alignment against an
  // empty string. We will initialize the whole matrix to zeroes.
  vector<vector<int>> align_matrix;
  for (int ii = 0; ii < sequence_1_.sequence_str.size() + 1; ++ ii) {
    align_matrix.push_back(vector<int>(
      sequence_2_.sequence_str.size() + 1, 0));
  }
  return align_matrix;
}

//-----------------------------------------------------------------------------

vector<vector<char>> Alignment::CreateTracebackMatrix() const {
  // Matrix containing chars that will act as crumbs to traceback from any
  // point.
  vector<vector<char>> trace_matrix;
  for (int ii = 0; ii < sequence_1_.sequence_str.size() + 1; ++ ii) {
    trace_matrix.push_back(vector<char>(
      sequence_2_.sequence_str.size() + 1, 'r'));
  }
  return trace_matrix;
}

//-----------------------------------------------------------------------------

void Alignment::PrintSequences() const {
  cout << "### Sequence 1 ###" << endl;
  PrintWithColumns(sequence_1_.sequence_str);
  cout << "### Sequence 2 ###" << endl;
  PrintWithColumns(sequence_2_.sequence_str);
}

//-----------------------------------------------------------------------------

void Alignment::PrintWithColumns(const string& str) const {
  for (int ii = 0; ii < str.size(); ++ii) {
    if (ii % 60 == 0 && ii != 0) {
      cout << endl;
    }
    cout << str[ii];
  }
  cout << endl;
}

//-----------------------------------------------------------------------------

void Alignment::PrintMatrix(const vector<vector<int>>& matrix) const {
  for (int ii = 0; ii < matrix.size(); ++ii) {
    for (int kk = 0; kk < matrix[ii].size(); ++kk) {
      if (matrix[ii][kk] >= 0) {
        cout << " ";
      }
      if ((matrix[ii][kk] / 10) == 0) {
        cout << " ";
      }
      cout << matrix[ii][kk] << " ";
    }
    cout << endl;
  }
}

//-----------------------------------------------------------------------------
