/*
 * Driver for sequence alignment.
 * Usage:
 * ./a.out <global|local> <file1_path> < file2_path>
 *
 * Files provided should contain valid amino acids and be in FASTA format.
 *
 * Author: Dan Lin, danielin@uw.edu
 */

#include <cstdlib>
#include <iostream>
#include <time.h>
#include "alignment.h"

using namespace std;

// Number of permutations that will be performed in order to calculate pvalue.
static const int kNumPermutations = 999;

// Print help for command line args.
void PrintHelp(const char *const help_str);

string Permute(string permutation);

int main(int argc, char *argv[]) {
  if (argc < 4) {
    PrintHelp("Invalid number of args");
    return -1;
  }

  string file_1(argv[2]);
  string file_2(argv[3]);
  Alignment align(file_1, file_2);

  int match_score;
  if (strncmp(argv[1], "local", 5) == 0) {
    match_score = align.GetLocalAlignment();
  } else if (strncmp(argv[1], "global", 6) == 0) {
    match_score = align.GetGlobalAlignment();
  } else {
    PrintHelp(argv[1]);
  }

  int n = kNumPermutations;
  if (argc == 5) {
    n = stoi(string(argv[4]));
  }

  int k = 1;
  srand(time(NULL));
  for (int ii = 0; ii < n; ++ii) {
    align.set_sequence_str2(Permute(align.sequence_str2()));
    int score;
    if (strncmp(argv[1], "local", 5) == 0) {
      score = align.GetLocalAlignment(false /* print */);
    } else if (strncmp(argv[1], "global", 6) == 0) {
      score = align.GetGlobalAlignment(false /* print */);
    }

    if (score >= match_score) {
      ++k;
    }
  }
  double pvalue = (1.0 * k) / (n + 1);
  cout << "P-Value: " << pvalue << " after " << n
       << " permutations" << endl;
}

void PrintHelp(const char *const invalid_arg) {
  cout << "Invalid argument provided: " << invalid_arg << endl
       << "Expected use: " << endl
       << "./a.out <global|local> <file1> <file2> [print_pvalue]" << endl;
}

string Permute(string permutation) {
  for (int ii = 0; ii < permutation.size(); ++ii) {
    int rand_index = rand() % (ii + 1);
    swap(permutation[ii], permutation[rand_index]);
  }
  return permutation;
}
