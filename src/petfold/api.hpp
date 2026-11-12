#pragma once
/*
 * api.hpp - high level PETfold / PETcofold API used by the command-line
 * tools and the Python binding.
 */
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "types.hpp"

namespace petfold {

inline constexpr const char *kVersion = "3.0.0";  // C++ edition (PETfold 2.2 / PETcofold 3.4 models)

/* ------------------------------------------------------------------ */
/* PETfold                                                             */

struct PetfoldInput {
  Alignment alignment;
  std::string tree;       // Newick tree ("" = neighbour joining + ML branch lengths)
  std::string structure;  // structure to score ("" = predict)
  int suboptimal = 0;     // number of sampled structures
  std::FILE *ppfold_input = nullptr;  // read Pfold reliabilities from here (PPfold)
};

struct PetfoldResult {
  int alignment_length = 0;          // original number of columns
  std::vector<int> gap_columns;      // removed columns (0-based)
  Alignment alignment_gf;            // gap-column-free alignment
  std::string tree;                  // tree used (ML tree if estimated)
  bool has_evo = false;
  std::string pfold_structure_gf;    // Pfold MEA structure (gap-free)
  std::string constraints_gf;        // evolutionary constraints '(' ')' 'x'
  Matrix paired_tree, paired_seq, paired_pet;          // gap-free, see core.hpp
  std::vector<double> single_tree, single_seq, single_pet;
  std::string structure_gf;          // PETfold structure (gap-free)
  std::string structure;             // PETfold structure with gap columns
  double score = 0;
  double ensemble_diversity = 0;
  std::vector<std::pair<std::string, double>> suboptimals;  // with gap columns
};

PetfoldResult run_petfold(Context &ctx, const PetfoldInput &in);

/* ------------------------------------------------------------------ */
/* PETcofold                                                           */

struct PetcofoldInput {
  Alignment alignment1, alignment2;
  std::string tree;
  std::string structure;   // joint structure "...&..." to score
  std::string structure1;  // known structure of alignment 1
  std::string structure2;  // known structure of alignment 2
  float delta = 0.9f;      // --setdelta (petcon)
  float gamma = 0.1f;      // --setgamma (partprob)
  bool extstem = false;
  bool legacy_cutpoint = false;  // reproduce the PETcofold 3.4 cut point handling
};

struct PetcofoldResult {
  std::vector<std::string> common_ids;
  std::vector<int> gap_columns1, gap_columns2;
  int length1 = 0, length2 = 0;  // original alignment lengths
  float delta1 = 0, delta2 = 0;
  double partial_prob1[2] = {0, 0}, partial_prob2[2] = {0, 0};  // {thermodynamic, evolutionary}
  std::string partial_structure1, partial_structure2;          // with gap columns
  std::string petfold_structure1, petfold_structure2;          // step-1 structures, with gap columns
  std::string structure1, structure2;                          // joint structure parts, with gap columns
  std::string intermolecular1, intermolecular2;
  double score = 0;
  bool has_evo = false;
  Matrix duplex_paired_pet;            // gap-free duplex reliabilities (len1+len2)
  std::vector<double> duplex_single_pet;
  std::vector<std::string> fasta_ids, fasta_seq1, fasta_seq2;  // for --war output
};

PetcofoldResult run_petcofold(Context &ctx, const PetcofoldInput &in, const std::string &ppfile = "");

/* ------------------------------------------------------------------ */
/* command line front ends (identical options and output to 2.2/3.4)   */

int petfold_main(int argc, char **argv);
int petcofold_main(int argc, char **argv);

}  // namespace petfold
