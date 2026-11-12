#pragma once
/*
 * core.hpp - building blocks of the PETfold model (formerly petfoldlibs.c,
 * thermodynamic.c and evolutionary.c).
 */
#include <string>
#include <vector>

#include "types.hpp"

namespace petfold {

/* ---------- alignment handling (alignment.cpp) ---------- */

Alignment read_fasta(const std::string &path);
Alignment parse_fasta(const std::string &text, const std::string &source = "<string>");
/* Returns the first line of `arg` if it names a readable file, else `arg`. */
std::string read_file_or_string(const std::string &arg);
/* Removes columns whose gap fraction is >= gap; gap_col receives the removed
   column indices in ascending order. */
Alignment delete_gap_columns(const Context &ctx, const Alignment &align, std::vector<int> &gap_col, float gap);
std::string include_gap_columns(const std::string &db, const std::vector<int> &gap_col, int full_len);
/* ungapped sequences plus map to alignment columns */
struct SeqList {
  std::vector<std::string> seq;
  std::vector<std::vector<int>> origid;
};
SeqList get_align_without_gaps(const Alignment &align);

/* ---------- structures (structure.cpp) ---------- */

std::string get_dot_bracket(const std::vector<int> &struc);
PartStruc get_partial_struc(const Matrix &paired, const std::vector<double> &single,
                            const std::vector<int> &struc, int len, float evocon_bp, float evocon_ss);
std::string get_constraint_string(const PartStruc *evocon, const std::vector<int> &origid, int len_gf);
/* Parses a dot-bracket structure given in alignment coordinates, checks it and
   removes gap columns.  Accepts ()<>[] and .- */
PartStruc mod_setstruct(std::string setstruc, const std::vector<int> &gap_col, int len, int gapnr);
std::vector<int> partstruc2list(const PartStruc *evocon, int len);
std::vector<int> partstruc2str(const PartStruc *evocon, int len);
std::string partstruc2string(const PartStruc *evocon, int len);
PartStruc str2partstruc(const std::vector<int> &struc, int len);

/* ---------- PET model and MEA folding (mea.cpp) ---------- */

std::vector<double> get_pet_single(const std::vector<double> &single_tree, const std::vector<double> &single_seq, float beta);
Matrix get_pet_paired(const Matrix &paired_tree, const Matrix &paired_seq, float beta);

/* MEA dynamic programming table over the Knudsen-Hein grammar S/F/L */
class MeaTable {
 public:
  MeaTable(const Matrix &paired_pet, const std::vector<double> &single_pet, const PartStruc *constraint, float alpha);
  double score() const { return rel(0, 0, n_ - 1); }
  std::string structure() const;  // MEA structure by traceback
  int size() const { return n_; }
 private:
  std::size_t idx(int t, int i, int j) const { return (static_cast<std::size_t>(t) * n_ + i) * n_ + j; }
  double &rel(int t, int i, int j) { return rel_[idx(t, i, j)]; }
  double rel(int t, int i, int j) const { return rel_[idx(t, i, j)]; }
  int n_;
  std::vector<double> rel_;
  std::vector<int> code_;  // rule used: low 2 bits = rule type, high bits = split point
};

double ensemble_diversity(const Matrix &paired, int len);

/* structure sampling (heuristic sampler of PETfold 2.2, with fixed RNG handling) */
class SuboptSampler {
 public:
  SuboptSampler(const Matrix &paired_pet, const std::vector<double> &single_pet, float alpha);
  double sample(std::mt19937_64 &rng, std::string &struc) const;
 private:
  double P(int i, int j, int t) const { return P_[(static_cast<std::size_t>(i) * n_ + j) * 3 + t]; }
  const Matrix &paired_;
  const std::vector<double> &single_;
  float alpha_;
  int n_;
  std::vector<double> P_;
};

/* expands gap-free reliabilities to alignment coordinates (zeros in gap
   columns, symmetric pair matrix) exactly as written to the pp-file */
void expand_reliabilities(const Matrix &paired, const std::vector<double> &single, const std::vector<int> &gap_col,
                          int full_len, Matrix &paired_out, std::vector<double> &single_out);
void create_ppfile(const Matrix &paired, const std::vector<double> &single, const std::vector<int> &gap_col,
                   int len, int gapnr, const std::string &ppfile);

/* ---------- thermodynamic model (thermodynamic.cpp) ---------- */

Matrix get_prob_paired_seq(Context &ctx, const Alignment &align_gf, const PartStruc *evocon);
std::vector<double> get_prob_unpaired_seq(const Context &ctx, const Alignment &align_gf, const Matrix &paired_seq);
double get_fold_thermo_partial_prob(double econstraint, double eensemble);
double get_constraint_energy(Context &ctx, const Alignment &align, const PartStruc *evocon);
double get_ensemble_energy(Context &ctx, const Alignment &align);
/* co-folding of concatenated sequences; len1 = gap-free length of alignment 1.
   With legacy_cutpoint the cut is placed as in PETcofold 3.4 (after position
   len1 of every ungapped sequence), otherwise at the true per-sequence
   boundary between the two molecules. */
void get_cofold_prob_paired_seq(Context &ctx, const Alignment &align, int len1, Matrix &paired_seq,
                                const PartStruc *pstruc, bool legacy_cutpoint);

/* ---------- evolutionary model (evolutionary.cpp) ---------- */

std::string grammar_file(const Context &ctx, const std::string &name);
/* last directory searched for the grammar files (set by the Python package) */
void set_default_data_dir(const std::string &dir);
/* Pfold: tree (NJ + ML branch lengths unless given), posterior pair and
   unpaired probabilities and the Pfold MEA structure.  `tree` is input
   (Newick, may be empty) and output (the tree actually used, as in 2.2). */
std::vector<int> get_prob_tree(Context &ctx, const Alignment &align_gf, Matrix &paired_tree,
                               std::vector<double> &single_tree, std::string &tree,
                               const PartStruc *constraint = nullptr);
std::vector<int> get_ppfold(int len, Matrix &paired_tree, std::vector<double> &single_tree, std::FILE *in);
double get_fold_evol_partial_prob(double evolconstraint, double evolensemble);
double get_constraint_evolprob(Context &ctx, const Alignment &align, const std::string &tree, const PartStruc &evocon);
double get_ensemble_evolprob(Context &ctx, const Alignment &align, const std::string &tree);

}  // namespace petfold
