#pragma once
/*
 * types.hpp - core data types of the PETfold / PETcofold C++ library.
 *
 * Numerical conventions are identical to PETfold 2.2 / PETcofold 3.4:
 * parameters that were `float` in the original C code are still `float`
 * so that results are bit-compatible, all matrices are indexed in gap-free
 * alignment coordinates, and INFINITE_D marks undefined entries.
 */
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace petfold {

inline constexpr double INFINITE_D = 1000000000.;  // value of undefined double fields
inline constexpr int INFINITE_I = 1000000000;      // value of undefined integer fields

/* Error raised by the library.  `to_stdout` preserves where the original
   command-line tools printed the corresponding message. */
class Error : public std::runtime_error {
 public:
  explicit Error(const std::string &msg, bool to_stdout = false)
      : std::runtime_error(msg), to_stdout_(to_stdout) {}
  bool to_stdout() const { return to_stdout_; }
 private:
  bool to_stdout_;
};

/* printf-style formatting into a std::string (keeps the exact output format
   of the original tools). */
std::string strf(const char *fmt, ...)
#if defined(__GNUC__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;

/* Dense n x n matrix of doubles, row-major. */
class Matrix {
 public:
  Matrix() = default;
  explicit Matrix(int n, double init = 0.0) : n_(n), v_(static_cast<std::size_t>(n) * n, init) {}
  double &operator()(int i, int j) { return v_[static_cast<std::size_t>(i) * n_ + j]; }
  double operator()(int i, int j) const { return v_[static_cast<std::size_t>(i) * n_ + j]; }
  int size() const { return n_; }
  const std::vector<double> &data() const { return v_; }
  std::vector<double> &data() { return v_; }
 private:
  int n_ = 0;
  std::vector<double> v_;
};

/* Multiple alignment.  `seqs` are upper case, T->U, gaps normalised to '-'.
   `raw` keeps the sequences as given (for FASTA-style output). */
struct Alignment {
  std::vector<std::string> ids;
  std::vector<std::string> seqs;
  std::vector<std::string> raw;
  int nr() const { return static_cast<int>(seqs.size()); }
  int len() const { return seqs.empty() ? 0 : static_cast<int>(seqs[0].size()); }
};

/* Partial (or complete) structure: constrained base pairs and unpaired bases. */
struct PartStruc {
  std::vector<int> bp_left, bp_right, ss;
  int bp_nr() const { return static_cast<int>(bp_left.size()); }
  int ss_nr() const { return static_cast<int>(ss.size()); }
};

/* Model parameters shared by PETfold and PETcofold. */
struct Options {
  float gap = 0.25f;        // max. fraction of gaps in a column (column removed if >=)
  float evocon_bp = 0.9f;   // reliability threshold for conserved base pairs
  float evocon_ss = 1.0f;   // reliability threshold for conserved unpaired bases
  float alpha = 0.2f;       // weight of unpaired reliabilities
  float beta = 1.0f;        // weight of the thermodynamic model
  bool use_evo = true;      // evolutionary model (Pfold)
  bool use_thermo = true;   // thermodynamic model (ViennaRNA)
  bool verbose = false;     // long output
  std::string grammar_dir;  // directory holding scfg.rate and article.grm ("" = search)
};

/* Message sink.  The command-line tools print to stdout/stderr, the Python
   binding collects the messages. */
struct Log {
  std::function<void(const std::string &)> out;  // normal/verbose output (stdout)
  std::function<void(const std::string &)> err;  // diagnostics (stderr)
  void print(const std::string &s) const { if (out) out(s); }
  void warn(const std::string &s) const { if (err) err(s); }
};

/* Per-run state. */
struct Context {
  Options opt;
  Log log;
  /* Replicates the global `pf_scale` of the deprecated ViennaRNA interface:
     it is set when a sequence longer than 200 nt is folded and then kept
     for all following partition function computations. */
  double pf_scale = -1.0;
  std::mt19937_64 rng{0};
};

}  // namespace petfold
