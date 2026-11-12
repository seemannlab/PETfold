/*
 * mea.cpp - PET reliabilities, MEA folding, sampling and pp-files
 * (from petfoldlibs.c).
 */
#include <cmath>
#include <cstdio>

#include "core.hpp"

namespace petfold {

namespace {
enum : int { S = 0, F = 1, L = 2 };
enum : int { R_NONE = 0, R_UNIT = 1, R_LS = 2, R_DFD = 3 };  // R_UNIT: S->L or L->s

inline double add(double x, double y) { return (x != INFINITE_D && y != INFINITE_D) ? x + y : INFINITE_D; }
/* true if x should replace y (same semantics as max() in 2.2) */
inline bool better(double x, double y) {
  if (x != INFINITE_D && y != INFINITE_D) return x >= y;
  return y == INFINITE_D;
}
}  // namespace

std::vector<double> get_pet_single(const std::vector<double> &single_tree, const std::vector<double> &single_seq,
                                   float beta) {
  std::vector<double> s(single_tree.size());
  for (std::size_t i = 0; i < s.size(); i++) s[i] = (single_tree[i] + beta * single_seq[i]) / 2;
  return s;
}

Matrix get_pet_paired(const Matrix &paired_tree, const Matrix &paired_seq, float beta) {
  const int len = paired_tree.size();
  Matrix p(len, INFINITE_D);
  for (int j = 1; j < len; j++)
    for (int i = 0; i < len - j; i++) p(i, i + j) = (paired_tree(i, i + j) + beta * paired_seq(i, i + j)) / 2;
  return p;
}

MeaTable::MeaTable(const Matrix &paired_pet, const std::vector<double> &single_pet, const PartStruc *constraint,
                   float alpha)
    : n_(paired_pet.size()) {
  const int len = n_;
  if (len <= 0) throw Error("The consensus RNA structure could not be predicted!");
  rel_.assign(static_cast<std::size_t>(3) * len * len, INFINITE_D);
  code_.assign(static_cast<std::size_t>(3) * len * len, R_NONE);
  const std::vector<int> consstr = partstruc2list(constraint, len);
  auto set = [&](int t, int i, int j, double v, int rule, int k = 0) {
    rel_[idx(t, i, j)] = v;
    code_[idx(t, i, j)] = rule | (k << 2);
  };

  for (int i = 0; i < len; i++) {
    if (consstr[i] != INFINITE_I && consstr[i] != -1) continue;
    set(L, i, 0, 2 * alpha * single_pet[i] / len, R_UNIT);  // L->s
    set(S, i, 0, rel(L, i, 0), R_UNIT);                      // S->L
  }

  for (int j = 0; j < len; j++) {
    for (int i = 0; i < len - j; i++) {
      for (int k = 1; k < j + 1; k++) {  // S->LS and F->LS
        const double l = rel(L, i, k - 1), s = rel(S, i + k, j - k);
        if (l == INFINITE_D || s == INFINITE_D) continue;
        const double r = add(l, s);
        if (better(r, rel(S, i, j))) set(S, i, j, r, R_LS, k);
        if (better(r, rel(F, i, j))) set(F, i, j, r, R_LS, k);
      }
      if (j > 3) {
        const double inner = rel(F, i + 1, j - 2);
        const bool allowed = (consstr[i] == INFINITE_I && consstr[i + j] == INFINITE_I) ||
                             (inner != INFINITE_D && consstr[i] == consstr[i + j] && consstr[i] != -1);
        if (allowed) {
          const double r = add(inner, 2 * paired_pet(i, i + j) / len);
          if (better(r, rel(F, i, j))) set(F, i, j, r, R_DFD);  // F->dFd
          if (better(r, rel(L, i, j))) set(L, i, j, r, R_DFD);  // L->dFd
        }
      }
      if (rel(L, i, j) != INFINITE_D) {  // S->L
        const double r = rel(L, i, j);
        if (better(r, rel(S, i, j))) set(S, i, j, r, R_UNIT);
      }
    }
  }
  if (rel(S, 0, len - 1) == INFINITE_D) throw Error("The consensus RNA structure could not be predicted!");
}

std::string MeaTable::structure() const {
  std::string db(n_, '.');
  struct Item { int t, i, j; };
  std::vector<Item> stack{{S, 0, n_ - 1}};
  while (!stack.empty()) {
    const Item it = stack.back();
    stack.pop_back();
    const int c = code_[idx(it.t, it.i, it.j)];
    const int rule = c & 3, k = c >> 2;
    switch (rule) {
      case R_LS:
        stack.push_back({S, it.i + k, it.j - k});
        stack.push_back({L, it.i, k - 1});
        break;
      case R_UNIT:
        if (it.t == S) stack.push_back({L, it.i, it.j});
        else db[it.i] = '.';  // L->s
        break;
      case R_DFD:
        db[it.i] = '(';
        db[it.i + it.j] = ')';
        stack.push_back({F, it.i + 1, it.j - 2});
        break;
      default:
        throw Error("internal error: MEA traceback reached an undefined cell");
    }
  }
  return db;
}

double ensemble_diversity(const Matrix &paired, int len) {
  int c = 0;
  double ensdiv = .0;
  for (int i = 0; i < len; i++)
    for (int j = i + 4; j < len; j++) {
      ensdiv += paired(i, j) * (1 - paired(i, j));
      c++;
    }
  return ensdiv / c;
}

/* ---------------------------------------------------------------------- */
/* structure sampling                                                      */

SuboptSampler::SuboptSampler(const Matrix &paired_pet, const std::vector<double> &single_pet, float alpha)
    : paired_(paired_pet), single_(single_pet), alpha_(alpha), n_(paired_pet.size()) {
  const int len = n_;
  P_.assign(static_cast<std::size_t>(len) * len * 3, INFINITE_D);
  auto p = [&](int i, int j, int t) -> double & { return P_[(static_cast<std::size_t>(i) * len + j) * 3 + t]; };
  for (int i = 0; i < len; i++) {
    p(i, 0, L) = 2 * alpha * single_pet[i] / len;
    p(i, 0, S) = p(i, 0, L);
  }
  for (int j = 0; j < len; j++)
    for (int i = 0; i < len - j; i++) {
      for (int k = 1; k < j + 1; k++) {
        if (p(i, k - 1, L) == INFINITE_D || p(i + k, j - k, S) == INFINITE_D) continue;
        if (p(i, j, S) == INFINITE_D) p(i, j, S) = 0.;
        p(i, j, S) = add(p(i, j, S), add(p(i, k - 1, L), p(i + k, j - k, S)));
        if (p(i, j, F) == INFINITE_D) p(i, j, F) = 0.;
        p(i, j, F) = add(p(i, j, F), add(p(i, k - 1, L), p(i + k, j - k, S)));
      }
      if (j > 3) {
        if (p(i, j, F) == INFINITE_D) p(i, j, F) = 0;
        p(i, j, F) = add(p(i, j, F), add(p(i + 1, j - 2, F), 2 * paired_pet(i, i + j) / len));
        if (p(i, j, L) == INFINITE_D) p(i, j, L) = 0;
        p(i, j, L) = add(p(i, j, L), add(p(i + 1, j - 2, F), 2 * paired_pet(i, i + j) / len));
      }
      if (p(i, j, L) != INFINITE_D) {
        if (p(i, j, S) == INFINITE_D) p(i, j, S) = 0;
        p(i, j, S) = add(p(i, j, S), p(i, j, L));
      }
    }
  if (p(0, len - 1, S) == INFINITE_D) throw Error("The consensus RNA structure could not be predicted!");
}

namespace {
/* 2.2: round(x) macro = (unsigned long)(x+1) for x >= 0 */
inline double weight(double x) { return x >= 0 ? std::floor(x + 1.0) : std::ceil(x - 1.0); }

/* FIX: 2.2 used int sums (overflow) and divided by zero for empty lists */
void change_list_sum_to_max(std::vector<double> &w, int start, int end) {
  double sum = 0, maxi = 0;
  for (int i = start - 1; i < end; i++) {
    sum += w[i];
    if (w[i] >= maxi) maxi = w[i];
  }
  if (sum <= 0) return;
  for (int i = start - 1; i < end; i++) w[i] = std::floor(maxi * w[i] / sum);
}

/* FIX: 2.2 used rand()%sum (division by zero for sum 0, modulo bias, and an
   off-by-one that could select zero-weight entries). */
int statistical_sampling(const std::vector<double> &w, int n, std::mt19937_64 &rng) {
  double sum = 0;
  for (int i = 0; i < n; i++) sum += w[i];
  if (sum <= 0) return static_cast<int>(std::uniform_int_distribution<int>(0, n - 1)(rng));
  const double sample = std::uniform_real_distribution<double>(0.0, sum)(rng);
  double cum = 0;
  for (int i = 0; i < n; i++) {
    cum += w[i];
    if (sample < cum) return i;
  }
  return n - 1;
}
}  // namespace

double SuboptSampler::sample(std::mt19937_64 &rng, std::string &struc) const {
  const int len = n_;
  /* SCFG rule probabilities from Knudsen & Hein (2003) */
  const double S_LS = 0.868534, S_L = 0.131466, F_dFd = 0.787640, F_LS = 0.212360, L_dFd = 0.105397;
  const double CONST = 1000000;
  std::vector<int> rid{S}, depth{0}, istack{0}, jstack{len - 1};
  std::vector<double> sam_prob(len + 1);
  std::vector<int> klist(len + 1);
  double score = 0;
  struc.assign(len, '\0');
  std::string form = "S";  // sentential form, rewritten left to right as in 2.2
  auto rewrite = [&](char nt, const char *rhs) {
    const std::size_t off = form.find(nt);
    form.replace(off, 1, rhs);
  };
  auto pop_depth = [&]() {
    for (int k = 0; k < depth.back(); k++) {
      istack.pop_back();
      rid.pop_back();
    }
    depth.pop_back();
  };
  while (true) {
    const int i = istack.back(), j = jstack.back();
    const int cur = rid.back();
    if (cur == L && j == 0) {  // L->s
      rewrite('L', ".");
      score = add(score, 2 * alpha_ * single_[i] / len);
      pop_depth();
    } else if (cur == S && j == 0) {  // S->L->s
      rewrite('S', ".");
      score = add(score, 2 * alpha_ * single_[i] / len);
      pop_depth();
    } else if (cur == L) {  // L->dFd
      if (j <= 3) throw Error("Unpaired bases were not reached!");
      rewrite('L', "(F)");
      score = add(score, 2 * paired_(i, i + j) / len);
      rid.push_back(F);
      istack.push_back(i + 1);
      jstack.push_back(j - 2);
      depth.back()++;
    } else {  // S or F: sample between X->LS (all splits) and X->dFd
      const bool isS = (cur == S);
      int sam_len = 0;
      for (int k = 1; k < j + 1; k++) {
        if (P(i, k - 1, L) != INFINITE_D && P(i + k, j - k, S) != INFINITE_D) {
          klist[sam_len] = k;
          sam_prob[sam_len++] = weight(CONST * (isS ? S_LS : F_LS) *
                                       add(P(i, k - 1, L) / std::pow(2, k), P(i + k, j - k, S) / std::pow(2, j - k + 1)));
        }
      }
      change_list_sum_to_max(sam_prob, 2, sam_len);
      if (j > 3 && P(i + 1, j - 2, F) != INFINITE_D)
        sam_prob[sam_len++] = weight(CONST * (isS ? S_L * L_dFd : F_dFd) *
                                     add(P(i + 1, j - 2, F) / std::pow(2, j - 1), 2 * paired_(i, i + j) / len));
      else
        sam_prob[sam_len++] = 0;
      const int smp = statistical_sampling(sam_prob, sam_len, rng);
      const char nt = isS ? 'S' : 'F';
      if (smp < sam_len - 1) {  // X->LS
        rewrite(nt, "LS");
        rid.push_back(S);
        istack.push_back(i + klist[smp]);
        jstack.push_back(j - klist[smp]);
        depth.back()++;
        rid.push_back(L);
        istack.push_back(i);
        jstack.push_back(klist[smp] - 1);
        depth.push_back(1);
      } else {  // X->dFd
        rewrite(nt, "(F)");
        score = add(score, 2 * paired_(i, i + j) / len);
        rid.push_back(F);
        istack.push_back(i + 1);
        jstack.push_back(j - 2);
        depth.back()++;
      }
    }
    jstack.resize(istack.size());
    if (rid.size() == 1) break;
  }
  struc = form;
  return score;
}

/* ---------------------------------------------------------------------- */

void expand_reliabilities(const Matrix &paired, const std::vector<double> &single, const std::vector<int> &gap_col,
                          int full_len, Matrix &paired_out, std::vector<double> &single_out) {
  std::vector<int> map(full_len, -1);
  std::size_t g = 0;
  int p = 0;
  for (int col = 0; col < full_len; col++) {
    if (g < gap_col.size() && gap_col[g] == col) g++;
    else map[col] = p++;
  }
  paired_out = Matrix(full_len, 0.0);
  single_out.assign(full_len, 0.0);
  for (int a = 0; a < full_len; a++) {
    if (map[a] < 0) continue;
    single_out[a] = single[map[a]];
    for (int b = 0; b < full_len; b++) {
      if (map[b] < 0) continue;
      const int i = map[a], j = map[b];
      double v = paired(i, j);
      if (v == INFINITE_D) v = paired(j, i);
      if (v == INFINITE_D) v = 0.0;
      paired_out(a, b) = v;
    }
  }
}

void create_ppfile(const Matrix &paired, const std::vector<double> &single, const std::vector<int> &gap_col, int len,
                   int gapnr, const std::string &ppfile) {
  std::FILE *fp = std::fopen(ppfile.c_str(), "w");
  if (!fp) throw Error("ppfile: Error in opening file '" + ppfile + "'");
  Matrix pm;
  std::vector<double> sv;
  expand_reliabilities(paired, single, gap_col, len + gapnr, pm, sv);
  const int n = len + gapnr;
  std::fprintf(fp, "%i\n", n);
  for (int a = 0; a < n; a++) {
    for (int b = 0; b < n; b++) std::fprintf(fp, "%10.9e ", pm(a, b));
    std::fprintf(fp, "\n");
  }
  std::fprintf(fp, "\n");
  for (int a = 0; a < n; a++) std::fprintf(fp, "%10.9e ", sv[a]);
  std::fprintf(fp, "\n");
  if (std::fclose(fp) != 0) throw Error("ppfile: Error in closing file '" + ppfile + "'");
}

}  // namespace petfold
