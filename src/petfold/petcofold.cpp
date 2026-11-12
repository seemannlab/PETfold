/*
 * petcofold.cpp - the PETcofold pipeline (from petcofold.c and
 * petcofoldlibs.c).
 */
#include <algorithm>
#include <cmath>

#include "api.hpp"
#include "core.hpp"

namespace petfold {

namespace {

inline double add(double x, double y) { return (x != INFINITE_D && y != INFINITE_D) ? x + y : INFINITE_D; }
double geometric_mean(double a, double b) { return std::sqrt(a * b); }

/* indices of identifiers common to both alignments (first match), in the
   order of alignment 1 */
void intersect_alignments(const Alignment &a1, const Alignment &a2, std::vector<int> &idx1, std::vector<int> &idx2) {
  for (int i = 0; i < a1.nr(); i++)
    for (int j = 0; j < a2.nr(); j++)
      if (a1.ids[i] == a2.ids[j]) {
        idx1.push_back(i);
        idx2.push_back(j);
        break;
      }
}

Alignment get_subalignment(const Alignment &a, const std::vector<int> &idx) {
  Alignment s;
  for (int i : idx) {
    s.ids.push_back(a.ids[i]);
    s.seqs.push_back(a.seqs[i]);
    s.raw.push_back(a.raw[i]);
  }
  return s;
}

/* PETfold for one alignment: reliabilities and MEA structure */
PartStruc run_petfold_step(Context &ctx, const Alignment &align, Matrix &paired_pet, std::vector<double> &single_pet,
                           std::string &tree, const PartStruc *setstruc) {
  const Options &opt = ctx.opt;
  const int len = align.len();
  Matrix paired_tree;
  std::vector<double> single_tree;
  const std::vector<int> pfold_struc = get_prob_tree(ctx, align, paired_tree, single_tree, tree);
  PartStruc evocon;
  const PartStruc *evocon_p = nullptr;
  if (opt.use_evo) {
    evocon = get_partial_struc(paired_tree, single_tree, pfold_struc, len, opt.evocon_bp, opt.evocon_ss);
    evocon_p = &evocon;
  }
  const Matrix paired_seq = get_prob_paired_seq(ctx, align, evocon_p);
  const std::vector<double> single_seq = get_prob_unpaired_seq(ctx, align, paired_seq);
  single_pet = get_pet_single(single_tree, single_seq, opt.beta);
  paired_pet = get_pet_paired(paired_tree, paired_seq, opt.beta);
  MeaTable T(paired_pet, single_pet, setstruc ? setstruc : evocon_p, opt.alpha);
  const std::string db = T.structure();
  return mod_setstruct(db, {}, len, 0);
}

PartStruc extend_partial_stems(const Context &ctx, const PartStruc &ps, const std::vector<int> &inputstruc,
                               const Matrix &paired_pet, int len, float petcon) {
  int newstem = 0, stem_start = -100, stem_l = 0;
  float stempet, ext_sum, stem_sum = 0.;
  std::vector<int> ext(len, INFINITE_I);
  const auto &L = ps.bp_left;
  const auto &R = ps.bp_right;
  for (int i = ps.bp_nr() - 1; i >= 0; i--) {
    ext[L[i]] = R[i];
    ext[R[i]] = L[i];
    if (i == ps.bp_nr() - 1 || R[i] != R[i + 1] + 1) {  // start of a stem
      newstem = 1;
      stem_start = i;
      stem_l = 1;
      stem_sum = paired_pet(L[i], R[i]);
    }
    if (i == 0 || R[i] != R[i - 1] - 1) {  // end of a stem
      ext_sum = stem_sum;                  // outer base pairs
      for (int k = 1; k < len - R[i]; k++) {
        const int outerleft = inputstruc[R[i] + k];
        if (outerleft != -1 && outerleft == L[i] - k) {
          stempet = (ext_sum + paired_pet(L[i] - k, R[i] + k)) / (stem_l + k);
          if (stempet >= petcon) {
            ext[L[i] - k] = R[i] + k;
            ext[R[i] + k] = L[i] - k;
            ext_sum += paired_pet(L[i] - k, R[i] + k);
          } else
            break;
        } else
          break;
      }
      ext_sum = stem_sum;  // inner base pairs
      for (int k = 1; k < R[stem_start]; k++) {
        const int innerleft = inputstruc[R[stem_start] - k];
        if (innerleft == L[stem_start] + k) {
          stempet = (ext_sum + paired_pet(L[stem_start] + k, R[stem_start] - k)) / (stem_l + k);
          if (stempet >= petcon) {
            ext[L[stem_start] + k] = R[stem_start] - k;
            ext[R[stem_start] - k] = L[stem_start] + k;
            ext_sum += paired_pet(L[stem_start] + k, R[stem_start] - k);
          } else
            break;
        } else
          break;
      }
    } else {  // stem continues
      if (!newstem) {
        stem_l++;
        stem_sum += paired_pet(L[i], R[i]);
      }
      newstem = 0;
    }
  }
  PartStruc out = str2partstruc(ext, len);
  if (ctx.opt.verbose) {
    ctx.log.print("partial structure:          " + partstruc2string(&ps, len) + "\n");
    ctx.log.print("extended partial structure: " + partstruc2string(&out, len) + "\n");
  }
  return out;
}

/* partial intra-molecular structure: raise delta until the partial structure
   is probable enough (from get_probable_partial_struc) */
PartStruc get_probable_partial_struc(Context &ctx, const Alignment &align, const Matrix &paired_pet,
                                     const PartStruc &struc_pet, const std::string &tree, float petcon,
                                     float partprob, bool extstem, double *partialprob, float &delta_out) {
  const int len = align.len();
  const std::vector<int> inputstruc = partstruc2str(&struc_pet, len);
  const std::vector<double> temp(len, 0.0);
  float mypetcon = petcon;
  PartStruc evocon;  // FIX: 3.4 returned an uninitialised pointer if the loop never ran (delta >= 1)

  const double thermoensemble = get_ensemble_energy(ctx, align);
  /* FIX (--noevo): without the evolutionary model the evolutionary
     probability is not defined; it is set to 0 (3.4 crashed) */
  const double evolensemble = ctx.opt.use_evo ? get_ensemble_evolprob(ctx, align, tree) : 0.0;

  if (!(mypetcon < 1)) {
    /* delta >= 1: no base pair can be selected; the empty partial structure
       has probability 1 in both models */
    partialprob[0] = 1.0;
    partialprob[1] = ctx.opt.use_evo ? 1.0 : 0.0;
  }
  while (partialprob[0] <= partprob && partialprob[1] <= partprob && mypetcon < 1) {
    petcon = mypetcon;
    evocon = get_partial_struc(paired_pet, temp, inputstruc, len, mypetcon, 1);
    if (extstem) evocon = extend_partial_stems(ctx, evocon, inputstruc, paired_pet, len, mypetcon);

    const double thermoconstraint = get_constraint_energy(ctx, align, &evocon);
    partialprob[0] = get_fold_thermo_partial_prob(thermoconstraint, thermoensemble);
    double evolconstraint = 0.0;
    if (ctx.opt.use_evo) {
      evolconstraint = get_constraint_evolprob(ctx, align, tree, evocon);
      partialprob[1] = get_fold_evol_partial_prob(evolconstraint, evolensemble);
    }
    if (ctx.opt.verbose)
      ctx.log.print(strf("Delta = %3.2f:\n"
                         " average free energy of constrained ensemble=%5.4f kcal/mol and of ensemble=%5.4f kcal/mol "
                         "=> thermodyn partial struct prob = %6.5f\n"
                         " evol score of constrained ensemble=%5.4f and of ensemble=%5.4f => evol partial struct prob = "
                         "%6.5f\n",
                         mypetcon, thermoconstraint, thermoensemble, partialprob[0], evolconstraint, evolensemble,
                         partialprob[1]));
    if (mypetcon >= static_cast<float>(0.9))
      mypetcon += 0.02;
    else
      mypetcon += 0.1;
  }
  ctx.log.print(strf("Delta = %3.2f; thermodynamic partial structure probability = %6.5f; evolutionary partial "
                     "structure probability = %6.5f\n",
                     petcon, partialprob[0], partialprob[1]));
  delta_out = petcon;
  return evocon;
}

/* ss constraint list: bases of both partial structures (+ optional linker) */
PartStruc duplex_blocked(const PartStruc &p1, const PartStruc &p2, int len1, bool linker) {
  PartStruc ps;
  for (int i = 0; i < p1.bp_nr(); i++) ps.ss.push_back(p1.bp_left[i]);
  for (int i = 0; i < p1.bp_nr(); i++) ps.ss.push_back(p1.bp_right[i]);
  const int off = linker ? len1 + 3 : len1;
  if (linker)
    for (int i = 0; i < 3; i++) ps.ss.push_back(len1 + i);
  for (int i = 0; i < p2.bp_nr(); i++) ps.ss.push_back(off + p2.bp_left[i]);
  for (int i = 0; i < p2.bp_nr(); i++) ps.ss.push_back(off + p2.bp_right[i]);
  std::sort(ps.ss.begin(), ps.ss.end());
  return ps;
}

Alignment concatenate(const Alignment &a1, const Alignment &a2, const std::string &linker) {
  Alignment a;
  for (int i = 0; i < a1.nr(); i++) {
    a.ids.push_back(a1.ids[i]);
    a.seqs.push_back(a1.seqs[i] + linker + a2.seqs[i]);
  }
  a.raw = a.seqs;
  return a;
}

/* evolutionary scoring of the duplex (one tree for both alignments) */
void get_duplex_prob_tree(Context &ctx, const Alignment &a1, const Alignment &a2, Matrix &dp, std::vector<double> &ds,
                          const PartStruc &p1, const PartStruc &p2, std::string &tree) {
  const int l1 = a1.len(), l2 = a2.len();
  const Alignment align = concatenate(a1, a2, "nnn");
  const PartStruc ps = duplex_blocked(p1, p2, l1, true);
  Matrix tp;
  std::vector<double> ts;
  get_prob_tree(ctx, align, tp, ts, tree, &ps);
  dp = Matrix(l1 + l2, 0.0);
  ds.assign(l1 + l2, 0.0);
  if (!ctx.opt.use_evo) return;
  auto m = [&](int i) { return i < l1 ? i : i + 3; };  // skip linker
  for (int i = 0; i < l1 + l2; i++) {
    ds[i] = ts[m(i)];
    for (int j = 0; j < l1 + l2; j++) dp(i, j) = tp(m(i), m(j));
  }
}

/* thermodynamic scoring of the duplex -> RNAcofold */
void get_duplex_prob_energy(Context &ctx, const Alignment &a1, const Alignment &a2, Matrix &dp, std::vector<double> &ds,
                            const PartStruc &p1, const PartStruc &p2, bool legacy_cutpoint) {
  const Alignment align = concatenate(a1, a2, "");
  const PartStruc ps = duplex_blocked(p1, p2, a1.len(), false);
  dp = Matrix(align.len(), 0.0);
  get_cofold_prob_paired_seq(ctx, align, a1.len(), dp, &ps, legacy_cutpoint);
  ds = get_prob_unpaired_seq(ctx, align, dp);
}

void adjust_duplex_paired(Matrix &m, int len, double p) {
  for (int i = 0; i < len; i++)
    for (int j = i + 1; j < len; j++) m(i, j) *= p;
}
void adjust_duplex_single(std::vector<double> &v, double p) {
  for (double &x : v) x *= p;
}

/* MEA structure of the duplex with a 3-column linker */
std::string nussinov_duplex(const Matrix &dp, const std::vector<double> &ds, const PartStruc &ps, int len1, int len2,
                            float alpha) {
  const int n = len1 + len2 + 3;
  Matrix tp(n, 0.0);
  std::vector<double> ts(n, 0.0);
  auto m = [&](int i) { return i < len1 ? i : i + 3; };
  for (int i = 0; i < len1 + len2; i++) {
    ts[m(i)] = ds[i];
    for (int j = 0; j < len1 + len2; j++) tp(m(i), m(j)) = dp(i, j);
  }
  MeaTable T(tp, ts, &ps, alpha);
  return T.structure();
}

/* split the duplex structure and mark inter-molecular pairs with [ and ] */
void remove_linker(const std::string &db, std::string &db1, std::string &db2, int len1, int len2) {
  std::vector<int> stack;
  db1 = db.substr(0, len1);
  for (int i = 0; i < len1; i++) {
    if (db[i] == '(') stack.push_back(i);
    else if (db[i] == ')' && !stack.empty()) stack.pop_back();
  }
  for (int i : stack) db1[i] = '[';
  stack.clear();
  db2 = db.substr(len1 + 3, len2);
  for (int i = len2 - 1; i >= 0; i--) {
    if (db[i + len1 + 3] == ')') stack.push_back(i);
    else if (db[i + len1 + 3] == '(' && !stack.empty()) stack.pop_back();
  }
  for (int i : stack) db2[i] = ']';
}

double get_score(const std::string &db1, const std::string &db2, const Matrix &dp, const std::vector<double> &ds,
                 int len1, int len2, float alpha) {
  double score = 0;
  const int len = len1 + len2;
  std::vector<int> s1, s2, s3;
  for (int i = 0; i < len1; i++) {
    switch (db1[i]) {
      case '.': score += 2 * alpha * ds[i] / len; break;
      case '(': s1.push_back(i); break;
      case ')': score += 2 * dp(s1.back(), i) / len; s1.pop_back(); break;
      case '{': s2.push_back(i); break;
      case '}': score += 2 * dp(s2.back(), i) / len; s2.pop_back(); break;
      case '[': s3.push_back(i); break;
    }
  }
  s1.clear();
  s2.clear();
  for (int i = 0; i < len2; i++) {
    switch (db2[i]) {
      case '.': score += 2 * alpha * ds[i + len1] / len; break;
      case '(': s1.push_back(i); break;
      case ')': score += 2 * dp(s1.back() + len1, i + len1) / len; s1.pop_back(); break;
      case '{': s2.push_back(i); break;
      case '}': score += 2 * dp(s2.back() + len1, i + len1) / len; s2.pop_back(); break;
      case ']': score += 2 * dp(s3.back(), i + len1) / len; s3.pop_back(); break;
    }
  }
  return score;
}

void check_alignment(const Alignment &a, const char *msg) {
  for (int i = 1; i < a.nr(); i++)
    if (static_cast<int>(a.seqs[i].size()) != a.len()) throw Error(msg, true);
}

}  // namespace

PetcofoldResult run_petcofold(Context &ctx, const PetcofoldInput &in, const std::string &ppfile) {
  const Options &opt = ctx.opt;
  PetcofoldResult r;
  float petcon = in.delta, partprob = in.gamma;

  if (opt.verbose) {
    ctx.log.print(strf("PARAMETERS:\n maximal percent of gaps per column (--setgap) = %3.2f\n"
                       " minimal base pair reliability for evolutionary constraints (--setevocon_bp) = %3.2f\n"
                       " minimal single stranded reliability for evolutionary constraints (--setevocon_ss) = %3.2f\n"
                       " weighting factor for thermodynamic overlap (--setbeta) = %3.2f\n"
                       " weighting factor for single stranded reliabilities divided by 2 (--setalpha) = %3.2f\n"
                       " maximal allowed intra-molecular base-paired reliability to be accessible for RNA duplex "
                       "(--setdelta) = %3.2f\n"
                       " minimal probability of single partial structure (--setgamma) = %3.2f\n",
                       opt.gap, opt.evocon_bp, opt.evocon_ss, opt.beta, opt.alpha, petcon, partprob));
  }
  const Alignment &align1 = in.alignment1, &align2 = in.alignment2;
  check_alignment(align1, "The first file doesn't look like a FASTA file of aligned RNA sequences - sorry.");
  check_alignment(align2, "This second file doesn't look like a FASTA file of aligned RNA sequences - sorry.");
  if (opt.verbose)
    ctx.log.print(strf("INPUT:\n ALIGNMENT 1:\n number sequences = %i\n sequence length = %i\n"
                       " ALIGNMENT 2:\n number sequences = %i\n sequence length = %i\n",
                       align1.nr(), align1.len(), align2.nr(), align2.len()));

  Alignment gf1 = delete_gap_columns(ctx, align1, r.gap_columns1, opt.gap);
  Alignment gf2 = delete_gap_columns(ctx, align2, r.gap_columns2, opt.gap);
  std::vector<int> idx1, idx2;
  intersect_alignments(gf1, gf2, idx1, idx2);
  const Alignment sub1 = get_subalignment(gf1, idx1);
  if (sub1.nr() == 0) throw Error("Error <fasta>: Intersection of fasta files results in zero common identifiers!");
  if (sub1.nr() < 3)
    throw Error(strf("WARNING: PETcofold is cancelled as your input consists of ONLY %i sequence(s) with common "
                     "identifiers.\nPlease provide two alignments with at least 3 common identifiers.",
                     sub1.nr()),
                true);
  const Alignment sub2 = get_subalignment(gf2, idx2);
  const int len1 = sub1.len(), len2 = sub2.len();
  r.length1 = align1.len();
  r.length2 = align2.len();
  const int gapnr1 = align1.len() - len1, gapnr2 = align2.len() - len2;
  r.common_ids = sub1.ids;
  if (len1 == 0 || len2 == 0) throw Error("No alignment columns left after removing gap columns.");

  if (opt.verbose) {
    std::string s = " Common identifiers =";
    for (const auto &id : sub1.ids) s += " " + id + ",";
    s += "\n ALIGNMENT 1:\n gap columns = ";
    for (int c : r.gap_columns1) s += strf("%i, ", c + 1);
    s += "\n gap column free alignment = \n";
    for (int i = 0; i < sub1.nr(); i++) s += "  " + sub1.ids[i] + "\t" + sub1.seqs[i] + "\n";
    s += " ALIGNMENT 2:\n gap columns = ";
    for (int c : r.gap_columns2) s += strf("%i, ", c + 1);
    s += "\n gap column free alignment = \n";
    for (int i = 0; i < sub2.nr(); i++) s += "  " + sub2.ids[i] + "\t" + sub2.seqs[i] + "\n";
    ctx.log.print(s);
  }

  const std::string &settree = in.tree;
  std::string tree = settree;

  /* given structures */
  PartStruc setstructure1, setstructure2, pstruc;
  bool has_set1 = false, has_set2 = false, has_setjoint = false;
  if (!in.structure.empty()) {
    const std::string &ss = in.structure;
    petcon = partprob = 0.;
    /* FIX: 3.4 wrote the terminator one byte past the buffer and did not
       check the form of the duplex structure */
    if (static_cast<int>(ss.size()) != align1.len() + 1 + align2.len() || ss[align1.len()] != '&')
      throw Error(strf("Error <setstruc>: The duplex structure must have the form <structure1>&<structure2> with "
                       "lengths %i and %i.",
                       align1.len(), align2.len()));
    std::string t1 = ss.substr(0, align1.len()), t2 = ss.substr(align1.len() + 1);
    for (char &c : t1)
      if (c == '[' || c == ']') c = '.';
    for (char &c : t2)
      if (c == '[' || c == ']') c = '.';
    setstructure1 = mod_setstruct(t1, r.gap_columns1, align1.len(), gapnr1);
    setstructure2 = mod_setstruct(t2, r.gap_columns2, align2.len(), gapnr2);
    has_set1 = has_set2 = true;
    std::string inter = ss.substr(0, align1.len()) + "..." + ss.substr(align1.len() + 1);
    int n_inter = 0;
    for (int i = 0; i < align1.len(); i++) {
      if (inter[i] == '[') n_inter++;
      else inter[i] = '.';
    }
    for (int i = 0; i < align2.len(); i++) {
      if (inter[i + align1.len() + 3] == ']') n_inter--;
      else inter[i + align1.len() + 3] = '.';
    }
    if (n_inter)
      throw Error("Error <setstruc>: Number of left and right parenthesis are unequal for inter-molecular base pairs!");
    std::vector<int> gap_col = r.gap_columns1;
    for (int c : r.gap_columns2) gap_col.push_back(c + align1.len() + 3);
    pstruc = mod_setstruct(inter, gap_col, align1.len() + align2.len() + 3, gapnr1 + gapnr2);
    has_setjoint = true;
  }

  /* step 1: PETfold for both alignments */
  Matrix paired_pet1, paired_pet2;
  std::vector<double> single_pet1, single_pet2;

  ctx.log.print("PETfold Alignment 1:\n");
  if (!in.structure1.empty() && !has_setjoint) {
    setstructure1 = mod_setstruct(in.structure1, r.gap_columns1, align1.len(), gapnr1);
    has_set1 = true;
  }
  const PartStruc struc_pet1 =
      run_petfold_step(ctx, sub1, paired_pet1, single_pet1, tree, has_set1 ? &setstructure1 : nullptr);
  if (opt.verbose) ctx.log.print("PARTIAL STRUCTURE 1:\n");
  const PartStruc ss_partial1 = get_probable_partial_struc(ctx, sub1, paired_pet1, struc_pet1, tree, petcon, partprob,
                                                           in.extstem, r.partial_prob1, r.delta1);
  r.partial_structure1 = include_gap_columns(partstruc2string(&ss_partial1, len1), r.gap_columns1, align1.len());
  ctx.log.print("Partial Structure 1:\t" + r.partial_structure1 + "\n");

  ctx.log.print("PETfold Alignment 2:\n");
  if (!in.structure2.empty() && !has_setjoint) {
    setstructure2 = mod_setstruct(in.structure2, r.gap_columns2, align2.len(), gapnr2);
    has_set2 = true;
  }
  /* each alignment gets its own tree unless a complete tree was given
     (FIX: 3.4 cleared the buffer with memset(tree, 0, sizeof(pointer))) */
  if (settree.empty()) tree.clear();
  else if (settree.find(':') == std::string::npos) tree = settree;
  const PartStruc struc_pet2 =
      run_petfold_step(ctx, sub2, paired_pet2, single_pet2, tree, has_set2 ? &setstructure2 : nullptr);
  if (opt.verbose) ctx.log.print("PARTIAL STRUCTURE 2:\n");
  const PartStruc ss_partial2 = get_probable_partial_struc(ctx, sub2, paired_pet2, struc_pet2, tree, petcon, partprob,
                                                           in.extstem, r.partial_prob2, r.delta2);
  r.partial_structure2 = include_gap_columns(partstruc2string(&ss_partial2, len2), r.gap_columns2, align2.len());
  ctx.log.print("Partial Structure 2:\t" + r.partial_structure2 + "\n");

  /* step 2: PETcofold for the RNA duplex */
  if (opt.verbose) ctx.log.print("COFOLDING:\n");
  if (settree.empty()) tree.clear();
  else if (settree.find(':') == std::string::npos) tree = settree;

  const int L = len1 + len2;
  Matrix dpt, dpe;
  std::vector<double> dst, dse;
  get_duplex_prob_tree(ctx, sub1, sub2, dpt, dst, ss_partial1, ss_partial2, tree);
  adjust_duplex_paired(dpt, L, geometric_mean(r.partial_prob1[1], r.partial_prob2[1]));
  adjust_duplex_single(dst, geometric_mean(r.partial_prob1[1], r.partial_prob2[1]));
  get_duplex_prob_energy(ctx, sub1, sub2, dpe, dse, ss_partial1, ss_partial2, in.legacy_cutpoint);
  adjust_duplex_paired(dpe, L, geometric_mean(r.partial_prob1[0], r.partial_prob2[0]));
  adjust_duplex_single(dse, geometric_mean(r.partial_prob1[0], r.partial_prob2[0]));

  std::vector<double> dsp = get_pet_single(dst, dse, opt.beta);
  Matrix dpp = get_pet_paired(dpt, dpe, opt.beta);
  for (int i = 0; i < ss_partial1.bp_nr(); i++) {
    dsp[ss_partial1.bp_left[i]] = 0.;
    dsp[ss_partial1.bp_right[i]] = 0.;
  }
  for (int i = 0; i < ss_partial2.bp_nr(); i++) {
    dsp[len1 + ss_partial2.bp_left[i]] = 0.;
    dsp[len1 + ss_partial2.bp_right[i]] = 0.;
  }
  for (int i = 0; i < ss_partial1.bp_nr(); i++)
    dpp(ss_partial1.bp_left[i], ss_partial1.bp_right[i]) = add(r.partial_prob1[1], opt.beta * r.partial_prob1[0]) / 2;
  for (int i = 0; i < ss_partial2.bp_nr(); i++)
    dpp(len1 + ss_partial2.bp_left[i], len1 + ss_partial2.bp_right[i]) =
        add(r.partial_prob2[1], opt.beta * r.partial_prob2[0]) / 2;

  if (!has_setjoint) pstruc = duplex_blocked(ss_partial1, ss_partial2, len1, true);
  const std::string db = nussinov_duplex(dpp, dsp, pstruc, len1, len2, opt.alpha);
  std::string db1, db2;
  remove_linker(db, db1, db2, len1, len2);
  for (int i = 0; i < ss_partial1.bp_nr(); i++) {
    db1[ss_partial1.bp_left[i]] = '{';
    db1[ss_partial1.bp_right[i]] = '}';
  }
  for (int i = 0; i < ss_partial2.bp_nr(); i++) {
    db2[ss_partial2.bp_left[i]] = '{';
    db2[ss_partial2.bp_right[i]] = '}';
  }
  r.structure1 = include_gap_columns(db1, r.gap_columns1, align1.len());
  r.structure2 = include_gap_columns(db2, r.gap_columns2, align2.len());

  if (!ppfile.empty()) {
    std::vector<int> gap_col = r.gap_columns1;
    for (int c : r.gap_columns2) gap_col.push_back(c + align1.len());
    create_ppfile(dpp, dsp, gap_col, L, gapnr1 + gapnr2, ppfile);
  }
  r.score = get_score(db1, db2, dpp, dsp, len1, len2, opt.alpha);

  r.intermolecular1 = r.structure1;
  for (char &c : r.intermolecular1)
    if (c != '.' && c != '-' && c != '[') c = '.';
  r.intermolecular2 = r.structure2;
  for (char &c : r.intermolecular2)
    if (c != '.' && c != '-' && c != ']') c = '.';

  r.has_evo = opt.use_evo;
  auto no_x = [](std::string s) {
    for (char &c : s)
      if (c == 'x') c = '.';
    return s;
  };
  r.petfold_structure1 = include_gap_columns(no_x(partstruc2string(&struc_pet1, len1)), r.gap_columns1, align1.len());
  r.petfold_structure2 = include_gap_columns(no_x(partstruc2string(&struc_pet2, len2)), r.gap_columns2, align2.len());

  /* FIX: 3.4 paired the sequences of the two files by position; they are
     now paired by identifier (common identifiers in the order of file 1) */
  for (const auto &id : r.common_ids) {
    auto find = [&](const Alignment &a) {
      for (int i = 0; i < a.nr(); i++)
        if (a.ids[i] == id) return a.raw[i];
      return std::string();
    };
    r.fasta_ids.push_back(id);
    r.fasta_seq1.push_back(find(align1));
    r.fasta_seq2.push_back(find(align2));
  }
  r.duplex_paired_pet = std::move(dpp);
  r.duplex_single_pet = std::move(dsp);
  return r;
}

}  // namespace petfold
