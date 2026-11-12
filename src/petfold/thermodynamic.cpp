/*
 * thermodynamic.cpp - thermodynamic model via the ViennaRNA 2.x API.
 *
 * PETfold 2.2 used the deprecated global-variable interface (pf_fold, pr,
 * iindx, cut_point, pf_scale, ...), which is not thread safe.  This file uses
 * vrna_fold_compound_t objects instead but reproduces exactly what the
 * deprecated wrappers did (model settings, constraint options, handling of
 * the pf_scale scaling factor), so the probabilities are unchanged.
 */
#include <cmath>
#include <memory>

extern "C" {
#include <ViennaRNA/constraints/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/part_func_co.h>
}

#include "core.hpp"

namespace petfold {

namespace {

constexpr int SCALELENGTH = 200;
constexpr double SFACT = 1.07;

struct FcDeleter {
  void operator()(vrna_fold_compound_t *fc) const { vrna_fold_compound_free(fc); }
};
using FoldCompound = std::unique_ptr<vrna_fold_compound_t, FcDeleter>;

/* model details of the deprecated interface with dangles = 2, noLP = 1 */
vrna_md_t legacy_md(bool dimer) {
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.dangles = 2;
  md.noLP = 1;
  md.circ = 0;
  md.compute_bpp = 1;
  if (dimer) md.min_loop_size = 0;  // as in the deprecated co-folding wrappers
  return md;
}

constexpr unsigned int kConstraintOptions = VRNA_CONSTRAINT_DB | VRNA_CONSTRAINT_DB_PIPE | VRNA_CONSTRAINT_DB_DOT |
                                            VRNA_CONSTRAINT_DB_X | VRNA_CONSTRAINT_DB_ANG_BRACK |
                                            VRNA_CONSTRAINT_DB_RND_BRACK;

double kT_kcal() {
  vrna_md_t md;
  vrna_md_set_default(&md);
  return (md.temperature + 273.15) * 1.98717 / 1000.;  // kT in kcal/mol
}

/* equivalent of fold(): unconstrained MFE */
double mfe_single(const std::string &seq) {
  vrna_md_t md = legacy_md(false);
  FoldCompound fc(vrna_fold_compound(seq.c_str(), &md, VRNA_OPTION_DEFAULT));
  if (!fc) throw Error("ViennaRNA: cannot create fold compound");
  return vrna_mfe(fc.get(), nullptr);
}

/* equivalent of cofold(): MFE of a dimer, min loop size 0 */
double mfe_dimer(const std::string &seq_with_cut) {
  vrna_md_t md = legacy_md(true);
  FoldCompound fc(vrna_fold_compound(seq_with_cut.c_str(), &md, 0));
  if (!fc) throw Error("ViennaRNA: cannot create fold compound");
  return vrna_mfe(fc.get(), nullptr);
}

/* equivalent of pf_fold(seq, constraint) with fold_constrained = 1 (or 0 if
   constraint is null); returns the fold compound for reading probabilities */
FoldCompound pf_single(Context &ctx, const std::string &seq, const std::string *constraint, float &energy,
                       std::string *symbols = nullptr) {
  vrna_md_t md = legacy_md(false);
  FoldCompound fc(vrna_fold_compound(seq.c_str(), &md, VRNA_OPTION_DEFAULT));
  if (!fc) throw Error("ViennaRNA: cannot create fold compound");
  if (fc->exp_params) free(fc->exp_params);
  fc->exp_params = vrna_exp_params(&(fc->params->model_details));
  fc->exp_params->pf_scale = ctx.pf_scale;
  if (constraint) vrna_constraints_add(fc.get(), constraint->c_str(), kConstraintOptions);
  /* 2.2 passed the constraint buffer, which vrna_pf overwrites with the pair
     probability symbol string; --verbose printed that string */
  std::string buf(seq.size() + 2, '\0');
  energy = vrna_pf(fc.get(), buf.data());
  if (symbols) *symbols = buf.c_str();
  return fc;
}

void maybe_rescale(Context &ctx, const std::string &seq, double min_en, int len) {
  (void)seq;
  const double kT = kT_kcal();
  ctx.pf_scale = std::exp(-(SFACT * min_en) / kT / len);
}

}  // namespace

Matrix get_prob_paired_seq(Context &ctx, const Alignment &align_gf, const PartStruc *evocon) {
  const int len = align_gf.len();
  Matrix paired(len, 0.0);
  const SeqList sl = get_align_without_gaps(align_gf);
  if (ctx.opt.verbose) ctx.log.print("RNAFOLD:\n");
  if (!ctx.opt.use_thermo) return paired;

  for (int i = 0; i < align_gf.nr(); i++) {
    const int n = static_cast<int>(sl.seq[i].size());
    const std::string cons = get_constraint_string(evocon, sl.origid[i], n);
    if (ctx.opt.verbose) ctx.log.print(" Constrained structure = " + cons + "\n");
    if (n > SCALELENGTH) maybe_rescale(ctx, sl.seq[i], mfe_single(sl.seq[i]), n);
    float e = 0;
    std::string symbols;
    FoldCompound fc = pf_single(ctx, sl.seq[i], &cons, e, &symbols);
    const FLT_OR_DBL *pr = fc->exp_matrices->probs;
    const int *iindx = fc->iindx;
    for (int d = 1; d <= n; d++)
      for (int k = d + 1; k <= n; k++) paired(sl.origid[i][d - 1], sl.origid[i][k - 1]) += pr[iindx[d] - k];
    if (ctx.opt.verbose)
      ctx.log.print(strf(" Sequence %2i structure = %s\n free energy of ensemble=%5.2f kcal/mol\n", i + 1,
                         symbols.c_str(), e));
  }
  /* divide the sum of base pair probabilities by the number of gap-free base pairs */
  for (int d = 0; d < len; d++)
    for (int k = d + 1; k < len; k++) {
      int gappos = 0;
      for (int m = 0; m < align_gf.nr(); m++)
        if (align_gf.seqs[m][d] == '-' || align_gf.seqs[m][k] == '-') gappos++;
      /* FIX: 2.2 divided by zero (NaN) if every sequence had a gap at d or k */
      if (align_gf.nr() - gappos > 0) paired(d, k) /= align_gf.nr() - gappos;
    }
  return paired;
}

std::vector<double> get_prob_unpaired_seq(const Context &ctx, const Alignment &align_gf, const Matrix &paired_seq) {
  const int len = align_gf.len();
  std::vector<double> single(len, 0.0);
  if (!ctx.opt.use_thermo) return single;
  for (int i = 0; i < len; i++) {
    double sum = 0;
    for (int j = i + 1; j < len; j++) sum += paired_seq(i, j);
    for (int j = 0; j < i; j++) sum += paired_seq(j, i);
    single[i] = 1 - sum;
  }
  return single;
}

double get_fold_thermo_partial_prob(double econstraint, double eensemble) {
  return std::exp((eensemble - econstraint) / kT_kcal());
}

double get_constraint_energy(Context &ctx, const Alignment &align, const PartStruc *evocon) {
  const SeqList sl = get_align_without_gaps(align);
  double e = 0;
  for (int i = 0; i < align.nr(); i++) {
    const int n = static_cast<int>(sl.seq[i].size());
    const std::string cons = get_constraint_string(evocon, sl.origid[i], n);
    if (n > SCALELENGTH) maybe_rescale(ctx, sl.seq[i], mfe_single(sl.seq[i]), n);
    float ei = 0;
    pf_single(ctx, sl.seq[i], &cons, ei);
    e += static_cast<double>(ei);
  }
  return e / align.nr();
}

double get_ensemble_energy(Context &ctx, const Alignment &align) {
  const SeqList sl = get_align_without_gaps(align);
  double e = 0.;
  for (int i = 0; i < align.nr(); i++) {
    const int n = static_cast<int>(sl.seq[i].size());
    if (n > SCALELENGTH) maybe_rescale(ctx, sl.seq[i], mfe_single(sl.seq[i]), n);
    float ei = 0;
    pf_single(ctx, sl.seq[i], nullptr, ei);
    e += static_cast<double>(ei);
  }
  return e / align.nr();
}

void get_cofold_prob_paired_seq(Context &ctx, const Alignment &align, int len1, Matrix &paired_seq,
                                const PartStruc *pstruc, bool legacy_cutpoint) {
  const SeqList sl = get_align_without_gaps(align);
  if (ctx.opt.verbose) ctx.log.print("RNACOFOLD:\n");
  if (!ctx.opt.use_thermo) return;

  for (int i = 0; i < align.nr(); i++) {
    const int n = static_cast<int>(sl.seq[i].size());
    const std::string cons = get_constraint_string(pstruc, sl.origid[i], n);
    if (ctx.opt.verbose) ctx.log.print(" Constrained structure = " + cons + "\n");
    /* position (0-based) of the first nucleotide of the second molecule */
    int cut;
    if (legacy_cutpoint) {
      cut = len1;  // 3.4: cut_point = len1 + 1 for every sequence, ignoring its gaps
    } else {
      /* FIX: cut at the true boundary = number of nucleotides of molecule 1 */
      cut = 0;
      while (cut < n && sl.origid[i][cut] < len1) cut++;
    }
    std::string seq = sl.seq[i];
    if (cut > 0 && cut < n) seq.insert(static_cast<std::size_t>(cut), 1, '&');
    if (n > SCALELENGTH) maybe_rescale(ctx, seq, mfe_dimer(seq), n);

    vrna_md_t md = legacy_md(true);
    FoldCompound fc(vrna_fold_compound(seq.c_str(), &md, VRNA_OPTION_DEFAULT));
    if (!fc) throw Error("ViennaRNA: cannot create fold compound");
    if (fc->exp_params) free(fc->exp_params);
    fc->exp_params = vrna_exp_params(&(fc->params->model_details));
    fc->exp_params->pf_scale = ctx.pf_scale;
    vrna_constraints_add(fc.get(), cons.c_str(), kConstraintOptions);
    std::string symbols(seq.size() + 2, '\0');
    const vrna_dimer_pf_t e = vrna_pf_dimer(fc.get(), symbols.data());

    const FLT_OR_DBL *pr = fc->exp_matrices->probs;
    const int *iindx = fc->iindx;
    for (int d = 1; d <= n; d++)
      for (int k = d + 1; k <= n; k++) paired_seq(sl.origid[i][d - 1], sl.origid[i][k - 1]) += pr[iindx[d] - k];
    if (ctx.opt.verbose)
      ctx.log.print(strf(" Sequence %2i structure = %s\n free energy of ensemble: F0AB=%5.2f kcal/mol; FA=%5.2f "
                         "kcal/mol; FAB=%5.2f kcal/mol; FB=%5.2f kcal/mol; FcAB=%5.2f kcal/mol\n",
                         i + 1, symbols.c_str(), e.F0AB, e.FA, e.FAB, e.FB, e.FcAB));
  }
  /* arithmetic mean of base pair probabilities over sequences without gap in column d */
  std::vector<int> gappos(align.len(), 0);
  for (int d = 0; d < align.len(); d++)
    for (int k = 0; k < align.nr(); k++)
      if (align.seqs[k][d] == '-') gappos[d]++;
  for (int d = 0; d < align.len(); d++)
    for (int k = d + 1; k < align.len(); k++)
      if (align.nr() - gappos[d] > 0) paired_seq(d, k) /= align.nr() - gappos[d];
}

}  // namespace petfold
