/*
 * petfold.cpp - the PETfold pipeline (from petfold.c).
 */
#include "api.hpp"
#include "core.hpp"

namespace petfold {

namespace {
void check_alignment(const Alignment &a, const char *msg) {
  for (int i = 1; i < a.nr(); i++)
    if (static_cast<int>(a.seqs[i].size()) != a.len()) throw Error(msg, true);
}
}  // namespace

PetfoldResult run_petfold(Context &ctx, const PetfoldInput &in) {
  const Options &opt = ctx.opt;
  const Alignment &align = in.alignment;
  PetfoldResult r;

  check_alignment(align, "This doesn't look like a FASTA file of aligned RNA sequences - sorry.");
  if (opt.verbose)
    ctx.log.print(strf("INPUT:\n number sequences = %i\n sequence length = %i\n", align.nr(), align.len()));
  if (align.nr() < 3)
    throw Error(strf("WARNING: PETfold is cancelled as your input consists of ONLY %i sequence(s).\n"
                     "Please provide an alignment of at least 3 sequences.",
                     align.nr()),
                true);

  r.alignment_length = align.len();
  r.alignment_gf = delete_gap_columns(ctx, align, r.gap_columns, opt.gap);
  const Alignment &gf = r.alignment_gf;
  const int len = gf.len();
  const int gapnr = align.len() - len;
  if (opt.verbose) {
    std::string s = strf(" maximal percent of gaps per column = %3.2f\n", opt.gap);
    s += " gap columns = ";
    for (int c : r.gap_columns) s += strf("%i, ", c + 1);
    s += "\n gap column free alignment = \n";
    for (int i = 0; i < gf.nr(); i++) s += "  " + gf.ids[i] + "\t" + gf.seqs[i] + "\n";
    ctx.log.print(s);
  }
  if (len == 0 || gf.nr() == 0) throw Error("No alignment columns left after removing gap columns.");

  /* evolutionary model -> Pfold (or reliabilities computed by PPfold) */
  r.tree = in.tree;
  std::vector<int> pfold_struc;
  if (in.ppfold_input)
    pfold_struc = get_ppfold(len, r.paired_tree, r.single_tree, in.ppfold_input);
  else
    pfold_struc = get_prob_tree(ctx, gf, r.paired_tree, r.single_tree, r.tree);

  PartStruc evocon;
  const PartStruc *evocon_p = nullptr;
  r.has_evo = opt.use_evo;
  if (opt.use_evo) {
    r.pfold_structure_gf = get_dot_bracket(pfold_struc);
    evocon = get_partial_struc(r.paired_tree, r.single_tree, pfold_struc, len, opt.evocon_bp, opt.evocon_ss);
    evocon_p = &evocon;
  }
  r.constraints_gf = partstruc2string(evocon_p, len);

  /* thermodynamic model -> RNAfold */
  r.paired_seq = get_prob_paired_seq(ctx, gf, evocon_p);
  r.single_seq = get_prob_unpaired_seq(ctx, gf, r.paired_seq);

  /* joint reliabilities */
  r.single_pet = get_pet_single(r.single_tree, r.single_seq, opt.beta);
  r.paired_pet = get_pet_paired(r.paired_tree, r.paired_seq, opt.beta);

  /* structure to score */
  PartStruc setstruc;
  const PartStruc *dp_constraint = evocon_p;
  if (!in.structure.empty()) {
    setstruc = mod_setstruct(in.structure, r.gap_columns, align.len(), gapnr);
    dp_constraint = &setstruc;
  }

  /* MEA structure (Nussinov-style recursion with sums instead of products) */
  MeaTable T(r.paired_pet, r.single_pet, dp_constraint, opt.alpha);
  r.structure_gf = T.structure();
  r.score = T.score();
  r.ensemble_diversity = ensemble_diversity(r.paired_pet, len);
  r.structure = include_gap_columns(r.structure_gf, r.gap_columns, align.len());

  if (in.suboptimal > 0) {
    SuboptSampler sampler(r.paired_pet, r.single_pet, opt.alpha);
    std::string s;
    for (int i = 0; i < in.suboptimal; i++) {
      const double sc = sampler.sample(ctx.rng, s);
      r.suboptimals.emplace_back(include_gap_columns(s, r.gap_columns, align.len()), sc);
    }
  }
  return r;
}

}  // namespace petfold
