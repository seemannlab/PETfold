/*
 * evolutionary.cpp - evolutionary model: drives the Pfold engine
 * (from evolutionary.c).
 */
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>

#include "../pfold/inout.hpp"
#include "../pfold/newcolprob.hpp"
#include "../pfold/pf_alloc.hpp"
#include "../pfold/phyl.hpp"
#include "core.hpp"

#ifndef PETFOLD_DATA_DIR
#define PETFOLD_DATA_DIR ""
#endif

namespace petfold {

using namespace pfold;

namespace {

using pfold::Matrix;

/* copied from Pfold:mltree.c */
struct NodeInfo {
  int seqnum;
  Matrix *pmatrix;
  Matrix **p_top_up;
  Matrix **p_top_down;
  Matrix **p_bot_up;
  Matrix **p_bot_down;
  Grammar *grammar;
};

int find_sym_cb(char c, void *g) { return FindSym(c, static_cast<Grammar *>(g)); }

Grammar *read_grammar(const Context &ctx, const char *name) {
  const std::string file = grammar_file(ctx, name);
  std::FILE *fp = std::fopen(file.c_str(), "r");
  if (!fp) throw petfold::Error("Pfold: Error in opening file '" + file + "'");
  Grammar *g = nullptr;
  try {
    g = ReadGrammar(fp);
  } catch (...) {
    std::fclose(fp);
    throw;
  }
  std::fclose(fp);
  return g;
}

/* Makes an entry, ready for filling out (adapted from Pfold:col.c) */
Entry *MyNewEntry(const char *type, const char *name, int len) {
  Entry *entry = MakeEntry();
  entry->textlen = 7;
  entry->text = (char **)pf_malloc(entry->textlen * sizeof(char *));
  const std::size_t tl = std::strlen(type), nl = std::strlen(name);
  auto line = [&](int k, std::size_t size, const std::string &s) {
    entry->text[k] = (char *)pf_malloc(size);
    std::snprintf(entry->text[k], size, "%s", s.c_str());
  };
  line(0, 32 + tl, std::string("; TYPE              ") + type + "\n");
  line(1, 32, "; COL 1             label\n");
  line(2, 32, "; COL 2             residue\n");
  line(3, 32, "; COL 3             seqpos\n");
  line(4, 32, "; COL 4             alignpos\n");
  line(5, 32 + nl, std::string("; ENTRY             ") + name + "\n");
  line(6, 16, "; ----------\n");
  entry->len = len;
  entry->line = (char **)pf_malloc(len * sizeof(char *));
  for (int i = 0; i < len; i++) {
    /* FIX: 2.2 allocated 16 bytes per line, which overflowed for alignments
       with more than 9999 columns */
    entry->line[i] = (char *)pf_malloc(48);
    std::strcpy(entry->line[i], "\n");
  }
  entry->endline = (char *)pf_malloc(16);
  std::strcpy(entry->endline, "; **********\n");
  return entry;
}

/* Converts an alignment into a NULL terminated list of col-format entries */
Entry **Aln2EntryList(const petfold::Alignment &align) {
  Entry **list = (Entry **)pf_malloc((align.nr() + 1) * sizeof(Entry *));
  for (int i = 0; i < align.nr(); i++) {
    Entry *entry = MyNewEntry("RNA", align.ids[i].c_str(), align.len());
    int k = 0;
    for (int j = 0; j < align.len(); j++) {
      const char c = align.seqs[i][j];
      if (c != '-')
        std::snprintf(entry->line[j], 48, "N %c %i %i\n", c, ++k, j + 1);
      else
        std::snprintf(entry->line[j], 48, "N %c . %i\n", c, j + 1);
    }
    list[i] = entry;
  }
  list[align.nr()] = nullptr;
  return list;
}

void newick_rec(std::string &out, PhylNode *pnode) {
  if (pnode->child == nullptr) {
    out += static_cast<char *>(pnode->elm);
  } else {
    out += '(';
    if (pnode->elm != nullptr) {
      out += static_cast<char *>(pnode->elm);
      out += ',';
    }
    for (PhylNode *c = pnode->child; c != nullptr; c = c->brother) {
      newick_rec(out, c);
      if (c->brother != nullptr) out += ',';
    }
    out += ')';
  }
  if (pnode->uplen != 0) out += strf(":%g", pnode->uplen);
}

/* Newick string of a tree, same format as PrintPhyl/Phyl2Char */
std::string phyl_to_newick(Phyl *phyl) {
  std::string s;
  newick_rec(s, phyl->root);
  return s + ";";
}

/* Sets alignment names and numbers (copied from Pfold:mltree.c) */
void initphyl(PhylNode *pnode, Grammar *grammar, Align *align) {
  const int numterm = grammar->term->count;
  char *name = static_cast<char *>(pnode->elm);
  pnode->elm = pf_malloc(sizeof(NodeInfo));
  if (pnode->uplen < 0.001) pnode->uplen = 0.001;
  NodeInfo *pn = static_cast<NodeInfo *>(pnode->elm);
  const int len = align->len;
  pn->seqnum = SeqNumber(align, name);
  pn->grammar = grammar;
  pn->pmatrix = nullptr;
  pn->p_bot_down = (Matrix **)pf_malloc(len * sizeof(Matrix *));
  for (int i = 0; i < len; i++) pn->p_bot_down[i] = MakeMatrix(1, numterm);
  pn->p_top_down = (Matrix **)pf_malloc(len * sizeof(Matrix *));
  for (int i = 0; i < len; i++) pn->p_top_down[i] = nullptr;
  pn->p_bot_up = (Matrix **)pf_malloc(len * sizeof(Matrix *));
  for (int i = 0; i < len; i++) pn->p_bot_up[i] = nullptr;
  pn->p_top_up = (Matrix **)pf_malloc(len * sizeof(Matrix *));
  for (int i = 0; i < len; i++) pn->p_top_up[i] = MakeMatrix(1, numterm);
  for (PhylNode *c = pnode->child; c != nullptr; c = c->brother) initphyl(c, grammar, align);
}

Sgrp *first_sgrp(Grammar *grammar) {
  LListCounter *lc = MakeCounter(grammar->sgrp, FIRST);
  Sgrp *s = static_cast<Sgrp *>(Next(lc));
  pf_free(lc);
  return s;
}

/* Calculates evolutionary matrices (copied from Pfold:mltree.c) */
void setmat(PhylNode *pnode, Grammar *grammar) {
  Sgrp *sgrp = first_sgrp(grammar);
  NodeInfo *pn = static_cast<NodeInfo *>(pnode->elm);
  FreeMatrix(pn->pmatrix);
  Matrix *temp = ExpMatrix(pnode->uplen, sgrp->eigen, sgrp->diag, sgrp->inveigen);
  pn->pmatrix = TransposeMatrix(temp);
  FreeMatrix(temp);
  for (PhylNode *c = pnode->child; c != nullptr; c = c->brother) setmat(c, grammar);
}

/* probability distributions downwards in the tree (copied from Pfold:mltree.c) */
void downcalc(PhylNode *pnode, Align *align) {
  NodeInfo *pn = static_cast<NodeInfo *>(pnode->elm);
  const int numterm = pn->grammar->term->count;
  for (PhylNode *c = Child(pnode); c != nullptr; c = Brother(c)) downcalc(c, align);
  for (int pos = 0; pos < align->len; pos++) {
    if (pn->seqnum == -1) {
      for (int i = 0; i < numterm; i++) pn->p_bot_down[pos]->entry[0][i] = Dbl2Edbl(1.);
    } else {
      const int symno = align->seq[pn->seqnum][pos];
      for (int i = 0; i < numterm; i++) pn->p_bot_down[pos]->entry[0][i] = pn->grammar->quickdist[symno][i];
    }
    for (PhylNode *c = Child(pnode); c != nullptr; c = Brother(c))
      for (int i = 0; i < numterm; i++)
        MulEdouble(&pn->p_bot_down[pos]->entry[0][i],
                   static_cast<NodeInfo *>(c->elm)->p_top_down[pos]->entry[0][i]);
    FreeMatrix(pn->p_top_down[pos]);
    pn->p_top_down[pos] = MulMatrix(pn->p_bot_down[pos], pn->pmatrix);
  }
}

/* probability distributions upwards in the tree (copied from Pfold:mltree.c) */
void upcalc(PhylNode *pnode, Align *align) {
  NodeInfo *pn = static_cast<NodeInfo *>(pnode->elm);
  const int numterm = pn->grammar->term->count;
  for (int pos = 0; pos < align->len; pos++) {
    for (int i = 0; i < numterm; i++) pn->p_top_up[pos]->entry[0][i] = Dbl2Edbl(1.);
    if (pnode->parent != nullptr) {
      for (int i = 0; i < numterm; i++)
        MulEdouble(&pn->p_top_up[pos]->entry[0][i],
                   static_cast<NodeInfo *>(pnode->parent->elm)->p_bot_up[pos]->entry[0][i]);
      for (PhylNode *c = Child(pnode->parent); c != nullptr; c = Brother(c))
        if (c != pnode)
          for (int i = 0; i < numterm; i++)
            MulEdouble(&pn->p_top_up[pos]->entry[0][i],
                       static_cast<NodeInfo *>(c->elm)->p_top_down[pos]->entry[0][i]);
    }
    FreeMatrix(pn->p_bot_up[pos]);
    pn->p_bot_up[pos] = MulMatrix(pn->p_top_up[pos], pn->pmatrix);
  }
  for (PhylNode *c = Child(pnode); c != nullptr; c = Brother(c)) upcalc(c, align);
}

/* branch length optimisation (copied from Pfold:mltree.c) */
double optimize(const Context &ctx, PhylNode *pnode, Align *align, Grammar *grammar, double *collikelhd) {
  double nll = 0;
  for (PhylNode *c = Child(pnode); c != nullptr; c = Brother(c)) nll = optimize(ctx, c, align, grammar, collikelhd);
  if (pnode->parent == nullptr) return nll;

  Sgrp *sgrp = first_sgrp(grammar);
  NodeInfo *pn = static_cast<NodeInfo *>(pnode->elm);
  const int numterm = pn->grammar->term->count;
  double x = std::log(pnode->uplen + 0.001);
  double l = std::exp(x);
  int finish = 0, error;
  InitMinimize(x - 0.2, x + 0.2, 0.001);
  while ((error = Minimize(&x, &nll)) == 0) {
    if (x < -10) {
      x = -10;
      finish = 10;
    } else if (x > 6) {
      x = 6;
      finish++;
    }
    if (finish > 4) break;
    l = std::exp(x);
    nll = 0;
    Matrix *temp = ExpMatrix(l, sgrp->eigen, sgrp->diag, sgrp->inveigen);
    Matrix *p_temp = TransposeMatrix(temp);
    FreeMatrix(temp);
    for (int pos = 0; pos < align->len; pos++) {
      temp = MulMatrix(pn->p_bot_down[pos], p_temp);
      Edouble sum = Dbl2Edbl(0.);
      for (int i = 0; i < numterm; i++) {
        Edouble e = ProdEdouble(sgrp->freq->entry[0][i], ProdEdouble(temp->entry[0][i], pn->p_top_up[pos]->entry[0][i]));
        AddEdouble(&sum, e);
      }
      FreeMatrix(temp);
      nll -= Edbl2Dbl(LogEdouble(sum));
      collikelhd[pos] = Edbl2Dbl(sum);
    }
    FreeMatrix(p_temp);
  }
  if (error == 1) throw petfold::Error("Pfold: problems in branch length optimisation");
  (void)ctx;

  Matrix *temp = ExpMatrix(l, sgrp->eigen, sgrp->diag, sgrp->inveigen);
  Matrix *p_temp = TransposeMatrix(temp);
  FreeMatrix(temp);
  pnode->uplen = std::exp(x);
  if (pnode->parent != nullptr) {
    for (int pos = 0; pos < align->len; pos++) {
      Matrix *old = pn->p_top_down[pos];
      Matrix *nw = MulMatrix(pn->p_bot_down[pos], p_temp);
      Matrix *mat = static_cast<NodeInfo *>(pnode->parent->elm)->p_bot_down[pos];
      for (int i = 0; i < numterm; i++) {
        if (Edbl2Dbl(old->entry[0][i]) != 0) DivEdouble(&mat->entry[0][i], old->entry[0][i]);
        MulEdouble(&mat->entry[0][i], nw->entry[0][i]);
      }
      for (PhylNode *c = Child(pnode->parent); c != nullptr; c = Brother(c))
        if (c != pnode) {
          mat = static_cast<NodeInfo *>(c->elm)->p_top_up[pos];
          for (int i = 0; i < numterm; i++) {
            if (Edbl2Dbl(old->entry[0][i]) != 0) DivEdouble(&mat->entry[0][i], old->entry[0][i]);
            MulEdouble(&mat->entry[0][i], nw->entry[0][i]);
          }
        }
      FreeMatrix(nw);
    }
  }
  FreeMatrix(p_temp);
  return nll;
}

/* Structure constraints for the column probabilities (Alignstr.str) */
void init_alignstr_str(Alignstr *alignstr, Entry **entry_list, const PartStruc *constraint, int len,
                       Grammar *grammar, bool singles_forced) {
  alignstr->str = (int *)pf_malloc(len * sizeof(int));
  if (constraint == nullptr || (constraint->ss_nr() == 0 && constraint->bp_nr() == 0)) {
    char field[MAXCOLW];
    const int align_bp_col = ReadColno(entry_list[0], "align_bp");
    const int alignpos_col = ReadColno(entry_list[0], "alignpos");
    const int seq_bp_col = ReadColno(entry_list[0], "seq_bp");
    const int seqpos_col = ReadColno(entry_list[0], "seqpos");
    if (align_bp_col == 0 && seq_bp_col == 0) {
      for (int i = 0; i < len; i++) alignstr->str[i] = -1;
    } else {
      for (int i = 0; i < len; i++) {
        GetField(field, entry_list[0], i + 1, align_bp_col != 0 ? align_bp_col : seq_bp_col);
        if (StrCmp(field, "s") == 0) alignstr->str[i] = SINGLE;
        else if (StrCmp(field, "d") == 0) alignstr->str[i] = DOUBLE;
        else alignstr->str[i] = FindPair(entry_list[0], i + 1, align_bp_col, alignpos_col, seq_bp_col, seqpos_col) - 1;
      }
    }
  } else {
    const std::vector<int> struc = partstruc2str(constraint, len);
    for (int i = 0; i < len; i++) {
      alignstr->str[i] = -1;
      if (struc[i] != INFINITE_I) {
        if (struc[i] == -1) {
          if (singles_forced) alignstr->str[i] = SINGLE;
        } else if (struc[i] >= i + grammar->mindist || struc[i] + grammar->mindist <= i) {
          alignstr->str[i] = struc[i];
        }
      }
    }
  }
  alignstr->map = (int *)pf_malloc((len + 1) * sizeof(int));
  for (int i = 0; i < len; i++) alignstr->map[i] = i + 1;
  alignstr->map[len] = len;
}

std::string &default_data_dir() {
  static std::string dir = PETFOLD_DATA_DIR;
  return dir;
}

}  // namespace

void set_default_data_dir(const std::string &dir) { default_data_dir() = dir; }

std::string grammar_file(const Context &ctx, const std::string &name) {
  namespace fs = std::filesystem;
  std::vector<std::string> dirs;
  if (!ctx.opt.grammar_dir.empty()) dirs.push_back(ctx.opt.grammar_dir);
  if (const char *env = std::getenv("PETFOLDBIN")) dirs.push_back(env);
  dirs.push_back(".");  // 2.2 default: current working directory
  if (!default_data_dir().empty()) dirs.push_back(default_data_dir());
  for (const auto &d : dirs) {
    fs::path p = fs::path(d) / name;
    std::error_code ec;
    if (fs::exists(p, ec)) return p.string();
  }
  /* not found: report the first candidate (as 2.2 did) */
  return (fs::path(dirs.front()) / name).string();
}

std::vector<int> get_prob_tree(Context &ctx, const petfold::Alignment &align_gf, petfold::Matrix &paired_tree,
                               std::vector<double> &single_tree, std::string &tree, const PartStruc *constraint) {
  const int len = align_gf.len();
  paired_tree = petfold::Matrix(len, 0.0);
  single_tree.assign(len, 0.0);
  if (!ctx.opt.use_evo) return {};

  const double logp_limit = -3, robust = 0.01;  // Pfold's default settings
  std::vector<int> struc;
  ArenaScope arena;

  Grammar *grammar = read_grammar(ctx, "scfg.rate");
  Entry **entry_list = Aln2EntryList(align_gf);
  Align *align = Col2Align(entry_list, find_sym_cb, grammar);

  Phyl *tmpphyl, *phyl;
  Entry *entry;
  if (tree.empty()) {
    /* neighbour joining tree (adapted from Pfold:findphyl.c) */
    Dist *dist = AlignDist(align, grammar);
    tmpphyl = Neighbour(dist);
    FixPhyl(tmpphyl);
    entry = PhylEntry(tmpphyl, "tree");
    phyl = ReadColPhyl(entry);
    FreeDist(dist);
  } else {
    std::string tmptree = tree + ":0";
    tmpphyl = ReadPhyl(tmptree.data(), 1);
    entry = PhylEntry(tmpphyl, "tree");
    phyl = ReadColPhyl(entry);
  }
  FreePhyl(tmpphyl);
  if (ctx.opt.verbose) ctx.log.print("PFOLD:\n NJ tree = " + phyl_to_newick(phyl) + "\n");

  /* branch lengths by maximum likelihood unless given (adapted from Pfold:mltree.c) */
  if (tree.empty() || tree.find(':') == std::string::npos) {
    initphyl(phyl->root, grammar, align);
    std::vector<double> collikelhd(align->len, 0.0);
    double oldnll = -1, nll = -1;
    while (oldnll == -1 || oldnll - nll > 0.1) {
      oldnll = nll;
      setmat(phyl->root, grammar);
      downcalc(phyl->root, align);
      upcalc(phyl->root, align);
      if (ctx.opt.verbose) ctx.log.warn(" Optimize tree:");
      nll = optimize(ctx, phyl->root, align, grammar, collikelhd.data());
      if (ctx.opt.verbose) ctx.log.warn(strf(" NLL = %f\n", nll));
    }
    FixPhyl(phyl);
    UpdatePhylEntry(phyl, entry);
    phyl = ReadColPhyl(entry);
    tree = phyl_to_newick(phyl);  // store the tree (2.2: Phyl2Char)
  }
  if (ctx.opt.verbose) ctx.log.print(" ML tree = " + phyl_to_newick(phyl) + "\n");

  /* structure reliabilities (adapted from Pfold:scfg.c) */
  grammar = read_grammar(ctx, "article.grm");
  grammar->mindist = 4;
  RobustGrammar(grammar, robust);

  Alignstr *alignstr = (Alignstr *)pf_malloc(sizeof(Alignstr));
  alignstr->align = align;
  init_alignstr_str(alignstr, entry_list, constraint, len, grammar, true);
  const int wlen = alignstr->align->len;

  Inside *e = MakeInside(grammar, wlen);
  Outside *f = MakeOutside(grammar, wlen);
  Colprob *cp = MakeColprob(grammar, wlen);
  Pairprob *pp = MakePairprob(wlen);
  PpCyk *ppcyk = MakePpCyk(wlen);

  InitCol(phyl, grammar, wlen, alignstr, logp_limit);
  InitColprob(cp, alignstr);
  InitInside(e, cp);
  InitOutside(e, f, cp);
  InitPairprob(pp, e, f, cp);
  InitPpCyk(ppcyk, pp);
  const int *s = PpCykStr(ppcyk);
  struc.assign(s, s + len);

  for (int i = 0; i < len; i++) {
    for (int j = 0; j < len; j++) {
      const int w = (i > j) ? i - j : j - i;
      const int p = (i > j) ? j : i;
      paired_tree(i, j) = Edbl2Dbl(pp->dbl[p][w]);
    }
    single_tree[i] = Edbl2Dbl(pp->sgl[i]);
  }
  return struc;  // arena releases all engine memory here
}

std::vector<int> get_ppfold(int len, petfold::Matrix &paired_tree, std::vector<double> &single_tree, std::FILE *in) {
  /* read STDIN: java -jar PPfold_forPETfold.jar aln.fasta --onlyCT -t tree.newick | PETfold -f aln.fasta --ppfold */
  const std::size_t nr1 = 16u * len * (len - 1) / 2 + 16u * len + 2;
  const std::size_t nr2 = 8u * len + 2;
  std::string s1(nr1, '\0'), s2(nr2, '\0');
  if (std::fgets(s1.data(), static_cast<int>(nr1), in) == nullptr) throw Error("PPfold input stream is not readable!");
  if (std::fgets(s2.data(), static_cast<int>(nr2), in) == nullptr) throw Error("PPfold input stream is not readable!");
  paired_tree = petfold::Matrix(len, 0.0);
  single_tree.assign(len, 0.0);
  std::size_t count = 0;
  auto field = [](const std::string &s, std::size_t pos, std::size_t w) {
    return pos < s.size() ? s.substr(pos, w) : std::string();
  };
  for (int i = 0; i < len; i++)
    for (int j = 1; j < len - i; j++) {
      const double v = std::atof(field(s1, count, 16).c_str());
      paired_tree(i, j + i) = v;
      paired_tree(j + i, i) = v;
      count += 16;
    }
  for (int i = 0; i < len; i++) {
    single_tree[i] = std::atof(field(s1, count, 16).c_str());
    count += 16;
  }
  std::vector<int> str(len);
  count = 0;
  for (int i = 0; i < len; i++) {
    str[i] = std::atoi(field(s2, count, 8).c_str());
    count += 8;
  }
  return str;
}

double get_fold_evol_partial_prob(double evolconstraint, double evolensemble) {
  return std::exp(evolconstraint - evolensemble);
}

namespace {
/* log-likelihood of the alignment under the SCFG and a fixed tree, optionally
   constrained by a partial structure (from __get_evolprob in 2.2) */
double get_evolprob(Context &ctx, const petfold::Alignment &align, const std::string &tree, const PartStruc &evocon) {
  const double logp_limit = -3, robust = 0.01;
  ArenaScope arena;
  Grammar *grammar = read_grammar(ctx, "scfg.rate");
  Entry **entry_list = Aln2EntryList(align);
  Align *alignment = Col2Align(entry_list, find_sym_cb, grammar);

  std::string tmptree = tree + ":0";
  Phyl *phyl = ReadPhyl(tmptree.data(), 1);
  Entry *entry = PhylEntry(phyl, "tree");
  phyl = ReadColPhyl(entry);

  grammar = read_grammar(ctx, "article.grm");
  grammar->mindist = 4;
  RobustGrammar(grammar, robust);

  Alignstr *alignstr = (Alignstr *)pf_malloc(sizeof(Alignstr));
  alignstr->align = alignment;
  /* 2.2 did not constrain unpaired bases here (only pairs) */
  init_alignstr_str(alignstr, entry_list, &evocon, align.len(), grammar, false);
  const int wlen = alignstr->align->len;

  Inside *e = MakeInside(grammar, wlen);
  Outside *f = MakeOutside(grammar, wlen);
  Colprob *cp = MakeColprob(grammar, wlen);
  Pairprob *pp = MakePairprob(wlen);
  InitCol(phyl, grammar, wlen, alignstr, logp_limit);
  InitColprob(cp, alignstr);
  InitInside(e, cp);
  InitOutside(e, f, cp);
  const Edouble totprob = InitPairprob(pp, e, f, cp);
  return Edbl2Dbl(LogEdouble(totprob));
}
}  // namespace

double get_constraint_evolprob(Context &ctx, const petfold::Alignment &align, const std::string &tree,
                               const PartStruc &evocon) {
  return get_evolprob(ctx, align, tree, evocon);
}

double get_ensemble_evolprob(Context &ctx, const petfold::Alignment &align, const std::string &tree) {
  return get_evolprob(ctx, align, tree, PartStruc{});
}

}  // namespace petfold
