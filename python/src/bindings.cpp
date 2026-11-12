/*
 * bindings.cpp - pybind11 module petfold._core
 */
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <chrono>
#include <cstdio>
#include <mutex>

#include "petfold/api.hpp"
#include "petfold/core.hpp"
#include "apps/tools.hpp"
#include "pfold/pf_alloc.hpp"

namespace py = pybind11;
using namespace petfold;

namespace {

Alignment make_alignment(const std::vector<std::string> &ids, const std::vector<std::string> &seqs) {
  if (ids.size() != seqs.size()) throw Error("ids and sequences must have the same length");
  std::string text;
  for (std::size_t i = 0; i < ids.size(); i++) text += ">" + ids[i] + "\n" + seqs[i] + "\n";
  return parse_fasta(text, "<python>");
}

/* pair matrix in alignment coordinates: upper triangle mirrored, zeros in
   removed columns and on the diagonal */
py::array_t<double> pair_matrix(const Matrix &m, const std::vector<int> &gap_col, int full_len) {
  std::vector<int> map(full_len, -1);
  std::size_t g = 0;
  for (int c = 0, p = 0; c < full_len; c++) {
    if (g < gap_col.size() && gap_col[g] == c) g++;
    else map[c] = p++;
  }
  py::array_t<double> a({full_len, full_len});
  auto r = a.mutable_unchecked<2>();
  for (int x = 0; x < full_len; x++)
    for (int y = 0; y < full_len; y++) {
      double v = 0.0;
      if (map[x] >= 0 && map[y] >= 0 && x != y) {
        const int i = std::min(map[x], map[y]), j = std::max(map[x], map[y]);
        v = m(i, j);
        if (v == INFINITE_D) v = 0.0;
      }
      r(x, y) = v;
    }
  return a;
}

py::array_t<double> single_vector(const std::vector<double> &v, const std::vector<int> &gap_col, int full_len) {
  py::array_t<double> a(full_len);
  auto r = a.mutable_unchecked<1>();
  std::size_t g = 0, p = 0;
  for (int c = 0; c < full_len; c++) {
    if (g < gap_col.size() && gap_col[g] == c) {
      r(c) = 0.0;
      g++;
    } else {
      r(c) = p < v.size() ? v[p++] : 0.0;
    }
  }
  return a;
}

Context make_context(float gap, float evocon_bp, float evocon_ss, float alpha, float beta, bool use_evo,
                     bool use_thermo, bool verbose, const std::string &grammar_dir, py::object seed,
                     std::vector<std::string> *messages) {
  Context ctx;
  ctx.opt.gap = gap;
  ctx.opt.evocon_bp = evocon_bp;
  ctx.opt.evocon_ss = evocon_ss;
  ctx.opt.alpha = alpha;
  ctx.opt.beta = beta;
  ctx.opt.use_evo = use_evo;
  ctx.opt.use_thermo = use_thermo;
  ctx.opt.verbose = verbose;
  ctx.opt.grammar_dir = grammar_dir;
  ctx.rng.seed(seed.is_none() ? static_cast<unsigned long long>(
                                    std::chrono::high_resolution_clock::now().time_since_epoch().count())
                              : seed.cast<unsigned long long>());
  ctx.log.out = [messages](const std::string &s) { messages->push_back(s); };
  ctx.log.err = [messages](const std::string &s) { messages->push_back(s); };
  return ctx;
}

std::string join(const std::vector<std::string> &v) {
  std::string s;
  for (const auto &x : v) s += x;
  return s;
}

py::dict fold(const std::vector<std::string> &ids, const std::vector<std::string> &seqs, const std::string &tree,
              const std::string &structure, int suboptimal, py::object seed, float gap, float evocon_bp,
              float evocon_ss, float alpha, float beta, bool use_evo, bool use_thermo, bool verbose,
              const std::string &grammar_dir) {
  std::vector<std::string> messages;
  Context ctx = make_context(gap, evocon_bp, evocon_ss, alpha, beta, use_evo, use_thermo, verbose, grammar_dir, seed,
                             &messages);
  PetfoldInput in;
  in.alignment = make_alignment(ids, seqs);
  in.tree = tree;
  in.structure = structure;
  in.suboptimal = suboptimal;
  PetfoldResult r;
  {
    py::gil_scoped_release release;
    r = run_petfold(ctx, in);
  }
  const int L = r.alignment_length;
  const auto &gc = r.gap_columns;
  py::dict d;
  d["structure"] = r.structure;
  d["score"] = r.score;
  d["ensemble_diversity"] = r.ensemble_diversity;
  d["pfold_structure"] = r.has_evo ? include_gap_columns(r.pfold_structure_gf, gc, L) : std::string();
  d["constraints"] = r.has_evo ? include_gap_columns(r.constraints_gf, gc, L) : std::string();
  d["tree"] = r.tree;
  d["gap_columns"] = gc;
  d["ids"] = in.alignment.ids;
  d["sequences"] = in.alignment.raw;
  d["paired_pet"] = pair_matrix(r.paired_pet, gc, L);
  d["single_pet"] = single_vector(r.single_pet, gc, L);
  d["paired_evo"] = pair_matrix(r.paired_tree, gc, L);
  d["single_evo"] = single_vector(r.single_tree, gc, L);
  d["paired_thermo"] = pair_matrix(r.paired_seq, gc, L);
  d["single_thermo"] = single_vector(r.single_seq, gc, L);
  py::list sub;
  for (const auto &[s, sc] : r.suboptimals) sub.append(py::make_tuple(s, sc));
  d["suboptimals"] = sub;
  d["log"] = join(messages);
  return d;
}

py::dict cofold(const std::vector<std::string> &ids1, const std::vector<std::string> &seqs1,
                const std::vector<std::string> &ids2, const std::vector<std::string> &seqs2, const std::string &tree,
                const std::string &structure, const std::string &structure1, const std::string &structure2,
                float delta, float gamma, bool extstem, bool legacy_cutpoint, float gap, float evocon_bp,
                float evocon_ss, float alpha, float beta, bool use_evo, bool use_thermo, bool verbose,
                const std::string &grammar_dir) {
  std::vector<std::string> messages;
  Context ctx = make_context(gap, evocon_bp, evocon_ss, alpha, beta, use_evo, use_thermo, verbose, grammar_dir,
                             py::none(), &messages);
  PetcofoldInput in;
  in.alignment1 = make_alignment(ids1, seqs1);
  in.alignment2 = make_alignment(ids2, seqs2);
  in.tree = tree;
  in.structure = structure;
  in.structure1 = structure1;
  in.structure2 = structure2;
  in.delta = delta;
  in.gamma = gamma;
  in.extstem = extstem;
  in.legacy_cutpoint = legacy_cutpoint;
  PetcofoldResult r;
  {
    py::gil_scoped_release release;
    r = run_petcofold(ctx, in);
  }
  std::vector<int> gc = r.gap_columns1;
  for (int c : r.gap_columns2) gc.push_back(c + r.length1);
  const int L = r.length1 + r.length2;
  py::dict d;
  d["structure"] = r.structure1 + "&" + r.structure2;
  d["structure1"] = r.structure1;
  d["structure2"] = r.structure2;
  d["intermolecular"] = r.intermolecular1 + "&" + r.intermolecular2;
  d["partial_structure1"] = r.partial_structure1;
  d["partial_structure2"] = r.partial_structure2;
  d["petfold_structure1"] = r.petfold_structure1;
  d["petfold_structure2"] = r.petfold_structure2;
  d["delta1"] = r.delta1;
  d["delta2"] = r.delta2;
  d["partial_prob1"] = py::make_tuple(r.partial_prob1[0], r.partial_prob1[1]);
  d["partial_prob2"] = py::make_tuple(r.partial_prob2[0], r.partial_prob2[1]);
  d["score"] = r.score;
  d["common_ids"] = r.common_ids;
  d["gap_columns1"] = r.gap_columns1;
  d["gap_columns2"] = r.gap_columns2;
  d["paired_pet"] = pair_matrix(r.duplex_paired_pet, gc, L);
  d["single_pet"] = single_vector(r.duplex_single_pet, gc, L);
  d["log"] = join(messages);
  return d;
}

int call_main(int (*fn)(int, char **), const std::string &prog, const std::vector<std::string> &args) {
  std::vector<std::string> storage{prog};
  storage.insert(storage.end(), args.begin(), args.end());
  std::vector<char *> argv;
  for (auto &s : storage) argv.push_back(s.data());
  argv.push_back(nullptr);
  int rc;
  {
    py::gil_scoped_release release;
    rc = fn(static_cast<int>(storage.size()), argv.data());
    std::fflush(stdout);
    std::fflush(stderr);
  }
  return rc;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "C++ core of PETfold / PETcofold";
  m.attr("__version__") = kVersion;

  py::register_exception<Error>(m, "PetfoldError", PyExc_RuntimeError);
  py::register_exception<pfold::EngineError>(m, "EngineError", PyExc_RuntimeError);

  m.def("set_default_data_dir", &set_default_data_dir, py::arg("path"));
  m.def("fold", &fold, py::arg("ids"), py::arg("sequences"), py::arg("tree") = "", py::arg("structure") = "",
        py::arg("suboptimal") = 0, py::arg("seed") = py::none(), py::arg("gap") = 0.25f,
        py::arg("evocon_bp") = 0.9f, py::arg("evocon_ss") = 1.0f, py::arg("alpha") = 0.2f, py::arg("beta") = 1.0f,
        py::arg("use_evo") = true, py::arg("use_thermo") = true, py::arg("verbose") = false,
        py::arg("grammar_dir") = "");
  m.def("cofold", &cofold, py::arg("ids1"), py::arg("sequences1"), py::arg("ids2"), py::arg("sequences2"),
        py::arg("tree") = "", py::arg("structure") = "", py::arg("structure1") = "", py::arg("structure2") = "",
        py::arg("delta") = 0.9f, py::arg("gamma") = 0.1f, py::arg("extstem") = false,
        py::arg("legacy_cutpoint") = false, py::arg("gap") = 0.25f, py::arg("evocon_bp") = 0.9f,
        py::arg("evocon_ss") = 1.0f, py::arg("alpha") = 0.2f, py::arg("beta") = 1.0f, py::arg("use_evo") = true,
        py::arg("use_thermo") = true, py::arg("verbose") = false, py::arg("grammar_dir") = "");
  m.def("petfold_main", [](const std::vector<std::string> &a) { return call_main(petfold_main, "PETfold", a); });
  m.def("petcofold_main", [](const std::vector<std::string> &a) { return call_main(petcofold_main, "PETcofold", a); });
  m.def("drawdot_main", [](const std::vector<std::string> &a) { return call_main(drawdot_main, "drawdot", a); });
  m.def("drawphyl_main", [](const std::vector<std::string> &a) { return call_main(drawphyl_main, "drawphyl", a); });
}
