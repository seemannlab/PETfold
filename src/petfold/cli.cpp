/*
 * cli.cpp - command-line front ends of PETfold and PETcofold.
 *
 * Options and output are identical to PETfold 2.2 / PETcofold 3.4; new
 * options are listed at the end of the help texts.
 */
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "../pfold/pf_alloc.hpp"
#include "api.hpp"
#include "core.hpp"

namespace petfold {

namespace {

/* ------------------------------------------------------------------ */
/* minimal portable replacement for getopt_long (supports -x v, -xv,    */
/* --long v, --long=v, unambiguous prefixes of long options)           */

struct LongOpt {
  const char *name;
  bool has_arg;
  int key;          // returned key (0 if flag)
  int *flag;        // set to `val` if not null
  int val;
};

struct ParsedOpt {
  int key;
  std::string arg;
};

class UsageError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

std::vector<ParsedOpt> parse_options(int argc, char **argv, const char *shortopts, const std::vector<LongOpt> &longopts) {
  std::vector<ParsedOpt> out;
  for (int a = 1; a < argc; a++) {
    std::string s = argv[a];
    if (s == "--") break;
    if (s.rfind("--", 0) == 0) {
      std::string name = s.substr(2), value;
      bool has_value = false;
      if (auto eq = name.find('='); eq != std::string::npos) {
        value = name.substr(eq + 1);
        name = name.substr(0, eq);
        has_value = true;
      }
      const LongOpt *match = nullptr;
      int nmatch = 0;
      for (const auto &lo : longopts) {
        if (name == lo.name) {
          match = &lo;
          nmatch = 1;
          break;
        }
        if (std::strncmp(lo.name, name.c_str(), name.size()) == 0) {
          match = &lo;
          nmatch++;
        }
      }
      if (nmatch != 1 || name.empty())
        throw UsageError(std::string(argv[0]) + ": " + (nmatch > 1 ? "option is ambiguous" : "unrecognized option") +
                         " '" + s + "'");
      if (match->has_arg) {
        if (!has_value) {
          if (a + 1 >= argc) throw UsageError(std::string(argv[0]) + ": option '" + s + "' requires an argument");
          value = argv[++a];
        }
      } else if (has_value) {
        throw UsageError(std::string(argv[0]) + ": option '--" + name + "' doesn't allow an argument");
      }
      if (match->flag) *match->flag = match->val;
      else out.push_back({match->key, value});
    } else if (s.size() > 1 && s[0] == '-') {
      for (std::size_t p = 1; p < s.size(); p++) {
        const char c = s[p];
        const char *q = std::strchr(shortopts, c);
        if (!q || c == ':') throw UsageError(std::string(argv[0]) + ": invalid option -- '" + c + "'");
        if (q[1] == ':') {
          std::string value = s.substr(p + 1);
          if (value.empty()) {
            if (a + 1 >= argc)
              throw UsageError(std::string(argv[0]) + ": option requires an argument -- '" + c + "'");
            value = argv[++a];
          }
          out.push_back({c, value});
          break;
        }
        out.push_back({c, ""});
      }
    }
    /* non-option arguments are ignored (as in 2.2) */
  }
  return out;
}

Log stdio_log() {
  Log l;
  l.out = [](const std::string &s) { std::fputs(s.c_str(), stdout); };
  l.err = [](const std::string &s) { std::fputs(s.c_str(), stderr); };
  return l;
}

int report(const std::exception &e) {
  std::fflush(stdout);
  if (auto *pe = dynamic_cast<const Error *>(&e); pe && pe->to_stdout()) {
    std::printf("%s\n", e.what());
    std::fflush(stdout);
  } else {
    std::fprintf(stderr, "%s\n", e.what());
  }
  return 1;
}

/* ------------------------------------------------------------------ */

void usage_petfold() {
  std::printf(
      "PETfold v2.2 (C++ edition %s)\n"
      "============\n"
      "by Stefan E Seemann (seemann@rth.dk)\n"
      "Reference: Seemann et al. Nucleic Acids Res. 36(20):6355-62, 2008\n"
      "Web service: http://rth.dk/resources/petfold\n"
      "\n"
      "Usage:\n"
      "  PETfold -f <file> [ options ] [ parameter settings ]\n"
      "\n"
      "  -f --fasta <file>          ... alignment in fasta format\n"
      "Options:\n"
      "  -s --setstruc <struc|file> ... calculates score for given structure in dot-bracket notation\n"
      "  -t --settree <tree|file>   ... calculates score for given tree in Newick tree format\n"
      "  --war                      ... fasta format output\n"
      "  -r --ppfile <file>         ... writes PET reliabilities in file\n"
      "  -o --suboptimal <nr>       ... number of alternative structures found by sampling\n"
      "  --ppfold                   ... reads ppfold calculated reliabilities from STDIN\n"
      "  --verbose                  ... writes long output\n"
      "  --help                     ... this output\n"
      "Parameter settings:\n"
      "  -p --setevocon_bp <reliab> ... reliab.threshold for conserved base pairs (default: 0.9)\n"
      "  -u --setevocon_ss <reliab> ... reliab.threshold for conserved unpaired bases (default: 1)\n"
      "  -a --setalpha <nr>         ... weighting factor for unpaired reliabilities (default: 0.2)\n"
      "  -b --setbeta <nr>          ... weighting factor for thermodynamic overlap (default: 1)\n"
      "  -g --setgap <nr>           ... max. percent of gaps in alignment column (default:0.25)\n"
      "New in the C++ edition:\n"
      "  --seed <nr>                ... random seed for --suboptimal (default: time)\n"
      "  --grammar-dir <dir>        ... directory with scfg.rate and article.grm\n"
      "  --noevo / --nothermo       ... switch off the evolutionary / thermodynamic model\n"
      "  --version                  ... print version\n"
      "\n"
      "You may set the environment variable PETFOLDBIN to path of files scfg.rate and article.grm:\n"
      "export PETFOLDBIN=<PATH>\n"
      "\n"
      "Parse sampling of suboptimal structures:\n"
      "PETfold -f <file> -o 10000 | grep Subopt | sort -k 4 -ur | awk '{print $3,$4}' | less\n\n",
      kVersion);
}

void usage_petcofold() {
  std::printf(
      "PETcofold v3.4 (C++ edition %s)\n"
      "==============\n"
      "by Stefan E Seemann (seemann@rth.dk)\n"
      "Reference: Seemann, Richter et al. Algorithms Mol Biol. 5:22, 2010\n"
      "           Seemann, Richter et al. Bioinformatics 27(2):211-9, 2011\n"
      "Web service: http://rth.dk/resources/petcofold\n"
      "\n"
      "Usage:\n"
      "  PETcofold -f <file1> -f <file2> [ options ] [ parameter settings ]\n"
      "\n"
      "  -f --fasta <file1>         ... 1st alignment in fasta format\n"
      "  -f --fasta <file2>         ... 2nd alignment in fasta format with same organisms than 1st\n"
      "Options:\n"
      "  -s --setstruc <struc|file> ... calculates score for given duplex structure in dot-bracket notation\n"
      "                             ... step 1 constraints intramolecular ('(',')') and step 2 intermolecular structure "
      "('[',']')\n"
      "  -n --setstruc1 <struc|file>... finds reliable base pairs in given structure of 1st alignment in dot-bracket "
      "notation\n"
      "  -m --setstruc2 <struc|file>... finds reliable base pairs in given structure of 2nd alignment in dot-bracket "
      "notation\n"
      "  -t --settree <tree|file>   ... calculates score for given tree in Newick tree format\n"
      "  --war                      ... fasta format output\n"
      "  --intermol                 ... structure output of intermolecular base pairs\n"
      "  -r --ppfile <file>         ... writes PET reliabilities in file\n"
      "  --verbose                  ... writes long output\n"
      "  --help                     ... this output\n"
      "Parameter settings:\n"
      "  -p --setevocon_bp <reliab> ... reliab.threshold for conserved base pairs (default: 0.9)\n"
      "  -u --setevocon_ss <reliab> ... reliab.threshold for conserved unpaired bases (default: 1)\n"
      "  -a --setalpha <nr>         ... weighting factor for unpaired reliabilities (default: 0.2)\n"
      "  -b --setbeta <nr>          ... weighting factor for thermodynamic overlap (default: 1)\n"
      "  -g --setgap <nr>           ... max. percent of gaps in alignment column (default:0.25)\n"
      "  -d --setdelta <reliab>     ... max. intramolecular base pair reliability\n"
      "                                 to be free for RNA-RNA interaction (former petcon; default: 0.9)\n"
      "  -i --setgamma <probab>     ... minimal partial structure probability (former partprob; default: 0.1)\n"
      "  --extstem                  ... constrained stems get extended by inner and outer base pairs\n"
      "New in the C++ edition:\n"
      "  --legacy-cutpoint          ... RNAcofold cut point as in PETcofold 3.4 (ignores gaps; for comparisons)\n"
      "  --grammar-dir <dir>        ... directory with scfg.rate and article.grm\n"
      "  --noevo / --nothermo       ... switch off the evolutionary / thermodynamic model\n"
      "  --version                  ... print version\n"
      "\n"
      "You may set the environment variable PETFOLDBIN to path of files scfg.rate and article.grm:\n"
      "export PETFOLDBIN=<PATH>\n\n",
      kVersion);
}

float to_float(const std::string &s) { return static_cast<float>(std::atof(s.c_str())); }

}  // namespace

int petfold_main(int argc, char **argv) {
  try {
    Context ctx;
    ctx.log = stdio_log();
    int setevo = 1, setthermo = 1, war = 0, help = 0, verbose = 0, ppfold = 0, version = 0;
    std::string fasta, settree, setstruc, ppfile, seed;
    int suboptimal = 0;
    enum { K_SEED = 1000, K_GRAMMAR };
    const std::vector<LongOpt> longopts = {
        {"setevo", false, 0, &setevo, 1},       {"noevo", false, 0, &setevo, 0},
        {"setthermo", false, 0, &setthermo, 1}, {"nothermo", false, 0, &setthermo, 0},
        {"war", false, 0, &war, 1},             {"help", false, 0, &help, 1},
        {"verbose", false, 0, &verbose, 1},     {"settree", true, 't', nullptr, 0},
        {"fasta", true, 'f', nullptr, 0},       {"setgap", true, 'g', nullptr, 0},
        {"setevocon_bp", true, 'p', nullptr, 0}, {"setevocon_ss", true, 'u', nullptr, 0},
        {"setbeta", true, 'b', nullptr, 0},     {"setstruc", true, 's', nullptr, 0},
        {"setalpha", true, 'a', nullptr, 0},    {"ppfile", true, 'r', nullptr, 0},
        {"suboptimal", true, 'o', nullptr, 0},  {"ppfold", false, 0, &ppfold, 1},
        {"seed", true, K_SEED, nullptr, 0},     {"grammar-dir", true, K_GRAMMAR, nullptr, 0},
        {"version", false, 0, &version, 1},
    };
    std::vector<ParsedOpt> opts;
    try {
      opts = parse_options(argc, argv, "f:g:t:p:u:b:s:a:r:o:", longopts);
    } catch (const UsageError &e) {
      std::fprintf(stderr, "%s\n", e.what());
      usage_petfold();
      return 1;
    }
    for (const auto &o : opts) {
      switch (o.key) {
        case 'f': fasta = o.arg; break;
        case 'g': ctx.opt.gap = to_float(o.arg); break;
        case 't': settree = o.arg; break;
        case 'p': ctx.opt.evocon_bp = to_float(o.arg); break;
        case 'u': ctx.opt.evocon_ss = to_float(o.arg); break;
        case 'b': ctx.opt.beta = to_float(o.arg); break;
        case 's': setstruc = o.arg; break;
        case 'a': ctx.opt.alpha = to_float(o.arg); break;
        case 'r': ppfile = o.arg; break;
        case 'o': suboptimal = std::atoi(o.arg.c_str()); break;
        case K_SEED: seed = o.arg; break;
        case K_GRAMMAR: ctx.opt.grammar_dir = o.arg; break;
      }
    }
    if (version) {
      std::printf("PETfold %s (C++ edition of PETfold 2.2)\n", kVersion);
      return 0;
    }
    if (help || fasta.empty()) {
      usage_petfold();
      return 1;
    }
    ctx.opt.use_evo = setevo;
    ctx.opt.use_thermo = setthermo;
    ctx.opt.verbose = verbose;
    ctx.rng.seed(seed.empty() ? static_cast<unsigned long long>(
                                    std::chrono::high_resolution_clock::now().time_since_epoch().count())
                              : std::strtoull(seed.c_str(), nullptr, 10));

    PetfoldInput in;
    in.alignment = read_fasta(fasta);
    if (!settree.empty()) in.tree = read_file_or_string(settree);
    if (!setstruc.empty()) in.structure = read_file_or_string(setstruc);
    in.suboptimal = suboptimal;
    if (ppfold) in.ppfold_input = stdin;

    const PetfoldResult r = run_petfold(ctx, in);
    const int gapnr = r.alignment_length - r.alignment_gf.len();

    if (!ppfile.empty()) {
      try {
        create_ppfile(r.paired_pet, r.single_pet, r.gap_columns, r.alignment_gf.len(), gapnr, ppfile);
      } catch (const Error &e) {
        std::fprintf(stderr, "%s\n", e.what());  // 2.2 reported the error and continued
      }
    }

    if (war) {
      for (int i = 0; i < in.alignment.nr(); i++)
        std::printf(">%s\n%s\n", in.alignment.ids[i].c_str(), in.alignment.raw[i].c_str());
      std::printf(">structure\n%s\n", r.structure.c_str());
    } else {
      if (r.has_evo) {
        std::printf("Pfold RNA structure:\t%s\n",
                    include_gap_columns(r.pfold_structure_gf, r.gap_columns, r.alignment_length).c_str());
        std::printf("Constraints:\t\t%s\n",
                    include_gap_columns(r.constraints_gf, r.gap_columns, r.alignment_length).c_str());
      }
      if (!r.structure.empty()) std::printf("PETfold RNA structure:\t%s\n", r.structure.c_str());
      std::printf("Score_{model,structure}{tree,alignment} = %7.6f\n", r.score);
      std::printf("Length-normalized ensemble diversity = %7.6f\n", r.ensemble_diversity);
    }
    for (const auto &[s, sc] : r.suboptimals) std::printf("Suboptimal structure:   %s\t%7.6f\n", s.c_str(), sc);
    std::fflush(stdout);
    return 0;
  } catch (const std::exception &e) {
    return report(e);
  }
}

int petcofold_main(int argc, char **argv) {
  try {
    Context ctx;
    ctx.log = stdio_log();
    int setevo = 1, setthermo = 1, war = 0, help = 0, verbose = 0, intermol = 0, extstem = 0, legacy = 0,
        version = 0;
    std::string fasta1, fasta2, settree, setstruc, setstruc1, setstruc2, ppfile;
    PetcofoldInput in;
    enum { K_GRAMMAR = 1000 };
    const std::vector<LongOpt> longopts = {
        {"setevo", false, 0, &setevo, 1},        {"noevo", false, 0, &setevo, 0},
        {"setthermo", false, 0, &setthermo, 1},  {"nothermo", false, 0, &setthermo, 0},
        {"war", false, 0, &war, 1},              {"help", false, 0, &help, 1},
        {"verbose", false, 0, &verbose, 1},      {"intermol", false, 0, &intermol, 1},
        {"settree", true, 't', nullptr, 0},      {"fasta", true, 'f', nullptr, 0},
        {"setgap", true, 'g', nullptr, 0},       {"setevocon_bp", true, 'p', nullptr, 0},
        {"setevocon_ss", true, 'u', nullptr, 0}, {"setbeta", true, 'b', nullptr, 0},
        {"setstruc", true, 's', nullptr, 0},     {"setstruc1", true, 'n', nullptr, 0},
        {"setstruc2", true, 'm', nullptr, 0},    {"setalpha", true, 'a', nullptr, 0},
        {"setdelta", true, 'd', nullptr, 0},     {"setgamma", true, 'i', nullptr, 0},
        {"extstem", false, 0, &extstem, 1},      {"ppfile", true, 'r', nullptr, 0},
        {"legacy-cutpoint", false, 0, &legacy, 1}, {"grammar-dir", true, K_GRAMMAR, nullptr, 0},
        {"version", false, 0, &version, 1},
    };
    std::vector<ParsedOpt> opts;
    try {
      opts = parse_options(argc, argv, "f:g:t:p:u:b:s:n:m:a:d:i:r:", longopts);
    } catch (const UsageError &e) {
      std::fprintf(stderr, "%s\n", e.what());
      usage_petcofold();
      return 1;
    }
    for (const auto &o : opts) {
      switch (o.key) {
        case 'f':
          if (fasta1.empty()) fasta1 = o.arg;
          else if (fasta2.empty()) fasta2 = o.arg;
          else {
            usage_petcofold();
            return 1;
          }
          break;
        case 'g': ctx.opt.gap = to_float(o.arg); break;
        case 't': settree = o.arg; break;
        case 'p': ctx.opt.evocon_bp = to_float(o.arg); break;
        case 'u': ctx.opt.evocon_ss = to_float(o.arg); break;
        case 'b': ctx.opt.beta = to_float(o.arg); break;
        case 's': setstruc = o.arg; break;
        case 'n': setstruc1 = o.arg; break;
        case 'm': setstruc2 = o.arg; break;
        case 'a': ctx.opt.alpha = to_float(o.arg); break;
        case 'd': in.delta = to_float(o.arg); break;
        case 'i': in.gamma = to_float(o.arg); break;
        case 'r': ppfile = o.arg; break;
        case K_GRAMMAR: ctx.opt.grammar_dir = o.arg; break;
      }
    }
    if (version) {
      std::printf("PETcofold %s (C++ edition of PETcofold 3.4)\n", kVersion);
      return 0;
    }
    if (help || fasta1.empty() || fasta2.empty()) {
      usage_petcofold();
      return 1;
    }
    ctx.opt.use_evo = setevo;
    ctx.opt.use_thermo = setthermo;
    ctx.opt.verbose = verbose;
    in.extstem = extstem;
    in.legacy_cutpoint = legacy;

    in.alignment1 = read_fasta(fasta1);
    in.alignment2 = read_fasta(fasta2);
    if (!settree.empty()) in.tree = read_file_or_string(settree);
    if (!setstruc.empty()) in.structure = read_file_or_string(setstruc);
    else {
      if (!setstruc1.empty()) in.structure1 = read_file_or_string(setstruc1);
      if (!setstruc2.empty()) in.structure2 = read_file_or_string(setstruc2);
    }

    const PetcofoldResult r = run_petcofold(ctx, in, ppfile);
    if (war) {
      for (std::size_t i = 0; i < r.fasta_ids.size(); i++)
        std::printf(">%s\n%s&%s\n", r.fasta_ids[i].c_str(), r.fasta_seq1[i].c_str(), r.fasta_seq2[i].c_str());
      if (intermol)
        std::printf(">structure\n%s&%s\n", r.intermolecular1.c_str(), r.intermolecular2.c_str());
      else
        std::printf(">structure\n%s&%s\n", r.structure1.c_str(), r.structure2.c_str());
    } else {
      if (r.has_evo)
        std::printf("PETfold RNA structure:\t\t%s %s\n", r.petfold_structure1.c_str(), r.petfold_structure2.c_str());
      if (r.structure1.size() + r.structure2.size())
        std::printf("PETcofold RNA structure:\t%s&%s\n", r.structure1.c_str(), r.structure2.c_str());
      if (intermol)
        std::printf("Intermolecular RNA structure:\t%s&%s\n", r.intermolecular1.c_str(), r.intermolecular2.c_str());
      std::printf("Score_{model,structure}{tree,alignment} = %7.6f\n", r.score);
    }
    std::fflush(stdout);
    return 0;
  } catch (const std::exception &e) {
    return report(e);
  }
}

}  // namespace petfold
