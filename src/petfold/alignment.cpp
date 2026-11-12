/*
 * alignment.cpp - FASTA input and gap handling (from petfoldlibs.c).
 */
#include <cctype>
#include <cstdarg>
#include <fstream>
#include <sstream>

#include "core.hpp"

namespace petfold {

std::string strf(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  va_list ap2;
  va_copy(ap2, ap);
  int n = std::vsnprintf(nullptr, 0, fmt, ap);
  va_end(ap);
  std::string s(n > 0 ? n : 0, '\0');
  if (n > 0) std::vsnprintf(s.data(), n + 1, fmt, ap2);
  va_end(ap2);
  return s;
}

Alignment parse_fasta(const std::string &text, const std::string &source) {
  Alignment a;
  std::istringstream in(text);
  std::string line;
  bool seen_header = false;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();  // tolerate CRLF files
    if (!line.empty() && line[0] == '>') {
      std::string id = line.substr(1);
      for (char &c : id)
        if (c == ' ') c = '_';  // spaces in identifiers become underscores (as in 2.2)
      a.ids.push_back(id);
      a.raw.emplace_back();
      seen_header = true;
    } else {
      if (!seen_header) {
        if (line.find_first_not_of(" \t") == std::string::npos) continue;
        throw Error("'" + source + "' does not look like a FASTA file (sequence before first '>').");
      }
      for (char &c : line) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));  // RNAfold needs upper case
        if (c == 'T') c = 'U';
      }
      a.raw.back() += line;
    }
  }
  if (a.ids.empty()) throw Error("'" + source + "' does not contain any FASTA record.");
  a.seqs = a.raw;
  /* FIX: '.' was counted as a gap when removing columns but treated as a
     nucleotide everywhere else; gaps are now normalised to '-'. */
  for (auto &s : a.seqs)
    for (char &c : s)
      if (c == '.') c = '-';
  return a;
}

Alignment read_fasta(const std::string &path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw Error("Can't open the file '" + path + "'");
  std::stringstream ss;
  ss << f.rdbuf();
  return parse_fasta(ss.str(), path);
}

std::string read_file_or_string(const std::string &arg) {
  std::ifstream f(arg, std::ios::binary);
  if (f) {
    std::string line;
    if (std::getline(f, line)) {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      return line;
    }
  }
  return arg;
}

Alignment delete_gap_columns(const Context &ctx, const Alignment &align, std::vector<int> &gap_col, float gap) {
  gap_col.clear();
  std::vector<int> k;
  const int len = align.len(), nr = align.nr();
  for (int j = 0; j < len; j++) {
    int n = 0;
    for (int i = 0; i < nr; i++)
      if (align.seqs[i][j] == '-') n++;
    if (static_cast<double>(n) / nr >= gap)
      gap_col.push_back(j);
    else
      k.push_back(j);
  }
  Alignment gf;
  for (int i = 0; i < nr; i++) {
    std::string s;
    s.reserve(k.size());
    int ngap = 0;
    for (int col : k) {
      s += align.seqs[i][col];
      if (align.seqs[i][col] == '-') ngap++;
    }
    if (ngap == static_cast<int>(k.size())) {
      ctx.log.print(strf("Sequence %s is removed because it has only gaps in gap-column-free alignment.\n",
                         align.ids[i].c_str()));
      continue;
    }
    gf.ids.push_back(align.ids[i]);
    gf.raw.push_back(s);
    gf.seqs.push_back(std::move(s));
  }
  return gf;
}

std::string include_gap_columns(const std::string &db, const std::vector<int> &gap_col, int full_len) {
  /* FIX: the 2.2 version shifted the string in place and read before the
     start of the buffer when column 0 was a gap column. */
  std::string out;
  out.reserve(full_len);
  std::size_t g = 0, p = 0;
  for (int col = 0; col < full_len; col++) {
    if (g < gap_col.size() && gap_col[g] == col) {
      out += '-';
      g++;
    } else {
      out += (p < db.size()) ? db[p++] : '-';
    }
  }
  return out;
}

SeqList get_align_without_gaps(const Alignment &align) {
  SeqList sl;
  sl.seq.resize(align.nr());
  sl.origid.resize(align.nr());
  for (int i = 0; i < align.nr(); i++)
    for (int j = 0; j < align.len(); j++)
      if (align.seqs[i][j] != '-') {
        sl.origid[i].push_back(j);
        sl.seq[i] += align.seqs[i][j];
      }
  return sl;
}

}  // namespace petfold
