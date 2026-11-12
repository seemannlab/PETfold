/*
 * structure.cpp - (partial) structure representations (from petfoldlibs.c).
 */
#include "core.hpp"

namespace petfold {

std::string get_dot_bracket(const std::vector<int> &s) {
  std::string db(s.size(), '.');
  for (std::size_t i = 0; i < s.size(); i++)
    if (s[i] != -1) db[i] = (s[i] > static_cast<int>(i)) ? '(' : ')';
  return db;
}

PartStruc get_partial_struc(const Matrix &paired, const std::vector<double> &single, const std::vector<int> &struc,
                            int len, float evocon_bp, float evocon_ss) {
  PartStruc p;
  for (int i = 0; i < len; i++) {
    if (struc[i] == -1) {
      if (static_cast<float>(single[i]) > evocon_ss) p.ss.push_back(i);
    } else if (struc[i] < i || struc[i] >= len) {
      continue;  // right base of a pair, or unconstrained entry (was out-of-bounds in 2.2)
    } else if (static_cast<float>(paired(i, struc[i])) > evocon_bp) {
      p.bp_left.push_back(i);
      p.bp_right.push_back(struc[i]);
    }
  }
  return p;
}

std::string get_constraint_string(const PartStruc *evocon, const std::vector<int> &origid, int len_gf) {
  std::string db(len_gf, '.');
  /* FIX: 2.2 ignored the constraints whenever --noevo was given (the test was
     on the global flag); now a constraint set is applied whenever present. */
  if (!evocon) return db;
  int k = 0;
  for (int i = 0; i < evocon->bp_nr(); i++) {
    for (int l = k; l < len_gf; l++) {
      if (origid[l] == evocon->bp_left[i]) {
        for (int m = l + 1; m < len_gf; m++) {
          if (origid[m] == evocon->bp_right[i]) {
            db[l] = '(';
            db[m] = ')';
            k = l + 1;
            break;
          }
        }
      } else if (origid[l] > evocon->bp_left[i]) {
        break;
      }
    }
  }
  k = 0;
  for (int i = 0; i < evocon->ss_nr(); i++)
    for (int l = k; l < len_gf; l++) {
      if (origid[l] == evocon->ss[i]) {
        db[l] = 'x';
        k = l + 1;
        break;
      } else if (origid[l] > evocon->ss[i]) {
        break;
      }
    }
  return db;
}

PartStruc mod_setstruct(std::string setstruc, const std::vector<int> &gap_col, int len, int gapnr) {
  if (static_cast<int>(setstruc.size()) != len)
    throw Error(strf("Error <setstruc>: Structure length (%i) is different from alignment (%i).",
                     static_cast<int>(setstruc.size()), len));
  for (char &c : setstruc) {
    if (c == '<' || c == '[') c = '(';
    else if (c == '>' || c == ']') c = ')';
    else if (c == '-') c = '.';
  }
  /* remove gap columns; a gap column must be unpaired */
  std::string tmp;
  int offset = 0;
  for (int i = 0; i < len; i++) {
    /* FIX: 2.2 read gap_col[offset] past the last gap column (uninitialised). */
    if (offset < gapnr && gap_col[offset] == i) {
      if (setstruc[i] != '.')
        throw Error(strf("Error <setstruc>: Position %i must be unpaired due to an alignment gap.", i + 1));
      offset++;
    } else {
      tmp += setstruc[i];
    }
  }
  PartStruc p;
  std::vector<int> stack;
  for (int i = 0; i < static_cast<int>(tmp.size()); i++) {
    switch (tmp[i]) {
      case '.':
        p.ss.push_back(i);
        break;
      case '(':
        stack.push_back(p.bp_nr());
        p.bp_left.push_back(i);
        p.bp_right.push_back(-1);
        break;
      case ')':
        /* FIX: unbalanced ')' overflowed the stack in 2.2 */
        if (stack.empty()) throw Error("Error <setstruc>: Number of left and right parenthesis are unequal!");
        p.bp_right[stack.back()] = i;
        stack.pop_back();
        break;
      default:
        throw Error("Error <setstruc>: The only allowed characters are '(', ')', '<', '>', '[', ']', '.', and '-'.");
    }
  }
  if (!stack.empty()) throw Error("Error <setstruc>: Number of left and right parenthesis are unequal!");
  return p;
}

std::vector<int> partstruc2list(const PartStruc *evocon, int len) {
  std::vector<int> c(len, INFINITE_I);
  if (evocon) {
    for (int i : evocon->ss) c[i] = -1;
    for (int i = 0; i < evocon->bp_nr(); i++) {
      c[evocon->bp_left[i]] = i;
      c[evocon->bp_right[i]] = i;
    }
  }
  return c;
}

std::vector<int> partstruc2str(const PartStruc *evocon, int len) {
  std::vector<int> c(len, INFINITE_I);
  if (evocon) {
    for (int i : evocon->ss) c[i] = -1;
    for (int i = 0; i < evocon->bp_nr(); i++) {
      c[evocon->bp_left[i]] = evocon->bp_right[i];
      c[evocon->bp_right[i]] = evocon->bp_left[i];
    }
  }
  return c;
}

std::string partstruc2string(const PartStruc *evocon, int len) {
  std::string db(len, '.');
  if (evocon) {
    for (int i = 0; i < evocon->bp_nr(); i++) {
      db[evocon->bp_left[i]] = '(';
      db[evocon->bp_right[i]] = ')';
    }
    for (int i : evocon->ss) db[i] = 'x';
  }
  return db;
}

PartStruc str2partstruc(const std::vector<int> &struc, int len) {
  PartStruc p;
  for (int i = 0; i < len; i++) {
    if (struc[i] == -1)
      p.ss.push_back(i);
    else if (struc[i] != INFINITE_I && i < struc[i]) {
      p.bp_left.push_back(i);
      p.bp_right.push_back(struc[i]);
    }
  }
  return p;
}

}  // namespace petfold
