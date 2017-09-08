/*
 * petcofoldlibs.h
 *
 *  Created on: 29.12.2012
 *      Author: Stefan Seemann
 */

#ifndef PETCOFOLDLIBS_H_
#define PETCOFOLDLIBS_H_

#ifdef MODUL0
    #define EXTERN
  #else
    #define EXTERN extern
  #endif

#include "petfoldlibs.h"

EXTERN float partprob;    /* lowest accepted probability of single partial structure */
EXTERN float gap;         /* maximal allowed percent of gaps in an alignment column; deletion of column if rate is equal or higher */
EXTERN float evocon_bp;   /* minimal base pair reliability for evolutionary constraints */
EXTERN float evocon_ss;   /* minimal single stranded reliability for evolutionary constraints */
EXTERN int extstem_flag; /* constrained stems get extended by base pairs if the average reliability of the extended stem is larger than PETCON */

int *intersect_alignments(Aln *align1, Aln *align2);
Aln *get_subalignment(Aln *align, int *cindex);
PartStruc *run_petfold(Aln *align, double **paired_pet, double *single_pet, char *settree, PartStruc *setstruc);
PartStruc *get_probable_partial_struc(Aln *align, double **paired_pet, PartStruc *struc_pet, char *tree, float petcon, double *partialprob);
int *get_duplex_prob_tree(Aln *align_gf1, Aln *align_gf2, double **duplex_paired_tree, double *duplex_single_tree, PartStruc *ss_partial1, PartStruc *ss_partial2, char *settree);
void get_duplex_prob_energy(Aln *align_gf1, Aln *align_gf2, double **duplex_paired_energy, double *duplex_single_energy, PartStruc *ss_partial1, PartStruc *ss_partial2);
void adjust_duplex_paired(double **duplex_paired, int len, double partialprob);
void adjust_duplex_single(double *duplex_single, int len, double partialprob);
void merge_single_fold_duplex_pet(double *duplex_single_pet, PartStruc *ss_partial1, PartStruc *ss_partial2, int len1);
void merge_paired_fold_duplex_pet(double **duplex_paired_pet, double *partialprob1, double *partialprob2, PartStruc *ss_partial1, PartStruc *ss_partial2, int len1);
DynEntry ***nussinov_duplex(double **duplex_paired_pet, double *duplex_single_pet, PartStruc *pstruc, int len1, int len2);
void remove_linker(char *pet_db_gf, char *pet_db_gf1, char *pet_db_gf2, int len1, int len2);
double get_score(char *pet_db_gf1, char *pet_db_gf2, double **duplex_paired_pet, double *duplex_single_pet, int len1, int len2);
PartStruc *extend_partial_stems(PartStruc *partstruc, int *inputstruc, double **paired_pet, int len, float petcon);

int sort(const void *x, const void *y);
double geometric_mean(double val1, double val2);
char *substitute_char(char *string, char a, char b, int len);

void usage_petcofold();

#endif /* PETCOFOLDLIBS_H_ */
