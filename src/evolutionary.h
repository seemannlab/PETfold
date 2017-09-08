/*
 * evolutionary.h
 *
 *  Created on: 06.02.2011
 *      Author: Stefan Seemann
 */

#ifndef EVOLUTIONARY_H_
#define EVOLUTIONARY_H_

#include "petfoldlibs.h"
#include "col.h"
#include "phyl.h"

#define BUFFERSIZE  10000


int *get_prob_tree(Aln * align, double **paired_tree, double *single_tree, char *settree);
int *get_prob_tree_constraint(Aln* align, double **paired_tree, double *single_tree, PartStruc *partstruc, char *settree);
int *__get_pfold(Aln *align_gf, double **paired_tree, double *single_tree, PartStruc *partstruc, char *settree);
//int *get_prob_tree(Aln *align_gf, double **paired_tree, double *single_tree, char *settree);
int *get_ppfold(int len, double **paired_tree, double *single_tree);
int Aln2EntryList(Aln *align, Entry **entry_list);
Entry *MyNewEntry(char *type, char *name, int len);

/* copied from Pfold:mltree.c */
typedef struct NodeInfo {
  int seqnum;               /* Sequence at this node */
  Matrix *pmatrix;          /* Matrix for up branch evolution */
  Matrix **p_top_up;         /* Probability distribution of tree above
			       top of up branch */
  Matrix **p_top_down;       /* Probability distribution of tree below
			       top of up branch */
  Matrix **p_bot_up;         /* Probability distribution of tree above
			       bottom of up branch */
  Matrix **p_bot_down;       /* Probability distribution of tree below
			       bottom of up branch */
  Grammar *grammar;
} NodeInfo;

void initphyl(PhylNode *pnode, Grammar *grammar, Align *align);
void setmat(PhylNode *pnode, Grammar *grammar);
void downcalc(PhylNode *pnode, Align *align);
void upcalc(PhylNode *pnode, Align *align);
double optimize(PhylNode *pnode, Align *align, Grammar *grammar, Phyl *phyl, double *collikelhd);

double get_fold_evol_partial_prob(double evolconstraint, double evolensemble);
double get_constraint_evolprob(Aln * align, char *tree, PartStruc *evocon);
double get_ensemble_evolprob(Aln* align, char *tree);
double __get_evolprob(Aln* align, char *tree, PartStruc *evocon);

void Phyl2Char(char *tree, Phyl *phyl);
int __Phyl2Char(char *tree, int offset, PhylNode *pnode);

#endif /* EVOLUTIONARY_H_ */
