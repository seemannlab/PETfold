/*
 * thermodynamics.h
 *
 *  Created on: 06.02.2011
 *      Author: Stefan Seemann
 */

#ifndef THERMODYNAMIC_H_
#define THERMODYNAMIC_H_

#include "petfoldlibs.h"

double **get_prob_paired_seq(Aln *align_gf, PartStruc *evocon);
double *get_prob_unpaired_seq(Aln *align_gf, double **paired_seq);
SeqList *get_align_without_gaps(Aln *align_gf);
double get_fold_thermo_partial_prob(double energyconstraint, double energyensemble);
double get_constraint_energy(Aln * align, PartStruc *evocon);
double get_ensemble_energy(Aln* align);
void get_cofold_prob_paired_seq(Aln *align, int len, double **paired_seq, PartStruc *pstruc);


#endif /* THERMODYNAMIC_H_ */
