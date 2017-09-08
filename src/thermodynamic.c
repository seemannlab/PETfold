/*
 * thermodynamic.c
 *
 *  Created on: 06.02.2011
 *      Author: Stefan Seemann
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "thermodynamic.h"
#include "utils.h"
#include "fold_vars.h"
#include "fold.h"
#include "cofold.h"
#include "part_func.h"
#include "part_func_co.h"
#include "petfoldlibs.h"

/*
 * Extracts base pair probabilities from thermodynamic partition function of RNAfold
 * Calculates energy distributions without gaps
 * Constraints an evolutionary conserved base pair if both bases aren't gaps
 * Uses for gaps the average of base paired probabilities where both bases aren't gaps
 */
double **get_prob_paired_seq(Aln *align_gf, PartStruc *evocon)
{
	double **paired;
	SeqList *seqlist_gf;
	char *evocon_db;
	float e;
	double kT, sfact=1.07, min_en;
	char *structure;
	int i, d, k;
	int gappos, m;

	static int SCALELENGTH = 200;

	/* list of sequences without gaps */
	seqlist_gf = get_align_without_gaps(align_gf);
	//for (i=0; i < seqlist_gf->nr; i++) printf("%s\t%i\n",seqlist_gf->sequence[i], seqlist_gf->len[i]);

	/* initialize base paired probability matrix of thermodynamic model */
	paired = (double **)malloc((align_gf->len) * sizeof(double *));
	for (i=0; i<align_gf->len; i++)
		paired[i] = (double *)malloc((align_gf->len) * sizeof(double));

	for (d=0; d<align_gf->len; d++)
		for (k=0; k<align_gf->len; k++)
			paired[d][k] = .0;

	if (verbose_flag)
		printf("RNAFOLD:\n");
	/* thermodynamic model is switched on */
	if( setthermo_flag )
	{
		/* folding parameters */
		//temperature = 30.;      /* fold at 30C instead of the default 37C */
		dangles = 2;
		kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */
		noLonelyPairs = 1; /* disallow all pairs which can only occur as lonely pairs */

		for (i=0; i<align_gf->nr; i++)
		{
			/* write partial structure as constraint string usable with 'RNAfold -C' */
			evocon_db = get_constraint_string(evocon, seqlist_gf->origid[i], seqlist_gf->len[i]);
			if (verbose_flag)
				printf(" Constrained structure = %s\n", evocon_db);

			/* for longer sequences one should also set a scaling factor for
			   partition function folding, e.g: */
			if( seqlist_gf->len[i]>SCALELENGTH )
			{
				structure = (char *) space((unsigned) align_gf->len+1);
				fold_constrained = 0;
				min_en = fold(seqlist_gf->sequence[i], structure);
				pf_scale = exp(-(sfact*min_en)/kT/seqlist_gf->len[i]);
				free(structure);
			}
			init_pf_fold(seqlist_gf->len[i]);

			/* calculate partition function and base pair probabilities by constrained folding */
			fold_constrained = 1;	/* the structure string is interpreted on 'pf_fold' input as a list of constraints for the folding */
			e = pf_fold(seqlist_gf->sequence[i], evocon_db);

			/* write probabilities in gap-including matrix */
			for (d=1; d<=seqlist_gf->len[i]; d++)
				for(k=d+1; k<=seqlist_gf->len[i]; k++)
					paired[ seqlist_gf->origid[i][d-1] ][ seqlist_gf->origid[i][k-1] ] += pr[iindx[d]-k]; //sqrt(pr[iindx[d]-k]);

			if (verbose_flag)
				printf(" Sequence %2i structure = %s\n free energy of ensemble=%5.2f kcal/mol\n", i+1, evocon_db, e);

			free_pf_arrays();  /* free space allocated for pf_fold() */
			free(evocon_db);
		}

		/* divide the sum of base pair probabilities by the number of gap-free base pairs */
		for (d=0; d<align_gf->len; d++)
		{
			for (k=d+1; k<align_gf->len; k++) {
				/* count number of sequences with at least one gap in the base pair */
				gappos = 0;
				for( m=0; m<align_gf->nr; m++ ) {
					if( align_gf->sequence[m][d] == '-' || align_gf->sequence[m][k] == '-' )
						gappos++;
				}
				paired[d][k] /= align_gf->nr - gappos;
			}
		}
	}

	FreeSeqList(seqlist_gf);

	return paired;
}


/*
 * Computes unpaired probabilities as
 * Prob_unpaired(i) = 1 - SUM_j{Prob_paired(i,j)}
 */
double *get_prob_unpaired_seq(Aln *align_gf, double **paired_seq)
{
	double *single;
	double sum;
	int i, j;

	/* initialize single stranded probability array of thermodynamic model */
	single = (double *)malloc((align_gf->len) * sizeof(double));

	for (i=0; i<align_gf->len; i++)
		single[i] = .0;

	/* thermodynamic model is switched on */
	if( setthermo_flag )
	{
		for (i=0; i<align_gf->len; i++)
		{
			sum = 0;
			for (j=i+1; j<align_gf->len; j++)
				sum += paired_seq[i][j];
			for (j=0; j<i; j++)
				sum += paired_seq[j][i];
			single[i] = 1 - sum;
		}
	}

	return single;
}


/*
 * Extracts gap-free sequences from an alignment
 * and keeps the index relation between original and gap free sequences
 */
SeqList *get_align_without_gaps(Aln *align_gf)
{
	SeqList *seqlist_gf;
	int i, j, k;

	/* Initialize SeqList */
	seqlist_gf = (SeqList *)malloc(sizeof(SeqList));
	seqlist_gf->sequence = (char **)malloc(align_gf->nr * sizeof(char *));
	for (i=0; i<align_gf->nr; i++)
		seqlist_gf->sequence[i] = (char *)malloc((align_gf->len+1) * sizeof(char));
	seqlist_gf->origid = (int **)malloc(align_gf->nr * sizeof(int *));
	for (i=0; i<align_gf->nr; i++)
		seqlist_gf->origid[i] = (int *)malloc(align_gf->len * sizeof(int));
	seqlist_gf->len = (int *)malloc(align_gf->nr * sizeof(int));

	seqlist_gf->nr = align_gf->nr;

	for (i=0; i<align_gf->nr; i++)
	{
		k = 0;
		seqlist_gf->len[i] = 0;
		for (j=0; j<align_gf->len; j++)
		{
			if( align_gf->sequence[i][j] != '-' )
			{
				seqlist_gf->origid[i][k] = j;
				seqlist_gf->sequence[i][k++] = align_gf->sequence[i][j];
				seqlist_gf->len[i]++;
			}
			seqlist_gf->sequence[i][k] = '\0';
		}
	}

	return seqlist_gf;
}


/*
 * calculate probability of constrained (partial) structure in the thermodynamic model
 */
double get_fold_thermo_partial_prob(double econstraint, double eensemble)
{
	double kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */

	/* probability of partial structure in ensemble */
	return exp((eensemble-econstraint)/kT);
}


/*
 * calculate free energy of constrained ensemble (partial structure)
 * return arithmetric mean of energies
 */
double get_constraint_energy(Aln * align, PartStruc *evocon)
{
	SeqList *seqlist_gf;
	char *evocon_db;
	double e=0;
	double kT, sfact=1.07, min_en;
	char *structure;
	int i;

	static int SCALELENGTH = 200;

	/* list of sequences without gaps */
	seqlist_gf = get_align_without_gaps(align);
	//for (i=0; i < seqlist_gf->nr; i++) printf("%s\t%i\n",seqlist_gf->sequence[i], seqlist_gf->len[i]);

	/* folding parameters */
	//temperature = 30.;      /* fold at 30C instead of the default 37C */
	dangles = 2;
	kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */
	noLonelyPairs = 1; /* disallow all pairs which can only occur as lonely pairs */

	for (i=0; i<align->nr; i++)
	{
		/* write partial structure as constraint string usable with 'RNAfold -C' */
		evocon_db = get_constraint_string(evocon, seqlist_gf->origid[i], seqlist_gf->len[i]);

		/* for longer sequences one should also set a scaling factor for
		   partition function folding, e.g: */
		if( seqlist_gf->len[i]>SCALELENGTH )
		{
			structure = (char *) space((unsigned) align->len+1);
			fold_constrained = 0;	/* structure string isn't interpreted on 'fold' input as a list of constraints */
			min_en = fold(seqlist_gf->sequence[i], structure);
			pf_scale = exp(-(sfact*min_en)/kT/seqlist_gf->len[i]);
			free(structure);
		}
		init_pf_fold(seqlist_gf->len[i]);

		/* calculate partition function and base pair probabilities by constrained folding */
		fold_constrained = 1;	/* evocon_db string is interpreted on 'pf_fold' input as a list of constraints for the folding */
		e += (double) pf_fold(seqlist_gf->sequence[i], evocon_db);
		//printf(" Constrained sequence %2i structure = %s\n free energy of ensemble=%5.2f kcal/mol\n", i+1, evocon_db, e);

		free_pf_arrays();  /* free space allocated for pf_fold() */
		free(evocon_db);
	}

	FreeSeqList(seqlist_gf);

	return e/align->nr;
}


/*
 * calculate free energy of ensemble
 * return arithmetric mean of energies
  */
double get_ensemble_energy(Aln* align)
{
	SeqList *seqlist_gf;
	double e=0.;
	double kT, sfact=1.07, min_en;
	char *structure;
	int i;

	static int SCALELENGTH = 200;

	/* list of sequences without gaps */
	seqlist_gf = get_align_without_gaps(align);
	//for (i=0; i < seqlist_gf->nr; i++) printf("%s\t%i\n",seqlist_gf->sequence[i], seqlist_gf->len[i]);

	/* folding parameters */
	//temperature = 30.;      /* fold at 30C instead of the default 37C */
	dangles = 2;
	kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */
	fold_constrained = 0;
	noLonelyPairs = 1; /* disallow all pairs which can only occur as lonely pairs */

	for (i=0; i<align->nr; i++)
	{
		structure = (char *) space((unsigned) align->len+1);

		/* for longer sequences one should also set a scaling factor for
		   partition function folding, e.g: */
		if( seqlist_gf->len[i]>SCALELENGTH )
		{
			min_en = fold(seqlist_gf->sequence[i], structure);
			pf_scale = exp(-(sfact*min_en)/kT/seqlist_gf->len[i]);
		}
		init_pf_fold(seqlist_gf->len[i]);

		/* calculate partition function and base pair probabilities by constrained folding */
		e += (double) pf_fold(seqlist_gf->sequence[i], structure);
		//printf(" Ensemble sequence %2i free energy of ensemble=%5.2f kcal/mol\n", i+1, e);

		free_pf_arrays();  /* free space allocated for pf_fold() */
		free(structure);
	}

	FreeSeqList(seqlist_gf);

	return e/align->nr;
}

/*
 * Partition Function Cofolding (RNAcofold)
 * As for folding one RNA molecule, this computes the partition function of all possible structures and the base pair probabilities
 * http://www.tbi.univie.ac.at/RNA/RNAlib/Cofolding.html
 * Calculates energy distributions without gaps
 * Constraints partial structure base pairs if both bases aren't gaps
 * Uses for gaps the average of probabilities in column
 */
void get_cofold_prob_paired_seq(Aln *align, int len, double **paired_seq, PartStruc *pstruc)
{
	SeqList *seqlist_gf;
	char *evocon_db;
	//cofoldF *e;
	double kT, sfact=1.07, min_en;
	char *structure;
	int *gappos;
	int i, d, k;

	static int SCALELENGTH = 200;

	/* list of sequences without gaps */
	seqlist_gf = get_align_without_gaps(align);
	//for (i=0; i < seqlist_gf->nr; i++) printf("%s\t%i\n",seqlist_gf->sequence[i], seqlist_gf->len[i]);

	if (verbose_flag)
		printf("RNACOFOLD:\n");
	/* thermodynamic model is switched on */
	if( setthermo_flag )
	{
		/* folding parameters */
		//temperature = 30.;      /* fold at 30C instead of the default 37C */
		dangles = 2;
		kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */
		noLonelyPairs = 1; /* disallow all pairs which can only occur as lonely pairs */
		cut_point = len+1; /* marks the position (starting from 1) of the first nucleotide of the second molecule within the concatenated sequence */

		for (i=0; i<align->nr; i++)
		{
			/* write partial structure as constraint string usable with 'RNAfold -C' */
			evocon_db = get_constraint_string(pstruc, seqlist_gf->origid[i], seqlist_gf->len[i]);
			if (verbose_flag)
				printf(" Constrained structure = %s\n", evocon_db);

			/* for longer sequences one should also set a scaling factor for
			   partition function folding, e.g: */
			if( seqlist_gf->len[i]>SCALELENGTH )
			{
				structure = (char *) space((unsigned) align->len+1);
				fold_constrained = 0;	/* structure string isn't interpreted on 'fold' input as a list of constraints */
				min_en = cofold(seqlist_gf->sequence[i], structure);
				pf_scale = exp(-(sfact*min_en)/kT/seqlist_gf->len[i]);
				free(structure);
			}
			init_co_pf_fold(seqlist_gf->len[i]);

			/* calculate partition function of two interacting RNA molecules as well as base pair probabilities by constrained folding */
			fold_constrained = 1;	/* the structure string is interpreted on 'pf_fold' input as a list of constraints for the folding */
			cofoldF e = co_pf_fold(seqlist_gf->sequence[i], evocon_db);

			/* write probabilities in gap-including matrix */
			for (d=1; d<=seqlist_gf->len[i]; d++)
				for(k=d+1; k<=seqlist_gf->len[i]; k++)
					paired_seq[ seqlist_gf->origid[i][d-1] ][ seqlist_gf->origid[i][k-1] ] += pr[iindx[d]-k]; //sqrt(pr[iindx[d]-k]);

			if (verbose_flag)
				printf(" Sequence %2i structure = %s\n free energy of ensemble: F0AB=%5.2f kcal/mol; FA=%5.2f kcal/mol; FAB=%5.2f kcal/mol; FB=%5.2f kcal/mol; FcAB=%5.2f kcal/mol\n",
						i+1, evocon_db, (&e)->F0AB, (&e)->FA, (&e)->FAB, (&e)->FB, (&e)->FcAB);

			free_co_pf_arrays();  /* free space allocated for co_pf_fold() */
			free_co_arrays();
			free(evocon_db);
		}

		/*
		 * arithmetic mean of base pair probabilities
		 */

		/* count for each column of alignment the number of gaps */
		gappos = (int *)malloc(align->len * sizeof(int));
		for (d=0; d<align->len; d++)
		{
			gappos[d] = 0;
			for (k=0; k<align->nr; k++)
				if( align->sequence[k][d] == '-')
					gappos[d]++;
		}

		/* divide the sum of base pair probabilities by the number of gap-free sequences */
		for (d=0; d<align->len; d++)
			for (k=d+1; k<align->len; k++)
				paired_seq[d][k] /= align->nr - gappos[d];

		free(gappos);
	}

	FreeSeqList(seqlist_gf);
}
