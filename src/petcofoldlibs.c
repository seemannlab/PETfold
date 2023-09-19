/*
 * petcofoldlibs.c
 *
 *  Created on: 29.12.2012
 *      Author: Stefan Seemann
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "petfoldlibs.h"
#include "thermodynamic.h"
#include "evolutionary.h"
#include "petcofoldlibs.h"

/*
 * scan both alignments and take maximal common amount of identifiers (organisms)
 * exit if less than 3 organisms are common
 *
 */
int *intersect_alignments(Aln *align1, Aln *align2)
{
	int i, j, l1;
	int k = 0;

	int *cindex = (int *)malloc((align1->nr+align2->nr) * sizeof(int));
	for (i=0; i < align1->nr+align2->nr; i++)
		cindex[i] = -1;
	int offset = align1->nr;

	for (i=0; i < align1->nr; i++) {
		l1 = strlen(align1->identifier[i]);
		for (j=0; j < align2->nr; j++) {
			if( l1 == strlen(align2->identifier[j]) ) {
				if( ! strncmp(align1->identifier[i], align2->identifier[j], l1) ) {
					cindex[k] = i;
					cindex[offset+k] = j;
					k++;
					break;
				}
			}
		}
	}

	return cindex;
}

/*
 * get sub-alignment
 */
Aln *get_subalignment(Aln *align, int *cindex)
{
	Aln *subalign;
	int maxlen = align->nr;      /* Start with 30 lines - Dynamic line allocation maximum */

	/* Initialization of Aln Object */
	subalign = (Aln *)malloc(sizeof(Aln));
	subalign->identifier = (char **)malloc(maxlen * sizeof(char *));
	subalign->sequence = (char **)malloc(maxlen * sizeof(char *));

	int ptr = 0;
	int i;
	for (i=0; i<align->nr; i++) {
		if( cindex[i] == -1)
			break;
		subalign->identifier[i] = (char *)malloc((strlen(align->identifier[cindex[i]])+1) * sizeof(char));
		strcpy(subalign->identifier[i], align->identifier[cindex[i]]);
		subalign->sequence[i] = (char *)malloc((strlen(align->sequence[cindex[i]])+1) * sizeof(char));
		strcpy(subalign->sequence[i], align->sequence[cindex[i]]);
		ptr++;
	}

	subalign->identifier = (char **)realloc(subalign->identifier, ptr * sizeof(char *));
	subalign->sequence = (char **)realloc(subalign->sequence, ptr * sizeof(char *));
	                                             /* Free excess space */
	subalign->nr = ptr;
	subalign->len = align->len;

	return subalign;
}

/*
 * run PETfold for single sequences to get reliabilities and MEA structure
 */
PartStruc *run_petfold(Aln *align, double **paired_pet, double *single_pet, char *settree, PartStruc *setstruc)
{
	double *single_tree, *single_seq;
	double **paired_tree, **paired_seq;
	int i;
	int *pfold_struc;
	char *pet_db;
	PartStruc *evocon = NULL;
	DynEntry ***T;

	/* initialize base paired reliability matrix of evolutionary model */
	paired_tree = (double **)malloc((align->len) * sizeof(double *));
	for (i=0; i<align->len; i++)
		paired_tree[i] = (double *)malloc((align->len) * sizeof(double));
	/* initialize single stranded reliability array of evolutionary model */
	single_tree = (double *)malloc((align->len) * sizeof(double));

	/* calculate phylogenetic tree and structure reliabilities in evolutionary model -> Pfold */
	pfold_struc = get_prob_tree(align, paired_tree, single_tree, settree);

	//for(i=0;i<align->len;i++){fprintf(stderr,"%i ",pfold_struc[i]);};fprintf(stderr,"\n");
	//for (i=0; i<align_gf->len; i++) printf("%11.9f ", single_tree[i]); printf("\n");
	//int j; for (i=1; i<=align_gf->len; i++) { for (j=1; j<=align_gf->len; j++) printf("%2d %2d %16.8e\n", i,j,paired_tree[i-1][j-1]); }

	if( setevo_flag )
	{
		/* write Pfold predicted RNA secondary structure in an string */
		//char *pfold_db = (char *)calloc(align->len+1, sizeof(char));
		//get_dot_bracket(pfold_struc, align->len, pfold_db);
		//printf("%s\n", pfold_db);
		//free(pfold_db);

		/* evolutionary highly reliable partial structure */
		evocon = get_partial_struc(paired_tree, single_tree,
				pfold_struc, align->len, evocon_bp, evocon_ss);
		//for (i=0; i<evocon->bp_nr; i++) printf("%i %i\n", evocon->bp_left[i], evocon->bp_right[i]);
		//for (i=0; i<evocon->ss_nr; i++) printf("%i\n", evocon->ss[i]);
	}

	/* write partial structure as string */
	//char *evocon_db = partstruc2string(evocon, align->len);

	/* calculate structure probabilities in thermodynamic model -> RNAfold */
	paired_seq = get_prob_paired_seq(align, evocon);
	single_seq = get_prob_unpaired_seq(align, paired_seq);
	//for (i=0; i<align->len; i++) printf("%11.9f ", single_seq[i]); printf("\n");
	//int j; for (i=1; i<=align->len; i++) { for (j=1; j<=align->len; j++) printf("%2d %2d %1.15f\n", i,j,paired_seq[i-1][j-1]); }

	/* joint reliabilities of evolutionary reliabilities and thermodynamic probabilities */
	double *single_tmp = get_pet_single(single_tree, single_seq, align->len);
	double **paired_tmp = get_pet_paired(paired_tree, paired_seq, align->len);
	memcpy(single_pet, single_tmp, sizeof(double)*align->len);
	for (i=0; i<align->len; i++)
		memcpy(paired_pet[i], paired_tmp[i], sizeof(double)*align->len);
	free(single_tmp);
	FreeArray(paired_tmp, align->len);
	//for (i=0; i<align->len; i++) printf("%11.9f ", single_pet[i]); printf("\n");
	//int j; for (i=0; i<align->len; i++) { for (j=0; j<align->len; j++) printf("%11.9f ", paired_pet[i][j]); printf("\n"); } printf("\n");

	/*
	 * get MEA structure using a Nussinov-style algorithm
	 * replacing products in CYK algorithm by sums
	 */
	if( setstruc != NULL )
		T = nussinov(paired_pet, single_pet, setstruc, align->len);
	else
		T = nussinov(paired_pet, single_pet, evocon, align->len);
	//int j; for (j=0; j<align->len; j++) { printf("This is j=%i\n", j); for (i=0; i<align->len; i++) printf("%11.9f %11.9f %11.9f\n", T[i][j][0]->rel, T[i][j][1]->rel, T[i][j][2]->rel); } printf("\n");

	/* backtracking to find the consensus MEA structure */
	//pet_db_gf = backtracking(T, align_gf->len, pet_db_gf);
	pet_db = (char *)calloc(align->len+1, sizeof(char));
	backtracking(T, align->len, pet_db);

	free(pfold_struc);
	FreeDynEntry(T, align->len);

	FreeArray(paired_seq, align->len);
	FreeArray(paired_tree, align->len);

	/* return secondary structure as PartProb object */
	int temp = -1;
	PartStruc *struc_pet = mod_setstruct(pet_db, &temp, strlen(pet_db), 0);
	//printf("%s\n", partstruc2string(struc_pet, align->len));

	/* free memory */
	free(single_seq);
	free(single_tree);
	FreePartStruc(evocon);
	free(pet_db);

	return struc_pet;
}

/*
 * find base pairs in the PETfold structure that should be constraint as intra-molecular base pairs -> partial structure
 * increase petcon (max. allowed intra-molecular base-paired reliability of bases to be free for duplex binding)
 * until probability of partial structure compared to all possible structures is greater than partprob
 */
PartStruc *get_probable_partial_struc(Aln *align, double **paired_pet, PartStruc *struc_pet, char *tree, float petcon, double *partialprob)
{
    /* partial structure - selection of base pairs with large reliability */
	int *inputstruc = partstruc2str(struc_pet, align->len);

	/* test probability of partial structure in the structure ensemble
	   if probability is too low then repeat with increased maximal intra-molecular base-paired reliability
	   => petcon (base-paired probability threshold of base-pairs in partial structure)
	*/
	double *temp = (double *)calloc(align->len, sizeof(double));
	float mypetcon = petcon;
	PartStruc *evocon;
	double thermoconstraint, evolconstraint;

	/* ensemble energy in thermodynamic model */
	double thermoensemble = get_ensemble_energy(align);

	/* ensemble probability in evolutionary model */
	double evolensemble = get_ensemble_evolprob(align, tree);

	while( partialprob[0] <= partprob && partialprob[1] <= partprob && mypetcon < 1 ) {

		petcon = mypetcon;
		evocon = get_partial_struc(paired_pet, temp, inputstruc, align->len, mypetcon, 1);
		//int *evoconlist = partstruc2str(evocon, align->len);
		//for(i=0; i<align->len; i++) { printf("%i ", evoconlist[i]); }; printf("\n");

		if( extstem_flag ) {
			PartStruc *tempstruc = extend_partial_stems(evocon, inputstruc, paired_pet, align->len, mypetcon);
			FreePartStruc(evocon);
			evocon = tempstruc;
		}

		/* calculate probability of partial structure compared to all possible structures */
		/* in the thermodynamic model */
		thermoconstraint = get_constraint_energy(align, evocon);
		partialprob[0] = get_fold_thermo_partial_prob(thermoconstraint, thermoensemble);
		/* in the evolutionary model */
		evolconstraint = get_constraint_evolprob(align, tree, evocon);
		partialprob[1] = get_fold_evol_partial_prob(evolconstraint, evolensemble);

		if (verbose_flag)
			printf("Delta = %3.2f:\n"
					" average free energy of constrained ensemble=%5.4f kcal/mol and of ensemble=%5.4f kcal/mol => thermodyn partial struct prob = %6.5f\n"
					" evol score of constrained ensemble=%5.4f and of ensemble=%5.4f => evol partial struct prob = %6.5f\n",
					mypetcon, thermoconstraint, thermoensemble, partialprob[0], evolconstraint, evolensemble, partialprob[1]);

		/* increase maximal allowed intra-molecular base-paired reliability of bases to be free for duplex binding */
		if( mypetcon >= (float) 0.9 )
			mypetcon += 0.02;
		else
			mypetcon += 0.1;
	}
	free(temp);
	free(inputstruc);

	printf("Delta = %3.2f; thermodynamic partial structure probability = %6.5f; evolutionary partial structure probability = %6.5f\n", petcon, partialprob[0], partialprob[1]);

	return evocon;
}

/*
 * evolutionary scoring of RNA duplex assuming one unique tree
 */
int *get_duplex_prob_tree(Aln *align_gf1, Aln *align_gf2, double **duplex_paired_tree, double *duplex_single_tree, PartStruc *ss_partial1, PartStruc *ss_partial2, char *settree)
{
	int i, j;
	Aln *align;
	PartStruc *pstruc;

	/* Initialization of duplex Aln Object and concatenate sequences with a linker 'nnn' */
	align = (Aln *)malloc(sizeof(Aln));
	align->identifier = (char **)malloc(align_gf1->nr * sizeof(char *));
	align->sequence = (char **)malloc(align_gf1->nr * sizeof(char *));
	for (i=0; i<align_gf1->nr; i++) {
		align->identifier[i] = (char *)malloc(strlen(align_gf1->identifier[i])+1 * sizeof(char));
		strcpy(align->identifier[i], align_gf1->identifier[i]);
		align->sequence[i] = (char *)malloc((align_gf1->len+align_gf2->len+3)+1 * sizeof(char));
		strcpy(align->sequence[i], align_gf1->sequence[i]);
		strcpy(align->sequence[i]+align_gf1->len, "nnn");
		strcpy(align->sequence[i]+align_gf1->len+3, align_gf2->sequence[i]);
	}
	align->nr = align_gf1->nr;
	align->len = align_gf1->len+align_gf2->len+3;

	/* Initialization of duplex PartStruc Object and concatenate partial structures with a linker '...' */
	pstruc = (PartStruc *)malloc(sizeof(PartStruc));
	pstruc->bp_left = NULL;
	pstruc->bp_right = NULL;
	pstruc->bp_nr = 0;
	pstruc->ss = (int *)malloc((2*ss_partial1->bp_nr+2*ss_partial2->bp_nr+3) * sizeof(int));
	pstruc->ss_nr = 2*ss_partial1->bp_nr+2*ss_partial2->bp_nr+3;
	for (i=0; i<ss_partial1->bp_nr; i++) {
		pstruc->ss[i] = ss_partial1->bp_left[i];
		pstruc->ss[i+ss_partial1->bp_nr] = ss_partial1->bp_right[i];
	}
	for (i=0; i<3; i++)
		pstruc->ss[i+2*ss_partial1->bp_nr] = align_gf1->len+i;
	for (i=0; i<ss_partial2->bp_nr; i++) {
		pstruc->ss[i+2*ss_partial1->bp_nr+3] = align_gf1->len+3+ss_partial2->bp_left[i];
		pstruc->ss[i+2*ss_partial1->bp_nr+3+ss_partial2->bp_nr] = align_gf1->len+3+ss_partial2->bp_right[i];
	}
	qsort(pstruc->ss, pstruc->ss_nr, sizeof(int), sort);
	//for (i=0; i < align->nr; i++) printf("  %s\t%s\n", align->identifier[i], align->sequence[i]);
	//for (i=0; i < pstruc->ss_nr; i++) printf("  %i ", pstruc->ss[i]); printf("\n");

	/* initialize temporary duplex base paired reliability matrix */
	double **temp_duplex_paired = (double **)malloc((align_gf1->len+align_gf2->len+3) * sizeof(double *));
	for (i=0; i < align_gf1->len+align_gf2->len+3; i++)
		temp_duplex_paired[i] = (double *)malloc((align_gf1->len+align_gf2->len+3) * sizeof(double));
	/* initialize temporary duplex single stranded reliability array */
	double *temp_duplex_single = (double *)malloc((align_gf1->len+align_gf2->len+3) * sizeof(double));

	int *temp_evocon = get_prob_tree_constraint(align, temp_duplex_paired, temp_duplex_single, pstruc, settree);

	/* copy reliability matrices without linker sequence */
	for (i=0; i < align_gf1->len; i++) {
		duplex_single_tree[i] = temp_duplex_single[i];
		for (j=0; j < align_gf1->len; j++)
			duplex_paired_tree[i][j] = temp_duplex_paired[i][j];
		for (j=0; j < align_gf2->len; j++)
			duplex_paired_tree[i][j+align_gf1->len] = temp_duplex_paired[i][j+align_gf1->len+3];
	}
	for (i=0; i < align_gf2->len; i++) {
		duplex_single_tree[i+align_gf1->len] = temp_duplex_single[i+align_gf1->len+3];
		for (j=0; j < align_gf1->len; j++)
			duplex_paired_tree[i+align_gf1->len][j] = temp_duplex_paired[i+align_gf1->len+3][j];
		for (j=0; j < align_gf2->len; j++)
			duplex_paired_tree[i+align_gf1->len][j+align_gf1->len] = temp_duplex_paired[i+align_gf1->len+3][j+align_gf1->len+3];
	}

	/* copy structure without linker sequence */
	int *evocon = (int *)malloc((align_gf1->len+align_gf2->len) * sizeof(int));
	for (i=0; i < align_gf1->len; i++)
		evocon[i] = temp_evocon[i];
	for (i=0; i < align_gf2->len; i++)
		evocon[i+align_gf1->len] = temp_evocon[i+align_gf1->len+3];

	FreeArray(temp_duplex_paired, align_gf1->len+align_gf2->len+3);
	free(temp_duplex_single);
	FreeAln(align);
	FreePartStruc(pstruc);
	free(temp_evocon);

	return evocon;
}

/*
 * thermodynamic scoring of RNA duplex -> RNAcofold
 * calculate energy distributions without gaps, use for gaps the average of probabilities in column
 */
void get_duplex_prob_energy(Aln *align_gf1, Aln *align_gf2, double **duplex_paired_energy, double *duplex_single_energy, PartStruc *ss_partial1, PartStruc *ss_partial2)
{
	int i;
	Aln *align;
	PartStruc *pstruc;

	/* Initialization of duplex Aln Object and concatenate sequences without linker */
	align = (Aln *)malloc(sizeof(Aln));
	align->identifier = (char **)malloc(align_gf1->nr * sizeof(char *));
	align->sequence = (char **)malloc(align_gf1->nr * sizeof(char *));
	for (i=0; i<align_gf1->nr; i++) {
		align->identifier[i] = (char *)malloc(strlen(align_gf1->identifier[i])+1 * sizeof(char));
		strcpy(align->identifier[i], align_gf1->identifier[i]);
		align->sequence[i] = (char *)malloc((align_gf1->len+align_gf2->len)+1 * sizeof(char));
		strcpy(align->sequence[i], align_gf1->sequence[i]);
		strcpy(align->sequence[i]+align_gf1->len, align_gf2->sequence[i]);
	}
	align->nr = align_gf1->nr;
	align->len = align_gf1->len+align_gf2->len;

	/* Initialization of duplex PartStruc Object and concatenate partial structures without linker */
	pstruc = (PartStruc *)malloc(sizeof(PartStruc));
	pstruc->bp_left = NULL;
	pstruc->bp_right = NULL;
	pstruc->bp_nr = 0;
	pstruc->ss = (int *)malloc((2*ss_partial1->bp_nr+2*ss_partial2->bp_nr) * sizeof(int));
	pstruc->ss_nr = 2*ss_partial1->bp_nr+2*ss_partial2->bp_nr;
	for (i=0; i<ss_partial1->bp_nr; i++) {
		pstruc->ss[i] = ss_partial1->bp_left[i];
		pstruc->ss[i+ss_partial1->bp_nr] = ss_partial1->bp_right[i];
	}
	for (i=0; i<ss_partial2->bp_nr; i++) {
		pstruc->ss[i+2*ss_partial1->bp_nr] = align_gf1->len+ss_partial2->bp_left[i];
		pstruc->ss[i+2*ss_partial1->bp_nr+ss_partial2->bp_nr] = align_gf1->len+ss_partial2->bp_right[i];
	}
	qsort(pstruc->ss, pstruc->ss_nr, sizeof(int), sort);
	//for (i=0; i < align->nr; i++)	printf("  %s\t%s\n", align->identifier[i], align->sequence[i]);
	//for (i=0; i < pstruc->ss_nr; i++)	printf("  %i ", pstruc->ss[i]);	printf("\n");

	/* calculate base-paired probabilities in thermodynamic model -> RNAcofold */
	get_cofold_prob_paired_seq(align, align_gf1->len, duplex_paired_energy, pstruc);

	/* calculate single-stranded probabilities in thermodynamic model */
	double *temp = get_prob_unpaired_seq(align, duplex_paired_energy);
	for(i=0; i < align->len; i++)
		duplex_single_energy[i] = temp[i];
	free(temp);

	FreeAln(align);
	FreePartStruc(pstruc);
}

/*
 * duplex base-paired reliabilities are adapted by the probability of partial structure (intra-molecular constraint structure)
 */
void adjust_duplex_paired(double **duplex_paired, int len, double partialprob)
{
	int i, j;
	for (i=0; i < len; i++)
		for (j=i+1; j < len; j++)
			duplex_paired[i][j] *= partialprob;
}

/*
 * duplex single-stranded reliabilities are adapted by the probability of partial structure (intra-molecular constraint structure)
 */
void adjust_duplex_single(double *duplex_single, int len, double partialprob)
{
	int i;
	for (i=0; i < len; i++)
		duplex_single[i] *= partialprob;
}

/*
 * merge reliability 0 for bases with large intra-molecular base-paired reliability
 * and single-stranded duplex reliabilities for bases with low intra-molecular base-paired reliability
 */
void merge_single_fold_duplex_pet(double *duplex_single_pet, PartStruc *ss_partial1, PartStruc *ss_partial2, int len1)
{
	int i;

    /* set partial structure (intra-molecular base-pair constraints) unpaired reliability to 0 */
    for (i=0; i < ss_partial1->bp_nr; i++) {
        duplex_single_pet[ss_partial1->bp_left[i]] = 0.;
        duplex_single_pet[ss_partial1->bp_right[i]] = 0.;
    }
    for (i=0; i < ss_partial2->bp_nr; i++) {
    	duplex_single_pet[len1+ss_partial2->bp_left[i]] = 0.;
    	duplex_single_pet[len1+ss_partial2->bp_right[i]] = 0.;
    }
}

/*
 * merge reliability of partial structures (OR step 1 reliability if structure is given) for bases with large intra-molecular base-paired reliability
 * and base-paired duplex reliabilities for bases with low intra-molecular base-paired reliability
 */
void merge_paired_fold_duplex_pet(double **duplex_paired_pet, double *partialprob1, double *partialprob2,
		PartStruc *ss_partial1, PartStruc *ss_partial2, int len1)
{
	int i;

	/* partial structure sequence 1 = partial structure probability of sequence 1 */
    for (i=0; i < ss_partial1->bp_nr; i++)
    	duplex_paired_pet[ss_partial1->bp_left[i]][ss_partial1->bp_right[i]] = add(partialprob1[1], beta * partialprob1[0]) / 2;
    for (i=0; i < ss_partial2->bp_nr; i++)
    	duplex_paired_pet[len1+ss_partial2->bp_left[i]][len1+ss_partial2->bp_right[i]] = add(partialprob2[1], beta * partialprob2[0]) / 2;
}

/*
 * get maximum expected accuracy (MEA) structure using a Nussinov-style recursion for finding the most reliable RNA duplex structure
 * call the function 'nussinov'
 */
DynEntry ***nussinov_duplex(double **duplex_paired_pet, double *duplex_single_pet, PartStruc *pstruc, int len1, int len2)
{
	int i, j;

	/* initialize temporary duplex base paired reliability matrix with linker*/
	double **temp_paired = (double **)malloc((len1+len2+3) * sizeof(double *));
	for (i=0; i < len1+len2+3; i++)
		temp_paired[i] = (double *)malloc((len1+len2+3) * sizeof(double));
	/* initialize temporary duplex single stranded reliability array with linker */
	double *temp_single = (double *)malloc((len1+len2+3) * sizeof(double));
	/* copy reliability matrices with linker */
	for (i=0; i < len1; i++) {
		temp_single[i] = duplex_single_pet[i];
		for (j=0; j < len1; j++)
			temp_paired[i][j] = duplex_paired_pet[i][j];
		for (j=0; j < len2; j++)
			temp_paired[i][j+len1+3] = duplex_paired_pet[i][j+len1];
	}
	/* add linker */
	for (i=0; i<3; i++) {
		temp_single[i+len1] = 0.;
		for (j=0; j<len1+len2+3; j++)
			temp_paired[i+len1][j] = temp_paired[j][i+len1] = 0.;
	}
	for (i=0; i < len2; i++) {
		temp_single[i+len1+3] = duplex_single_pet[i+len1];
		for (j=0; j < len1; j++)
			temp_paired[i+len1+3][j] = duplex_paired_pet[i+len1][j];
		for (j=0; j < len2; j++)
			temp_paired[i+len1+3][j+len1+3] = duplex_paired_pet[i+len1][j+len1];
	}

	/* run Nussinov-style algorithm */
	DynEntry ***T = nussinov(temp_paired, temp_single, pstruc, len1+len2+3);

	FreeArray(temp_paired, len1+len2+3);
	free(temp_single);

	return T;
}

/*
 * remove linker of 3 positions between 1st and 2nd sequence
 * find inter-molecular base-pairs and assign to them '[' and ']'
 */
void remove_linker(char *pet_db_gf, char *pet_db_gf1, char *pet_db_gf2, int len1, int len2)
{
	int i, k, mlen;
	int *stack;

	mlen = (len1 >= len2)? len1 : len2;
	stack = (int *)malloc(mlen * sizeof(int));

	k=-1;
	for (i=0; i<len1; i++) {
		pet_db_gf1[i] = pet_db_gf[i];
		if( pet_db_gf[i] == '(' )
			stack[++k] = i;
		else {
			if( pet_db_gf[i] == ')' )
				k--;
		}
	}
	for (i=0; i<=k; i++)
		pet_db_gf1[stack[i]] = '[';

	k=-1;
	for (i=len2-1; i>=0; i--) {
		pet_db_gf2[i] = pet_db_gf[i+len1+3];
		if( pet_db_gf[i+len1+3] == ')' )
			stack[++k] = i;
		else {
			if( pet_db_gf[i+len1+3] == '(' )
				k--;
		}
	}
	for (i=0; i<=k; i++)
		pet_db_gf2[stack[i]] = ']';

	free(stack);
}

/*
 * stems of the partial structure get extended by inner base pairs and outer base pairs
 * if the average reliability is larger than the threshold PETCON of the inner or outer extended stem, respectively
 */
PartStruc *extend_partial_stems(PartStruc *partstruc, int *inputstruc, double **paired_pet, int len, float petcon)
{
	int i, k, outerleft, innerleft, newstem = 0, stem_start=-100, stem_l=0;
	float stempet, ext_sum, stem_sum=0.;

	int *extstruc = (int *)malloc(len * sizeof(int));
	for (i=0; i<len; i++)
		extstruc[i] = INFINITE_I;

	for (i=partstruc->bp_nr-1; i>=0; i--) {
		extstruc[partstruc->bp_left[i]] = partstruc->bp_right[i];
		extstruc[partstruc->bp_right[i]] = partstruc->bp_left[i];

		/* start of a stem */
		if( i==partstruc->bp_nr-1 || partstruc->bp_right[i] != partstruc->bp_right[i+1]+1 ) {
			newstem = 1;
			stem_start = i;
			stem_l = 1;
			stem_sum = paired_pet[partstruc->bp_left[i]][partstruc->bp_right[i]];
		}

		if( i==0 || partstruc->bp_right[i] != partstruc->bp_right[i-1]-1 ) {
			/* end of a stem */
			/* check outer base pairs */
			ext_sum = stem_sum;
			for (k=1; k<len-partstruc->bp_right[i]; k++) {
				outerleft = inputstruc[ partstruc->bp_right[i] + k ];
				if( outerleft != -1 && outerleft == partstruc->bp_left[i] - k ) {
					/* stem get extended by outer base pair if average reliability of extended stem is higher as threshold PETCON */
					stempet = ( ext_sum + paired_pet[partstruc->bp_left[i]-k][partstruc->bp_right[i]+k] ) / ( stem_l+k );
					if( stempet >= petcon ) {
						/* add outer base pair to partial structure */
						extstruc[partstruc->bp_left[i]-k] = partstruc->bp_right[i]+k;
						extstruc[partstruc->bp_right[i]+k] = partstruc->bp_left[i]-k;
						//printf("Add outer bp(%i,%i)\n", partstruc->bp_left[i]-k, partstruc->bp_right[i]+k);
						ext_sum += paired_pet[partstruc->bp_left[i]-k][partstruc->bp_right[i]+k];
					}
					else
						break;
				}
				else
					break;
			}
			/* check inner base pairs */
			ext_sum = stem_sum;
			for (k=1; k<partstruc->bp_right[stem_start]; k++) {
				innerleft = inputstruc[ partstruc->bp_right[stem_start] - k ];
				if( innerleft == partstruc->bp_left[stem_start] + k ) {
					/* stem get extended by inner base pair if average reliability of extended stem is higher as threshold PETCON */
					stempet = ( ext_sum + paired_pet[partstruc->bp_left[stem_start]+k][partstruc->bp_right[stem_start]-k] ) / ( stem_l+k );
					if( stempet >= petcon ) {
						/* add inner base pair to partial structure */
						extstruc[partstruc->bp_left[stem_start]+k] = partstruc->bp_right[stem_start]-k;
						extstruc[partstruc->bp_right[stem_start]-k] = partstruc->bp_left[stem_start]+k;
						//printf("Add inner bp(%i,%i)\n", partstruc->bp_left[stem_start]+k, partstruc->bp_right[stem_start]-k);
						ext_sum += paired_pet[partstruc->bp_left[stem_start]+k][partstruc->bp_right[stem_start]-k];
					}
					else
						break;
				}
				else
					break;
			}
		}
		else {
			/* stem continues */
			if( !newstem ) {
				stem_l++;
				stem_sum += paired_pet[partstruc->bp_left[i]][partstruc->bp_right[i]];
			}
			newstem = 0;
		}
	}

	/* return extended partial structure as PartStruc object */
	PartStruc *extpartstruc = str2partstruc(extstruc, len);
	if( verbose_flag ) {
		printf("partial structure:          %s\n", partstruc2string(partstruc, len));
		printf("extended partial structure: %s\n", partstruc2string(extpartstruc, len));
	}
	free(extstruc);

	return extpartstruc;
}

/*
 * sum over all single-stranded and base-pair reliabilities of the MEA consensus structure
 */
double get_score(char *pet_db_gf1, char *pet_db_gf2, double **duplex_paired_pet, double *duplex_single_pet, int len1, int len2)
{
	double score=0;
	int i, k, l, m, mlen, len;

	len = len1+len2;
	mlen = (len1 >= len2)? len1 : len2;
	int *stack1 = (int *)malloc(mlen/2 * sizeof(int));
	int *stack2 = (int *)malloc(mlen/2 * sizeof(int));
	int *stack3 = (int *)malloc(mlen * sizeof(int));

	k = l = m = -1;
	for (i=0; i<len1; i++) {
		switch( pet_db_gf1[i] )
		{
			case '.': score += 2*alpha * duplex_single_pet[i]/len;
					  //printf("%i:%5.4f ", i, duplex_single_pet[i]);
					  break;
			case '(': stack1[++k] = i; break;
			case ')': score += 2 * duplex_paired_pet[stack1[k--]][i]/len;
					  //printf("%i:%5.4f ", i, duplex_paired_pet[stack1[k+1]][i]);
					  break;
			case '{': stack2[++l] = i; break;
			case '}': score += 2 * duplex_paired_pet[stack2[l--]][i]/len;
			          //printf("%i:%5.4f ", i, duplex_paired_pet[stack2[l+1]][i]);
					  break;
			case '[': stack3[++m] = i; break;
		}
	}
	k = l = -1;
	for (i=0; i<len2; i++) {
		switch( pet_db_gf2[i] )
		{
			case '.': score += 2*alpha * duplex_single_pet[i+len1]/len;
					  //printf("%i:%5.4f ", i+len1, duplex_single_pet[i+len1]);
				      break;
			case '(': stack1[++k] = i; break;
			case ')': score += 2 * duplex_paired_pet[stack1[k--]+len1][i+len1]/len;
					  //printf("%i:%5.4f ", i+len1, duplex_paired_pet[stack1[k+1]+len1][i+len1]);
					  break;
			case '{': stack2[++l] = i; break;
			case '}': score += 2 * duplex_paired_pet[stack2[l--]+len1][i+len1]/len;
					  //printf("%i:%5.4f ", i+len1, duplex_paired_pet[stack2[l+1]+len1][i+len1]);
					  break;
			case ']': score += 2 * duplex_paired_pet[stack3[m--]][i+len1]/len;
					  //printf("%i:%5.4f ", i+len1, duplex_paired_pet[stack3[m+1]][i+len1]);
					  break;
		}
	}
	//printf("\n");

	free(stack1);
	free(stack2);
	free(stack3);

	return score;
}

/*
 * sort function to sort integers in ascending order
 * used as 4th argument by function qsort
 */
int sort(const void *x, const void *y)
{
  return (*(int*)x - *(int*)y);
}

/*
 * geometric mean
 */
double geometric_mean(double val1, double val2)
{
	return sqrt(val1 * val2);
}

/*
 * change character 'a' to 'b' in string
 */
char *substitute_char(char *string, char a, char b, int len)
{
	int i;

	for (i=0; i<len; i++)
		if( string[i] == a )
			string[i] = b;

	return string;
}

/*
 * help output if help_flag is set or not exactly two FASTA-files as input
 */
void usage_petcofold()
{
	printf("PETcofold v3.4\n"
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
	 "                             ... step 1 constraints intramolecular ('(',')') and step 2 intermolecular structure ('[',']')\n"
     "  -n --setstruc1 <struc|file>... finds reliable base pairs in given structure of 1st alignment in dot-bracket notation\n"
     "  -m --setstruc2 <struc|file>... finds reliable base pairs in given structure of 2nd alignment in dot-bracket notation\n"
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
	 "\n"
	 "You may set the environment variable PETFOLDBIN to path of files scfg.rate and article.grm:\n"
	 "export PETFOLDBIN=<PATH>\n"
	 "\n"
	);

	exit(1);
}
