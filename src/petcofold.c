/*
 * petcofold.c
 *
 *  Created on: 17.06.2012
 *      Author: Stefan Seemann
 */

#define MODUL0
#include  <stdio.h>
#include  <math.h>
#include  <string.h>
#include  <stdlib.h>
#include  <getopt.h>
#include  <unistd.h>
#include  <time.h>
#include  "petcofoldlibs.h"
#include  "petfoldlibs.h"
#include  "thermodynamic.h"
#include  "evolutionary.h"

int main(int argc, char **argv)
{
	/* PETfold parameters */
	gap = 0.25;   /* maximal allowed percent of gaps in an alignment column; deletion of column if rate is equal or higher */
	evocon_bp = 0.9;   /* minimal base pair reliability for evolutionary constraints */
	evocon_ss = 1.0;   /* minimal single stranded reliablity for evolutionary constraints */
	beta = 1.;   /* weighting factor for thermodynamic overlap */
	alpha = 0.2;   /* weighting factor for single stranded reliabilities divided by 2 (alpha <= 0.5) */

	/* PETcofold additional parameters */
	float petcon = 0.9;   /* threshold for PET-model constraints to be free for RNA duplex: maximal allowed intra-molecular base-paired reliability */
	partprob = 0.1;   /* lowest accepted probability of single partial structure */
	extstem_flag = 0;   /* constrained stems get extended by base pairs if the average reliability of the extended stem is larger than PETCON */
	//static int nolp_flag = 0;   /* RNA(co)fold option: disallows pairs that can only occur isolated */

	/* PETcofold input and output */
	char *fasta1 = "";
	char *fasta2 = "";
	char *settree = "";  /* given phylogenetic tree in Newick tree format */
	char *setstruc = "";   /* given duplex structure in dot-bracket notation including '&' */
	char *setstruc1 = "";   /* given secondary structure of 1st alignment in dot-bracket notation */
	char *setstruc2 = "";   /* given secondary structure of 2nd alignment in dot-bracket notation */
	char *ppfile = "";   /* name of pp-file */
	static int war_flag = 0;   /* fasta format output */;
	static int help_flag = 0;   /* help */
	static int intermol_flag = 0;   /* structure output of intermolecular base pairs */
	verbose_flag = 0;   /* long output */

	/* model complexity */
	setevo_flag = 1;   /* usage of evolutionary reliabilities (performed by Pfold) */
	setthermo_flag = 1;   /* usage of thermodynamic probabilities (performed by RNAfold) */

	Aln *align1, *align2, *align_gf1, *align_gf2;
	int *gap_col1, *gap_col2, gapnr1, gapnr2;
	char *pfold_db1, *pfold_db2;
	PartStruc *setstructure1 = NULL;
	PartStruc *setstructure2 = NULL;
	PartStruc *pstruc = NULL;
	DynEntry ***T;
	char *pet_db_gf, *pet_db1, *pet_db2, *pet_db_gf1, *pet_db_gf2;

	int cmd, i;
	double score;

	/* parsing long options */
	while(1)
	{
		static struct option long_options[] =
		{
			{"setevo", no_argument, &setevo_flag, 1},
			{"noevo", no_argument, &setevo_flag, 0},
			{"setthermo", no_argument, &setthermo_flag, 1},
			{"nothermo", no_argument, &setthermo_flag, 0},
			{"war", no_argument, &war_flag, 1},
			{"help", no_argument, &help_flag, 1},
			{"verbose", no_argument, &verbose_flag, 1},
			{"intermol", no_argument, &intermol_flag, 1},
			{"settree", required_argument, 0, 't'},
			{"fasta", required_argument, 0, 'f'},
			{"setgap", required_argument, 0, 'g'},
			{"setevocon_bp", required_argument, 0, 'p'},
			{"setevocon_ss", required_argument, 0, 'u'},
			{"setbeta", required_argument, 0, 'b'},
			{"setstruc", required_argument, 0, 's'},
			{"setstruc1", required_argument, 0, 'n'},
			{"setstruc2", required_argument, 0, 'm'},
			{"setalpha", required_argument, 0, 'a'},
			{"setdelta", required_argument, 0, 'd'},
			{"setgamma", required_argument, 0, 'i'},
			{"extstem", no_argument, &extstem_flag, 1},
			//{"nolp", no_argument, &nolp_flag, 1},
			{"ppfile", required_argument, 0, 'r'},
			{0, 0, 0, 0}
		};

		/* print help text */
		if (help_flag)
			usage_petcofold();

		/* getopt_long stores the option index here. */
		int option_index = 0;

		cmd = getopt_long(argc, argv, "f:g:t:p:u:b:s:n:m:a:d:i:r:", long_options, &option_index);

		/* Detect the end of the options. */
		if (cmd == -1)
			break;

		switch(cmd)
		{
			case 0:   break;
			case 'f': if (!strlen(fasta1))
						  fasta1 = optarg;
			          else if(!strlen(fasta2))
						  fasta2 = optarg;
			          else
			        	  usage_petcofold();
			          break;
			case 'g': gap = atof(optarg); break;
			case 't': settree = optarg; break;
			case 'p': evocon_bp = atof(optarg); break;
			case 'u': evocon_ss = atof(optarg); break;
			case 'b': beta = atof(optarg); break;
			case 's': setstruc = optarg; break;
			case 'n': setstruc1 = optarg; break;
			case 'm': setstruc2 = optarg; break;
			case 'a': alpha = atof(optarg); break;
			case 'd': petcon = atof(optarg); break; /* delta */
			case 'i': partprob = atof(optarg); break; /* gamma */
			case 'r': ppfile = optarg; break;
			default:  abort();
			/* no break */
		}
	}

	/* parameters */
	if (verbose_flag) {
		printf("PARAMETERS:\n");
		printf(" maximal percent of gaps per column (--setgap) = %3.2f\n", gap);
		printf(" minimal base pair reliability for evolutionary constraints (--setevocon_bp) = %3.2f\n", evocon_bp);
		printf(" minimal single stranded reliability for evolutionary constraints (--setevocon_ss) = %3.2f\n", evocon_ss);
		printf(" weighting factor for thermodynamic overlap (--setbeta) = %3.2f\n", beta);
		printf(" weighting factor for single stranded reliabilities divided by 2 (--setalpha) = %3.2f\n", alpha);
		printf(" maximal allowed intra-molecular base-paired reliability to be accessible for RNA duplex (--setdelta) = %3.2f\n", petcon);
		printf(" minimal probability of single partial structure (--setgamma) = %3.2f\n", partprob);
	}

	/* read input alignment */
	if (!strlen(fasta1) || !strlen(fasta2))
		usage_petcofold();
	align1 = get_alignment(fasta1);
	/* check that all sequences have the same length -> RNA alignment */
	for (i = 1; i < align1->nr; i++)
		if( strlen(align1->sequence[i]) != align1->len ) {
			printf("The first file doesn't look like a FASTA file of aligned RNA sequences - sorry.\n");
			exit(1);
		}
	align2 = get_alignment(fasta2);
	/* check that all sequences have the same length -> RNA alignment */
	for (i = 1; i < align2->nr; i++)
		if( strlen(align2->sequence[i]) != align2->len ) {
			printf("This second file doesn't look like a FASTA file of aligned RNA sequences - sorry.\n");
			exit(1);
		}

	//for (i=0; i < align1->nr; i++) printf("%s\t%s\n", align1->identifier[i], align1->sequence[i]);
	//for (i=0; i < align2->nr; i++) printf("%s\t%s\n", align2->identifier[i], align2->sequence[i]);
	if (verbose_flag) {
		printf("INPUT:\n ALIGNMENT 1:\n number sequences = %i\n sequence length = %i\n", align1->nr, align1->len);
		printf(" ALIGNMENT 2:\n number sequences = %i\n sequence length = %i\n", align2->nr, align2->len);
	}

	/* remove gap columns from alignment */
	gap_col1 = (int *)malloc((align1->len) * sizeof(int));
	gap_col2 = (int *)malloc((align2->len) * sizeof(int));
	align_gf1 = delete_gap_columns(align1, gap_col1, gap);
	align_gf2 = delete_gap_columns(align2, gap_col2, gap);

	/* get intersection of identifiers in both alignments */
	int *cindex = intersect_alignments(align_gf1, align_gf2);
	Aln *subalign_gf1 = get_subalignment(align_gf1, cindex);
	if (!subalign_gf1->nr) {
		fprintf(stderr, "Error <fasta>: "
			"Intersection of fasta files results in zero common identifiers!\n");
		exit(1);
	}
	if (subalign_gf1->nr < 3) {
		printf("WARNING: PETcofold is cancelled as your input consists of ONLY %i sequence(s) with common identifiers.\n" ,
			subalign_gf1->nr);
		printf("Please provide two alignments with at least 3 common identifiers.\n");
		exit(1);
	}

	Aln *subalign_gf2 = get_subalignment(align_gf2, cindex+align_gf1->nr);
	free(cindex);
	FreeAln(align_gf1);
	FreeAln(align_gf2);

	gapnr1 = align1->len-subalign_gf1->len;
	gapnr2 = align2->len-subalign_gf2->len;

	if (verbose_flag)
	{
		printf(" Common identifiers =");
		for (i=0; i < subalign_gf1->nr; i++)
			printf(" %s,", subalign_gf1->identifier[i]);
		printf("\n");
		printf(" ALIGNMENT 1:\n gap columns = ");
		for (i=0; i < gapnr1; i++)
			printf("%i, ", gap_col1[i]+1);
		printf("\n gap column free alignment = \n");
		for (i=0; i < subalign_gf1->nr; i++)
			printf("  %s\t%s\n", subalign_gf1->identifier[i], subalign_gf1->sequence[i]);
		printf(" ALIGNMENT 2:\n gap columns = ");
		for (i=0; i < gapnr2; i++)
			printf("%i, ", gap_col2[i]+1);
		printf("\n gap column free alignment = \n");
		for (i=0; i < subalign_gf2->nr; i++)
			printf("  %s\t%s\n", subalign_gf2->identifier[i], subalign_gf2->sequence[i]);
	}

	/* open tree file in Newick format */
	char *tree = (char *)calloc(BUFFERSIZE, sizeof(char *));
	if (strlen(settree)) {
		settree = read_file(settree);
		strcpy(tree, settree);
		//tree = (char *)realloc(tree, strlen(tree)+3 * sizeof(char *));
	}
	/* open structure file in dot-bracket format */
	if (strlen(setstruc)) {
		setstruc = read_file(setstruc);
		petcon = partprob = 0.;
		char *temp_struc1 = (char *)malloc(align1->len+1);
		char *temp_struc2 = (char *)malloc(align2->len+1);
		strncpy(temp_struc1, setstruc, align1->len);
		temp_struc1[align1->len+1] = '\0';
		for( i=0; i<align1->len; i++ )
			if (temp_struc1[i] == '[' || temp_struc1[i] == ']')
				temp_struc1[i] = '.';
		strcpy(temp_struc2, &setstruc[align1->len+1]);
		for( i=0; i<align2->len; i++ )
			if (temp_struc2[i] == '[' || temp_struc2[i] == ']')
				temp_struc2[i] = '.';
		setstructure1 = mod_setstruct(temp_struc1, gap_col1, align1->len, gapnr1);
		setstructure2 = mod_setstruct(temp_struc2, gap_col2, align2->len, gapnr2);
		/* test for matching interaction */
		char *temp_inter = (char *)malloc(align1->len+align2->len+4);
		strncpy(temp_inter, setstruc, align1->len);
		strncpy(&temp_inter[align1->len],"...",3);
		strcpy(&temp_inter[align1->len+3],&setstruc[align1->len+1]);
		int inter = 0;
		for( i=0; i<align1->len; i++ )
			if (temp_inter[i] == '[')
				inter++;
			else
				temp_inter[i] = '.';
		for( i=0; i<align2->len; i++ )
			if (temp_inter[i+align1->len+3] == ']')
				inter--;
			else
				temp_inter[i+align1->len+3] = '.';
		if (inter) {
			fprintf(stderr, "Error <setstruc>: "
							"Number of left and right parenthesis are unequal for inter-molecular base pairs!\n");
			exit(1);
		}
		int *gap_col = (int *)calloc((gapnr1+gapnr2), sizeof(int *));
		memcpy(gap_col, gap_col1, gapnr1 * sizeof(int *));
		memcpy(gap_col + gapnr1, gap_col2, gapnr2 * sizeof(int *));
		for( i=0; i<gapnr2; i++ )
			gap_col[gapnr1+i] += align1->len+3;
		pstruc = mod_setstruct(temp_inter, gap_col, align1->len+align2->len+3, gapnr1+gapnr2);
		//printf("SETSTRUCINTER: %s\n", partstruc2string(pstruc,subalign_gf1->len+subalign_gf2->len+2));
		free(temp_struc1);
		free(temp_struc2);
		free(temp_inter);
		free(gap_col);
	}
	else {
		if (strlen(setstruc1) && !strlen(setstruc) )
			setstruc1 = read_file(setstruc1);
		if (strlen(setstruc2) && !strlen(setstruc) )
			setstruc2 = read_file(setstruc2);
	}

	/* PETfold for single sequences */

	/* initialize base paired reliability matrices */
	double **paired_pet1 = (double **)calloc(subalign_gf1->len, sizeof(double *));
	for (i=0; i<subalign_gf1->len; i++)
		paired_pet1[i] = (double *)calloc(subalign_gf1->len, sizeof(double));
	double **paired_pet2 = (double **)calloc(subalign_gf2->len, sizeof(double *));
	for (i=0; i<subalign_gf2->len; i++)
		paired_pet2[i] = (double *)calloc(subalign_gf2->len, sizeof(double));
	/* initialize single stranded reliability arrays */
	double *single_pet1 = (double *)calloc(subalign_gf1->len, sizeof(double));
	double *single_pet2 = (double *)calloc(subalign_gf2->len, sizeof(double));

	printf("PETfold Alignment 1:\n");
	/* parse input structure if given */
	if( strlen(setstruc1) && !strlen(setstruc) )
		setstructure1 = mod_setstruct(setstruc1, gap_col1, align1->len, gapnr1);

	/* run PETfold for 1st alignment */
	PartStruc *struc_pet1 = run_petfold(subalign_gf1, paired_pet1, single_pet1, tree, setstructure1);

	/* get partial intra-molecular structure that fulfills petcon and partprob */
	if (verbose_flag)
		printf("PARTIAL STRUCTURE 1:\n");
	double partialprob1[2] = { 0. };
	PartStruc *ss_partial1 = get_probable_partial_struc(subalign_gf1, paired_pet1, struc_pet1, tree, petcon, partialprob1);

	/* write partial probabilities and partial structures */
	char *inputstruc1 = partstruc2string(ss_partial1, subalign_gf1->len);
	char *partstruc1 = (char *)calloc(align1->len+1, sizeof(char));
	include_gap_columns(inputstruc1, gap_col1, subalign_gf1->len, gapnr1, partstruc1);
	//int *inputstruc1 = partstruc2str(ss_partial1, subalign_gf1->len);
	printf("Partial Structure 1:\t%s\n", partstruc1);
	free(inputstruc1);
	free(partstruc1);

	printf("PETfold Alignment 2:\n");
	/* parse input structure if given */
	if( strlen(setstruc2) && !strlen(setstruc) )
		setstructure2 = mod_setstruct(setstruc2, gap_col2, align2->len, gapnr2);

	/* run PETfold for 2nd alignment */
	if ( !strlen(settree) )
		memset(tree, '\0', sizeof(tree));
	else
		if ( strstr(settree, ":") == NULL )
			strcpy(tree, settree);

	PartStruc *struc_pet2 = run_petfold(subalign_gf2, paired_pet2, single_pet2, tree, setstructure2);

	/* get partial intra-molecular structure that fulfills petcon and partprob */
	if (verbose_flag)
		printf("PARTIAL STRUCTURE 2:\n");
	double partialprob2[2] = { 0. };
	PartStruc *ss_partial2 = get_probable_partial_struc(subalign_gf2, paired_pet2, struc_pet2, tree, petcon, partialprob2);

	/* write partial probabilities and partial structures */
	char *inputstruc2 = partstruc2string(ss_partial2, subalign_gf2->len);
	char *partstruc2 = (char *)calloc(align2->len+1, sizeof(char));
	include_gap_columns(inputstruc2, gap_col2, subalign_gf2->len, gapnr2, partstruc2);
	//int *inputstruc2 = partstruc2str(ss_partial2, subalign_gf2->len);
	printf("Partial Structure 2:\t%s\n", partstruc2);
	free(inputstruc2);
	free(partstruc2);

	/* PETcofold for RNA duplex binding */

	/* initialize duplex base paired reliability matrix of evolutionary model */
	double **duplex_paired_tree = (double **)calloc(subalign_gf1->len+subalign_gf2->len, sizeof(double *));
	for (i=0; i < subalign_gf1->len+subalign_gf2->len; i++)
		duplex_paired_tree[i] = (double *)calloc(subalign_gf1->len+subalign_gf2->len, sizeof(double));
	/* initialize single stranded reliability array of evolutionary model */
	double *duplex_single_tree = (double *)calloc(subalign_gf1->len+subalign_gf2->len, sizeof(double));

	if (verbose_flag)
		printf("COFOLDING:\n");
	/* evolutionary scoring adapts duplex reliabilities by evolutionary probability of partial structure
	 * (intra-molecular constrained structure):
	 * Pr_bp(step2) = P(E AND partprob(1,2)) = P(E | partprob(1,2)) x P(partprob(1,2))
	 */
	if ( !strlen(settree) )
		memset(tree, '\0', sizeof(tree));
	else
		if ( strstr(settree, ":") == NULL )
			strcpy(tree, settree);

	int *temp_struc = get_duplex_prob_tree(subalign_gf1, subalign_gf2, duplex_paired_tree, duplex_single_tree, ss_partial1, ss_partial2, tree);
	free(temp_struc);
	//if ( !strlen(settree) || strstr(settree, ":") == NULL )
	tree = (char *)realloc(tree, strlen(tree)+1 * sizeof(char *));

	//int *duplex_pfold_struc = get_duplex_prob_tree(subalign_gf1, subalign_gf2, duplex_paired_tree, duplex_single_tree, ss_partial1, ss_partial2, tree);
	//for (i=0; i<subalign_gf1->len+subalign_gf2->len; i++) printf("%i ", duplex_pfold_struc[i]); printf("\n");
	adjust_duplex_paired(duplex_paired_tree, subalign_gf1->len+subalign_gf2->len, geometric_mean(partialprob1[1], partialprob2[1]));
	adjust_duplex_single(duplex_single_tree, subalign_gf1->len+subalign_gf2->len, geometric_mean(partialprob1[1], partialprob2[1]));

	/* initialize duplex base paired reliability matrix of thermodynamic model */
	double **duplex_paired_energy = (double **)calloc(subalign_gf1->len+subalign_gf2->len, sizeof(double *));
	for (i=0; i < subalign_gf1->len+subalign_gf2->len; i++)
		duplex_paired_energy[i] = (double *)calloc(subalign_gf1->len+subalign_gf2->len, sizeof(double));
	/* initialize single stranded reliability array of thermodynamic model */
	double *duplex_single_energy = (double *)calloc(subalign_gf1->len+subalign_gf2->len, sizeof(double));

	/* thermodynamic scoring adapts duplex reliabilities by thermodynamic probability of partial structure
	 * (intra-molecular constrained structure):
	 * Pr_bp(step2) = P(E AND partprob(1,2)) = P(E | partprob(1,2)) x P(partprob(1,2))
	 */
	get_duplex_prob_energy(subalign_gf1, subalign_gf2, duplex_paired_energy, duplex_single_energy, ss_partial1, ss_partial2);
	adjust_duplex_paired(duplex_paired_energy, subalign_gf1->len+subalign_gf2->len, geometric_mean(partialprob1[0], partialprob2[0]));
	adjust_duplex_single(duplex_single_energy, subalign_gf1->len+subalign_gf2->len, geometric_mean(partialprob1[0], partialprob2[0]));

	/* calculate probabilities of the PETfold model */
	double *duplex_single_pet = get_pet_single(duplex_single_tree, duplex_single_energy, subalign_gf1->len+subalign_gf2->len);
	double **duplex_paired_pet = get_pet_paired(duplex_paired_tree, duplex_paired_energy, subalign_gf1->len+subalign_gf2->len);

	/* merging of partial structure probability for bases with large intra-molecular base paired reliabilities
	 * and inter-molecular reliabilities for the other bases:
	 * Pr_bp(step2,partial_struc_1) = partprob1
	 * Pr_bp(step2,partial_struc_2) = partprob2
	 */
	merge_single_fold_duplex_pet(duplex_single_pet, ss_partial1, ss_partial2, subalign_gf1->len);
	merge_paired_fold_duplex_pet(duplex_paired_pet, partialprob1, partialprob2, ss_partial1, ss_partial2, subalign_gf1->len);
	//int j; for( i=0; i < subalign_gf1->len+subalign_gf2->len; i++) printf("%5.4f %5.4f %5.4f :: ", duplex_single_pet[i], duplex_single_tree[i], duplex_single_energy[i]); printf("\n");
	//for( i=0; i < subalign_gf1->len+subalign_gf2->len; i++) { for( j=i+1; j < subalign_gf1->len+subalign_gf2->len; j++)	printf("%5.4f %5.4f %5.4f :: ", duplex_paired_pet[i][j], duplex_paired_tree[i][j], duplex_paired_energy[i][j]); printf("\n"); }

	/* MEA-CYK to find the most reliable RNA duplex structure */
	if( !strlen(setstruc) ) {
		/* Initialization of duplex PartStruc Object and concatenate partial structures with linker */
		pstruc = (PartStruc *)malloc(sizeof(PartStruc));
		pstruc->bp_left = NULL;
		pstruc->bp_right = NULL;
		pstruc->bp_nr = 0;
		pstruc->ss = (int *)malloc((2*ss_partial1->bp_nr+2*ss_partial2->bp_nr+3) * sizeof(int));
		pstruc->ss_nr = 2*ss_partial1->bp_nr+2*ss_partial2->bp_nr+3;
		/* constrain base pairs in partial structures as unpaired */
		for (i=0; i<ss_partial1->bp_nr; i++) {
			pstruc->ss[i] = ss_partial1->bp_left[i];
			pstruc->ss[i+ss_partial1->bp_nr] = ss_partial1->bp_right[i];
		}
		/* add linker */
		for (i=0; i<3; i++)
			pstruc->ss[i+2*ss_partial1->bp_nr] = subalign_gf1->len+i;
		for (i=0; i<ss_partial2->bp_nr; i++) {
			pstruc->ss[i+2*ss_partial1->bp_nr+3] = subalign_gf1->len+ss_partial2->bp_left[i]+3;
			pstruc->ss[i+2*ss_partial1->bp_nr+ss_partial2->bp_nr+3] = subalign_gf1->len+ss_partial2->bp_right[i]+3;
		}
		qsort(pstruc->ss, pstruc->ss_nr, sizeof(int), sort);
	}
	T = nussinov_duplex(duplex_paired_pet, duplex_single_pet, pstruc, subalign_gf1->len, subalign_gf2->len);
	//int j; for (j=0; j<subalign_gf1->len+subalign_gf2->len; j++) { printf("This is j=%i\n", j); for (i=0; i<subalign_gf1->len+subalign_gf2->len; i++) printf("%11.9f %11.9f %11.9f\n", T[i][j][0]->rel, T[i][j][1]->rel, T[i][j][2]->rel); } printf("\n");

	/* backtracking to find the consensus MEA structure */
	pet_db_gf = (char *)calloc(subalign_gf1->len+subalign_gf2->len+4, sizeof(char));
	score = backtracking(T, subalign_gf1->len+subalign_gf2->len+3, pet_db_gf);

	/* remove linker */
	pet_db_gf1 = (char *)calloc(subalign_gf1->len+1, sizeof(char));
	pet_db_gf2 = (char *)calloc(subalign_gf2->len+1, sizeof(char));
	remove_linker(pet_db_gf, pet_db_gf1, pet_db_gf2, subalign_gf1->len, subalign_gf2->len);

	/* merge partial structure and duplex structure */
	for( i=0; i<ss_partial1->bp_nr; i++ ) {
		pet_db_gf1[ss_partial1->bp_left[i]] = '{';
		pet_db_gf1[ss_partial1->bp_right[i]] = '}';
	}
	for( i=0; i<ss_partial2->bp_nr; i++ ) {
		pet_db_gf2[ss_partial2->bp_left[i]] = '{';
		pet_db_gf2[ss_partial2->bp_right[i]] = '}';
	}

	/* include gap columns */
	pet_db1 = (char *)calloc(align1->len+1, sizeof(char));
	pet_db2 = (char *)calloc(align2->len+1, sizeof(char));
	include_gap_columns(pet_db_gf1, gap_col1, subalign_gf1->len, gapnr1, pet_db1);
	include_gap_columns(pet_db_gf2, gap_col2, subalign_gf2->len, gapnr2, pet_db2);

	/* write pp-file */
	if( strlen(ppfile) ) {
		int *gap_col = (int *)malloc((align1->len+align2->len) * sizeof(int));
		for( i=0; i<gapnr1; i++ )
			gap_col[i] = gap_col1[i];
		for( i=0; i<gapnr2; i++ )
			gap_col[i+gapnr1] = gap_col2[i]+align1->len;
		create_ppfile(duplex_paired_pet, duplex_single_pet, gap_col, subalign_gf1->len+subalign_gf2->len, gapnr1+gapnr2, ppfile);
		free(gap_col);
	}

	/* calculate the score as sum of reliabilities in consensus duplex structure */
	score = get_score(pet_db_gf1, pet_db_gf2, duplex_paired_pet, duplex_single_pet, subalign_gf1->len, subalign_gf2->len);

	/* inter-molecular base pairs */
	char *intermol1 = (char *)calloc(align1->len+1, sizeof(char));
	for(i=0; i<align1->len; i++) {
		if (pet_db1[i]=='.' || pet_db1[i]=='-' || pet_db1[i]=='[')
			intermol1[i] = pet_db1[i];
		else
			intermol1[i] = '.';
	}
	intermol1[i] = '\0';
	char *intermol2 = (char *)calloc(align2->len+1, sizeof(char));
	for(i=0; i<align2->len; i++) {
		if (pet_db2[i]=='.' || pet_db2[i]=='-' || pet_db2[i]==']')
			intermol2[i] = pet_db2[i];
		else
			intermol2[i] = '.';
	}
	intermol2[i] = '\0';

	/* output */
	if( war_flag )
	{
		for (i=0; i < align1->nr; i++)
			printf(">%s\n%s&%s\n", align1->identifier[i], align1->sequence[i], align2->sequence[i]);
		if (intermol_flag)
			printf(">structure\n%s&%s\n", intermol1, intermol2);
		else
			printf(">structure\n%s&%s\n", pet_db1, pet_db2);
	}
	else
	{
		if (setevo_flag)
		{
			pfold_db1 = (char *)calloc(align1->len+1, sizeof(char));
			pfold_db2 = (char *)calloc(align2->len+1, sizeof(char));
			char *temp_pstruc1 = partstruc2string(struc_pet1, subalign_gf1->len);
			char *temp_pstruc2 = partstruc2string(struc_pet2, subalign_gf2->len);
			printf("PETfold RNA structure:\t\t%s %s\n",
			   include_gap_columns(substitute_char(temp_pstruc1, 'x', '.', subalign_gf1->len),
					   gap_col1, subalign_gf1->len, gapnr1, pfold_db1),
			   include_gap_columns(substitute_char(temp_pstruc2, 'x', '.', subalign_gf2->len),
					   gap_col2, subalign_gf2->len, gapnr2, pfold_db2));

			free(temp_pstruc1);
			free(temp_pstruc2);
			free(pfold_db1);
			free(pfold_db2);
		}
		if (strlen(pet_db1)+strlen(pet_db2))
			printf("PETcofold RNA structure:\t%s&%s\n", pet_db1, pet_db2);
		if (intermol_flag)
			printf("Intermolecular RNA structure:\t%s&%s\n", intermol1, intermol2);

		printf("Score_{model,structure}{tree,alignment} = %7.6f\n", score);
	}

	/* free memory */
	FreeDynEntry(T, subalign_gf1->len+subalign_gf2->len+3);

	FreeArray(paired_pet1, subalign_gf1->len);
	FreeArray(paired_pet2, subalign_gf2->len);
	free(single_pet1);
	free(single_pet2);

	FreeArray(duplex_paired_tree, subalign_gf1->len+subalign_gf2->len);
	free(duplex_single_tree);
	FreeArray(duplex_paired_energy, subalign_gf1->len+subalign_gf2->len);
	free(duplex_single_energy);
	FreeArray(duplex_paired_pet, subalign_gf1->len+subalign_gf2->len);
	free(duplex_single_pet);

	FreeAln(align1);
	FreeAln(align2);
	FreeAln(subalign_gf1);
	FreeAln(subalign_gf2);
	free(gap_col1);
	free(gap_col2);

	free(pet_db_gf);
	free(pet_db1);
	free(pet_db_gf1);
	free(pet_db2);
	free(pet_db_gf2);
	free(intermol1);
	free(intermol2);

	free(tree);
	FreePartStruc(struc_pet1);
	FreePartStruc(struc_pet2);
	FreePartStruc(ss_partial1);
	FreePartStruc(ss_partial2);
	if( setstructure1 != NULL )
		FreePartStruc(setstructure1);
	if( setstructure2 != NULL )
		FreePartStruc(setstructure2);
	if( pstruc != NULL )
		FreePartStruc(pstruc);

	return 0;
}
