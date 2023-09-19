/*
 * petfold.c
 *
 *  Created on: 06.02.2011
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
#include  "petfoldlibs.h"
#include  "thermodynamic.h"
#include  "evolutionary.h"

int main(int argc, char **argv) {
	/* PETfold parameters */
	float gap = 0.25; /* maximal allowed percent of gaps in an alignment column; deletion of column if rate is equal or higher */
	float evocon_bp = 0.9; /* minimal base pair reliability for evolutionary constraints */
	float evocon_ss = 1.0; /* minimal single stranded reliability for evolutionary constraints */
	beta = 1.; /* weighting factor for thermodynamic overlap */
	alpha = 0.2; /* weighting factor for single stranded reliabilities divided by 2 (alpha <= 0.5) */

	/* PETfold input and output */
	char *fasta = "";
	char *settree = ""; /* given phylogenetic tree in Newick tree format */
	char *setstruc = ""; /* given secondary structure */
	char *ppfile = ""; /* name of pp-file */
	static int war_flag = 0; /* fasta format output */
	static int help_flag = 0; /* help */
	verbose_flag = 0; /* long output */
	int suboptimal = 0; /* number of suboptimals to find by structure sampling */
	static int ppfold = 0; /* use pfold reliabilities calculated by parallelized ppfold */

	/* model complexity */
	setevo_flag = 1; /* usage of evolutionary reliabilities (performed by Pfold) */
	setthermo_flag = 1; /* usage of thermodynamic probabilities (performed by RNAfold) */

	Aln *align, *align_gf;
	int *gap_col, *pfold_struc;
	double *single_tree, *single_seq, *single_pet;
	double **paired_tree, **paired_seq, **paired_pet;
	char *pfold_db, *pet_db_gf, *pet_db, *evocon_db, *pet_subopt_gf,
			*pet_subopt;
	char *pfold_db_gf = NULL;
	PartStruc *evocon = NULL;
	DynEntry ***T;
	double ***P;
	int cmd, i;
	double score;
	double ensdiv;

	/* parsing long options */
	while (1)
	{
		static struct option long_options[] =
		{
				{ "setevo", no_argument, &setevo_flag, 1 },
				{ "noevo", no_argument, &setevo_flag, 0 },
				{ "setthermo", no_argument, &setthermo_flag, 1 },
				{ "nothermo", no_argument, &setthermo_flag, 0 },
				{ "war", no_argument, &war_flag, 1 },
				{ "help", no_argument, &help_flag, 1 },
				{ "verbose", no_argument, &verbose_flag, 1 },
				{ "settree", required_argument, 0, 't' },
				{ "fasta", required_argument, 0, 'f' },
				{ "setgap", required_argument, 0, 'g' },
				{ "setevocon_bp", required_argument, 0, 'p' },
				{ "setevocon_ss", required_argument, 0, 'u' },
				{ "setbeta", required_argument, 0, 'b' },
				{ "setstruc", required_argument, 0, 's' },
				{ "setalpha", required_argument, 0, 'a' },
				{ "ppfile", required_argument, 0, 'r' },
				{ "suboptimal", required_argument, 0, 'o' },
				{ "ppfold",	no_argument, &ppfold, 1 },
				{ 0, 0, 0, 0 }
		};

		/* print help text */
		if (help_flag)
			usage();

		/* getopt_long stores the option index here. */
		int option_index = 0;

		cmd = getopt_long(argc, argv, "f:g:t:p:u:b:s:a:r:o:", long_options,
				&option_index);

		/* Detect the end of the options. */
		if (cmd == -1)
			break;

		switch (cmd) {
		case 0:
			break;
		case 'f':
			fasta = optarg;
			break;
		case 'g':
			gap = atof(optarg);
			break;
		case 't':
			settree = optarg;
			break;
		case 'p':
			evocon_bp = atof(optarg);
			break;
		case 'u':
			evocon_ss = atof(optarg);
			break;
		case 'b':
			beta = atof(optarg);
			break;
		case 's':
			setstruc = optarg;
			break;
		case 'a':
			alpha = atof(optarg);
			break;
		case 'r':
			ppfile = optarg;
			break;
		case 'o':
			suboptimal = atoi(optarg);
			break;
		default:
			abort();
		}
	}

	/* read input alignment */
	if (!strlen(fasta))
		usage();
	align = get_alignment(fasta);
	//for (i=0; i < align->nr; i++) printf("%s\t%s\n", align->identifier[i], align->sequence[i]);
	/* check that all sequences have the same length -> RNA alignment */
	for (i = 1; i < align->nr; i++)
		if( strlen(align->sequence[i]) != align->len ) {
			printf("This doesn't look like a FASTA file of aligned RNA sequences - sorry.\n");
			exit(1);
		}

	if (verbose_flag)
		printf("INPUT:\n number sequences = %i\n sequence length = %i\n",
			align->nr, align->len);

	if (align->nr < 3) {
		printf("WARNING: PETfold is cancelled as your input consists of ONLY %i sequence(s).\n" ,
			align->nr);
		printf("Please provide an alignment of at least 3 sequences.\n");
		exit(1);
	}

	/* remove gap columns from alignment */
	gap_col = (int *) malloc((align->len) * sizeof(int));
	align_gf = delete_gap_columns(align, gap_col, gap);
	if (verbose_flag) {
		printf(" maximal percent of gaps per column = %3.2f\n", gap);
		printf(" gap columns = ");
		for (i = 0; i < align->len - align_gf->len; i++)
			printf("%i, ", gap_col[i]+1);
		printf("\n gap column free alignment = \n");
		for (i = 0; i < align_gf->nr; i++)
			printf("  %s\t%s\n", align_gf->identifier[i],
					align_gf->sequence[i]);
	}

	/* initialize base paired reliability matrix of evolutionary model */
	paired_tree = (double **) malloc((align_gf->len) * sizeof(double *));
	for (i = 0; i < align_gf->len; i++)
		paired_tree[i] = (double *) malloc((align_gf->len) * sizeof(double));
	/* initialize single stranded reliability array of evolutionary model */
	single_tree = (double *) malloc((align_gf->len) * sizeof(double));

	/* open tree file in Newick format */
	char *tree = (char *)calloc(BUFFERSIZE, sizeof(char *));
	if (strlen(settree)) {
		settree = read_file(settree);
		strcpy(tree, settree);
	}

	/* open structure file in dot-bracket format */
	if (strlen(setstruc))
		setstruc = read_file(setstruc);

	if (ppfold) {
		/* read bp reliabilities calculated by PPfold (multithreaded pfold) from STDIN */
		pfold_struc = get_ppfold(align_gf->len, paired_tree, single_tree);
	} else {
		/* calculate phylogenetic tree and structure reliabilities in evolutionary model -> Pfold */
		pfold_struc = get_prob_tree(align_gf, paired_tree, single_tree, tree);
	}
	//tree = (char *) realloc(tree, strlen(tree) + 1 * sizeof(char));

	//for(i=0;i<align_gf->len;i++){fprintf(stderr,"%i ",pfold_struc[i]);};fprintf(stderr,"\n");
	//for (i=0; i<align_gf->len; i++) printf("%11.9f ", single_tree[i]); printf("\n");
	//int j; for (i=1; i<=align_gf->len; i++) { for (j=1; j<=align_gf->len; j++) printf("%2d %2d %16.8e\n", i,j,paired_tree[i-1][j-1]); }
	//int j; for (i=1; i<=align_gf->len; i++) { for (j=i+1; j<=align_gf->len; j++) printf("%16.8e  ", paired_tree[i-1][j-1]); };

	if (setevo_flag) {
		/* write Pfold predicted RNA secondary structure in an string */
		pfold_db_gf = (char *) calloc(align_gf->len + 1, sizeof(char));
		get_dot_bracket(pfold_struc, align_gf->len, pfold_db_gf);
		//printf("%s\n", pfold_db_gf);

		/* evolutionary highly reliable partial structure */
		evocon = get_partial_struc(paired_tree, single_tree, pfold_struc,
				align_gf->len, evocon_bp, evocon_ss);
		//for (i=0; i<evocon->bp_nr; i++) printf("%i %i\n", evocon->bp_left[i], evocon->bp_right[i]);
		//for (i=0; i<evocon->ss_nr; i++) printf("%i\n", evocon->ss[i]);
	}

	/* write partial structure as string */
	evocon_db = partstruc2string(evocon, align_gf->len);

	/* calculate structure probabilities in thermodynamic model -> RNAfold */
	paired_seq = get_prob_paired_seq(align_gf, evocon);
	single_seq = get_prob_unpaired_seq(align_gf, paired_seq);
	//for (i=0; i<align_gf->len; i++) printf("%11.9f ", single_seq[i]); printf("\n");
	//int j; for (i=1; i<=align_gf->len; i++) { for (j=1; j<=align_gf->len; j++) printf("%2d %2d %1.15f\n", i,j,paired_seq[i-1][j-1]); }

	/* joint reliabilities of evolutionary reliabilities and thermodynamic probabilities */
	single_pet = get_pet_single(single_tree, single_seq, align_gf->len);
	paired_pet = get_pet_paired(paired_tree, paired_seq, align_gf->len);
	//for (i=0; i<align_gf->len; i++) printf("%11.9f ", single_pet[i]); printf("\n");
	//int j; for (i=0; i<align_gf->len; i++) { for (j=0; j<align_gf->len; j++) printf("%11.9f ", paired_pet[i][j]); printf("\n"); } printf("\n");

	/* parse input structure if given */
	if (strlen(setstruc))
		evocon = mod_setstruct(setstruc, gap_col, align->len,
				align->len - align_gf->len);

	/*
	 * get MEA structure using a Nussinov-style algorithm
	 * replacing products in CYK algorithm by sums
	 */
	T = nussinov(paired_pet, single_pet, evocon, align_gf->len);
	//int j; for (j=0; j<align_gf->len; j++) { printf("This is j=%i\n", j); for (i=0; i<align_gf->len; i++) printf("%11.9f %11.9f %11.9f\n", T[i][j][0]->rel, T[i][j][1]->rel, T[i][j][2]->rel); } printf("\n");

	/* backtracking to find the consensus MEA structure */
	pet_db_gf = (char *) calloc(align_gf->len + 1, sizeof(char));
	score = backtracking(T, align_gf->len, pet_db_gf);
	ensdiv = ensemble_diversity(paired_pet, align_gf->len);

	/* include gap columns */
	//pet_db = include_gap_columns(pet_db_gf, gap_col, align_gf->len, align->len-align_gf->len);
	pet_db = (char *) calloc(align->len + 1, sizeof(char));
	include_gap_columns(pet_db_gf, gap_col, align_gf->len,
			align->len - align_gf->len, pet_db);

	/* write pp-file */
	if (strlen(ppfile))
		create_ppfile(paired_pet, single_pet, gap_col, align_gf->len, align->len - align_gf->len, ppfile);
		//create_ppfile(paired_pet, single_pet, ppfile, align_gf->len);

	/* output */
	if (war_flag) {
		for (i = 0; i < align->nr; i++)
			printf(">%s\n%s\n", align->identifier[i], align->sequence[i]);
		printf(">structure\n%s\n", pet_db);
	} else {
		if (setevo_flag) {
			pfold_db = (char *) calloc(align->len + 1, sizeof(char));
			printf("Pfold RNA structure:\t%s\n",
					include_gap_columns(pfold_db_gf, gap_col, align_gf->len,
							align->len - align_gf->len, pfold_db));
			printf("Constraints:\t\t%s\n",
					include_gap_columns(evocon_db, gap_col, align_gf->len,
							align->len - align_gf->len, pfold_db));

			free(pfold_db);
			free(pfold_db_gf);
		}
		if (strlen(pet_db))
			printf("PETfold RNA structure:\t%s\n", pet_db);
		printf("Score_{model,structure}{tree,alignment} = %7.6f\n", score);
		printf("Length-normalized ensemble diversity = %7.6f\n", ensdiv);
	}

	/* structure sampling to find suboptimal structures */
	if (suboptimal) {
		int j;

		P = nussinov_subopt(paired_pet, single_pet, align_gf->len);
		//for (j=0; j<align_gf->len; j++) { printf("This is j=%i\n", j); for (i=0; i<align_gf->len; i++) printf("%11.9f %11.9f %11.9f\n", P[i][j][0], P[i][j][1], P[i][j][2]); } printf("\n");

		//float score;
		srand(time(NULL ));
		for (i = 0; i < suboptimal; i++) {
			pet_subopt_gf = (char *) calloc(align_gf->len + 1, sizeof(char));
			pet_subopt = (char *) calloc(align->len + 1, sizeof(char));

			//backtracking_sampling(P, paired_pet, align_gf->len, &score, pet_subopt_gf);
			score = backtracking_subopt(P, paired_pet, single_pet,
					align_gf->len, pet_subopt_gf);
			printf("Suboptimal structure:   %s\t%7.6f\n",
					include_gap_columns(pet_subopt_gf, gap_col, align_gf->len,
							align->len - align_gf->len, pet_subopt), score);

			free(pet_subopt_gf);
			free(pet_subopt);
		}

		for (i = 0; i < align_gf->len; i++)
			for (j = 0; j < align_gf->len; j++)
				free(P[i][j]);
		free(P);
	}

	free(pet_db_gf);
	free(pet_db);

	free(gap_col);
	free(pfold_struc);

	if (setevo_flag)
		FreePartStruc(evocon);

  free(tree);
	free(evocon_db);

	FreeDynEntry(T, align_gf->len);

	FreeArray(paired_pet, align_gf->len);
	FreeArray(paired_seq, align_gf->len);
	FreeArray(paired_tree, align_gf->len);

	FreeAln(align);
	FreeAln(align_gf);

	free(single_tree);
	free(single_seq);
	free(single_pet);

	return 0;
}
