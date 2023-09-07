/*
 * petfoldlibs.h
 *
 *  Created on: 06.02.2011
 *      Author: Stefan Seemann
 */

#ifndef PETFOLDLIBS_H_
#define PETFOLDLIBS_H_

#ifdef MODUL0
    #define EXTERN
  #else
    #define EXTERN extern
  #endif

#define INFINITE_D 1000000000.	/* value of undefined double fields */
#define INFINITE_I 1000000000	/* value of undefined integer fields */

EXTERN int setevo_flag;		/* usage of evolutionary reliabilities (performed by Pfold) */
EXTERN int setthermo_flag;	/* usage of thermodynamic probabilities (performed by RNAfold) */
EXTERN float beta;   		/* weighting factor for thermodynamic overlap */
EXTERN float alpha;			/* weighting factor for single stranded reliabilities (alpha <= 0.5) */
EXTERN int verbose_flag;	/* long output */
EXTERN int subopt_flag;		/* structure sampling to find suboptimal structures */

#define round(x) ((x)>=0?(unsigned long)((x)+1.0):(unsigned long)((x)-1.0))    /* round up to next larger integer */

typedef struct tagAln {
	char **identifier;      /* Text in entry */
	char **sequence;        /* Sequence info */
	int nr;					/* Number of sequences */
	int len;             	/* Sequence length */
} Aln;

typedef struct tagSeqList {
	char **sequence;		/* Sequence info */
	int nr;					/* Number of sequences */
	int *len;				/* Sequence length */
	int **origid;			/* Each index is attributed an original index,
							 * e.g, from the sequence with gaps */
} SeqList;

typedef struct tagPartStruc {
	int *bp_left;			/* left base of constrained base pairs */
	int *bp_right;			/* right base of constrained base pairs */
	int *ss;				/* constrained unpaired base */
	int bp_nr;				/* number of constrained base pairs */
	int ss_nr;				/* number of constrained unpaired bases */
} PartStruc;

struct tagDynEntry {
	double rel;				/* reliability */
	struct tagDynEntry *lchild;		/* reference to left child rule or infinite for rule L->s */
	struct tagDynEntry *rchild;		/* reference to right child rule or infinite */
};
typedef struct tagDynEntry* DynEntry;

/* stochastic context-free grammar (scfg) rule probabilities taken from Knudsen et al. (2003) */
typedef struct tagSCFG {
	double S_LS;   /* rule S->LS */
	double S_L;    /* rule S->L */
	double F_dFd;  /* rule F->dFd*/
	double F_LS;   /* rule F->LS */
	double L_s;    /* rule L->s */
	double L_dFd;  /* rule L->dFd */
} scfgProbs;


Aln *get_alignment(char *fasta);
char *read_file(char *file);
Aln *delete_gap_columns(Aln *align, int *gap_col, float gap);
char *get_dot_bracket(int *struct_coord, int len, char *db);
double adjustprecision(double val);
PartStruc *get_partial_struc(double **paired_tree, double *single_tree,
		int *pfold_struct, int len, float evocon_bp, float evocon_ss);
char *get_constraint_string(PartStruc *evocon, int *origid, int len_gf);
double *get_pet_single(double *single_tree, double *single_seq, int len);
double **get_pet_paired(double **paired_tree, double **paired_seq, int len);
PartStruc *mod_setstruct(char *setstruc, int *gap_col, int len, int gapnr);
DynEntry ***nussinov(double **paired_pet, double *single_pet, PartStruc *evocon, int len);
int *partstruc2list(PartStruc *evocon, int len);
int *partstruc2str(PartStruc *evocon, int len);
char *partstruc2string(PartStruc *evocon, int len);
PartStruc *str2partstruc(int *struc, int len);
double add(double x, double y);
_Bool max(double x, double y);
double backtracking(DynEntry ***T, int len, char *db);
char *include_gap_columns(char *db, int *gap_col, int len, int gapnr, char *newdb);
int create_ppfile(double **paired, double *single, int *gap_col, int len, int gapnr, char *ppfile);
double ensemble_diversity(double **paired, int len);
double ***nussinov_subopt(double **paired_pet, double *single_pet, int len);
double backtracking_subopt(double ***P, double **paired_pet, double *single_pet, int len, char *pet_subopt_gf);
int statistical_sampling(unsigned long *sam_prob, int sam_len);
void change_list_sum_to_max(unsigned long *sam_prob, int start, int end);
void chomp(char *s);
void usage();

void InitPartStruc(PartStruc *evocon);

void FreeAln(Aln *align);
void FreeSeqList(SeqList *seqlist);
void FreePartStruc(PartStruc *evocon);
void FreeDynEntry(DynEntry ***T, int len);
void FreeArray(double **array, int len);

#endif /* PETFOLDLIBS_H_ */
