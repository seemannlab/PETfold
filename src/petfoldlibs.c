/*
 * petfoldlibs.c
 *
 *  Created on: 06.02.2011
 *      Author: Stefan Seemann
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "petfoldlibs.h"
#include "file.h"

/*
 * extract the alignment from FASTA file in an Aln object
 */
Aln *get_alignment(char *fasta)
{
	char *s;         /* String for keeping lines */
	int ptr = 0;         /* Pointer to line number */
	int maxlen = 30;      /* Start with 30 lines - Dynamic line allocation maximum */
	Aln *align;
	FILE *fp;
	int d;
	int seqlen = 0;

	/* Initialization of Aln Object */
	align = (Aln *)malloc(sizeof(Aln));
	align->identifier = (char **)malloc(maxlen * sizeof(char *));
	align->sequence = (char **)malloc(maxlen * sizeof(char *));

	/* open alignment file in FASTA format */
	if( (fp = fopen(fasta,"rb")) == NULL )
	{
		fprintf(stderr, "Can't open the file '%s'\n", fasta);
		exit(1);
	}
	ptr = 0;

	while ( (s = GetLine(fp)) != NULL )
	{
		chomp(s);  // replace newline by string terminator
		if (strncmp(s, ">", 1) == 0)
		{
			ptr++;
			if (ptr == maxlen) {
				maxlen *= 2;
				align->identifier = (char **)realloc(align->identifier, maxlen * sizeof(char *));
				align->sequence = (char **)realloc(align->sequence, maxlen * sizeof(char *));
																			/* Get more space */
			}
			align->identifier[ptr-1] = (char *)malloc((strlen(s)+1) * sizeof(char));

			strcpy(align->identifier[ptr-1],s+1);
			//strcpy(align->identifier[ptr-1],"n");
			//strcat(align->identifier[ptr-1],s+1);

			/* replace spaces in sequence identifiers by underscore */
			unsigned char *tmp = s;
			s = strchr(align->identifier[ptr-1],' ');
			while( s != NULL )
			{
				*s = '_';
				s = strchr(align->identifier[ptr-1],' ');
			}
			s = tmp;

			seqlen = 0;
		}
		else {
            for( d=0; d<strlen(s); d++ )
            {
            	/* convert sequence to uppercase for hairpin loop lookups in RNAfold library */
                s[d] = toupper(s[d]);
                /* create RNA sequence T->U */
                if (s[d] == 'T') s[d] = 'U';
            }

            if( !seqlen )
            	align->sequence[ptr-1] = (char *)malloc((strlen(s)+1) * sizeof(char));
            else
            	align->sequence[ptr-1] = (char *)realloc(align->sequence[ptr-1], (seqlen + strlen(s) + 1) * sizeof(char));
            strcpy(align->sequence[ptr-1]+seqlen,s);
           	seqlen += strlen(s);
		}
		free(s);
	}
	fclose(fp);

	align->identifier = (char **)realloc(align->identifier, ptr * sizeof(char *));
	align->sequence = (char **)realloc(align->sequence, ptr * sizeof(char *));
	                                             /* Free excess space */
	align->nr = ptr;
	align->len = strlen(align->sequence[0]);

	return align;
}

/*
 * open and read file or copy string
 */
char *read_file(char *file) {
	char *s;
	FILE *fp;


	if( (fp = fopen(file,"rb")) != NULL )
	{
		if( (s = GetLine(fp)) != NULL )
		{
			chomp(s);  // replace newline by string terminator
		}
		else {
			s = file;
		}
	}
	else {
		s = file;
	}

	return s;
}

/*
 * delete columns in an alignment with more or equal than 'gap'x100% gaps (default 25% in Pfold)
 */
Aln *delete_gap_columns(Aln *align, int *gap_col, float gap)
{
	Aln *align_gf;
	int i, j;
	int k[align->len];

	/* find columns where gap rate is higher and equal as 'gap'x100% */
	int m = 0;
	int n = 0;
	for(j=0; j < align->len; j++)
	{
		int nr = 0;
		for(i=0; i < align->nr; i++)
		{
			if( align->sequence[i][j] == '-' || align->sequence[i][j] == '.')
				nr++;
		}
		if( (double) nr/(align->nr) >= gap )
			gap_col[m++] = j;	/* gap column */
		else
			k[n++] = j;	/* non-gap column */
	}

	/* Initialization of Aln Object */
	align_gf = (Aln *)malloc(sizeof(Aln));
	align_gf->identifier = (char **)malloc(align->nr * sizeof(char *));
	align_gf->sequence = (char **)malloc(align->nr * sizeof(char *));

	align_gf->len = n;
	int offset = 0;
	/* write Aln object without gaps */
	for(i=0; i < align->nr; i++)
	{
		/* check that sequence without gap-free columns has still valid nucleotides */
		int nr = 0;
		for(j=0; j < align_gf->len; j++) {
			if( align->sequence[i][k[j]] == '-' )
				nr++;
		}
		if( nr == align_gf->len ) {
			printf("Sequence %s is removed because it has only gaps in gap-column-free alignment.\n", align->identifier[i]);
			offset++;
			continue;
		}

		align_gf->identifier[i-offset] = (char *)malloc((strlen(align->identifier[i])+1) * sizeof(char));
		strcpy(align_gf->identifier[i-offset],align->identifier[i]);

		align_gf->sequence[i-offset] = (char *)malloc((align_gf->len+1) * sizeof(char));
		for(j=0; j < align_gf->len; j++) {
			align_gf->sequence[i-offset][j] = align->sequence[i][k[j]];
		}
		align_gf->sequence[i-offset][j] = '\0';
	}
	align_gf->nr = align->nr-offset;
	if( offset ) {
		align_gf->identifier = (char **)realloc(align_gf->identifier, align_gf->nr * sizeof(char *));
		align_gf->sequence = (char **)realloc(align_gf->sequence, align_gf->nr * sizeof(char *));
	}

	return align_gf;
}


/*
 * convert structure coordinates in dot-bracket string
 */
char *get_dot_bracket(int *struct_coord, int len, char *dotb)
{
	int i;

	for (i=0; i<len; i++)
	{
		if (struct_coord[i] == -1)
			dotb[i] = '.';
		else
			if (struct_coord[i] > i)
				dotb[i] = '(';
			else
				dotb[i] = ')';
	}
	dotb[i] = '\0';

	return dotb;
}


/*
 * test function to adjust precision of doubles
 */
double adjustprecision(double val)
{
	float tmp = val;
	double newval = tmp;

	fprintf(stderr,"%f %f %20.19f\n",val,tmp,newval);
	return newval;
}


/*
 * get indices of base pairs and unpaired bases in highly reliable partial structure
 */
PartStruc *get_partial_struc(double **paired_tree, double *single_tree,
		int *pfold_struct, int len, float evocon_bp, float evocon_ss)
{
	PartStruc *pstruc;
	int i;

	/* Initialization of PartStruc Object */
	pstruc = (PartStruc *)malloc(sizeof(PartStruc));
	pstruc->bp_left = (int *)malloc(len/2 * sizeof(int));
	pstruc->bp_right = (int *)malloc(len/2 * sizeof(int));
	pstruc->ss = (int *)malloc(len * sizeof(int));

	pstruc->bp_nr = 0;
	pstruc->ss_nr = 0;
	for (i=0; i<len; i++)
	{
		if( pfold_struct[i]==-1 )
		{
			if( (float)single_tree[i]>evocon_ss )
			{
				//printf("%10.9e %10.9e\n", single_tree[i], evocon_ss);
				pstruc->ss[pstruc->ss_nr++] = i;
			}
		}
		else
			if( pfold_struct[i]<i )
				continue;
			else
			{
				//printf("%i\t%i\t%5.4f\n",i,pfold_struct[i],paired_tree[i][pfold_struct[i]]);
				if( (float)paired_tree[i][pfold_struct[i]]>evocon_bp )
				{
					pstruc->bp_left[pstruc->bp_nr] = i;
					pstruc->bp_right[pstruc->bp_nr++] = pfold_struct[i];
				}
			}
	}

	/* Free excess space */
	pstruc->bp_left = (int *)realloc(pstruc->bp_left, pstruc->bp_nr * sizeof(int));
	pstruc->bp_right = (int *)realloc(pstruc->bp_right, pstruc->bp_nr * sizeof(int));
	pstruc->ss = (int *)realloc(pstruc->ss, pstruc->ss_nr * sizeof(int));

	return pstruc;
}


/*
 * returns a constraint string usable with 'RNAfold -C'
 */
char *get_constraint_string(PartStruc *evocon, int *origid, int len_gf)
{
	char *dotb;
	int i, k, l, m;

	dotb = (char *) malloc (len_gf+1 * sizeof(char));

	for (i=0; i<len_gf; i++)
		dotb[i] = '.';
	dotb[i] = '\0';

	if ( !setevo_flag )
		return dotb;

	k=0;
	for (i=0; i<evocon->bp_nr; i++)
	{
		for (l=k; l<len_gf; l++)
		{
			if( origid[l]==evocon->bp_left[i] )
				for (m=l+1; m<len_gf; m++)
				{
					if( origid[m]==evocon->bp_right[i])
					{
						dotb[l] = '(';
						dotb[m] = ')';
						k=l+1;
						break;
					}
				}
			else
				if( origid[l]>evocon->bp_left[i] )
					break;
		}
	}

	k=0;
	for (i=0; i<evocon->ss_nr; i++)
		for (l=k; l<len_gf; l++)
			if( origid[l]==evocon->ss[i] )
			{
				dotb[l] = 'x';
				k=l+1;
				break;
			}
			else
				if( origid[l]>evocon->ss[i] )
					break;

	return dotb;
}


/*
 * calculates single stranded reliabilities of the PET-model
 */
double *get_pet_single(double *single_tree, double *single_seq, int len)
{
	double *single;
	int i;

	/* initialize single stranded reliability array of PET-model */
	single = (double *)calloc(len, sizeof(double));

	for (i=0; i<len; i++)
		single[i] = (single_tree[i] + beta * single_seq[i]) / 2;

	return single;
}


/*
 * calculates base paired reliabilities of the PET-model
 */
double **get_pet_paired(double **paired_tree, double **paired_seq, int len)
{
	double **paired;
	int i, j;

	/* initialize base paired reliability matrix of PET-model */
	paired = (double **)calloc(len, sizeof(double *));
	for (i=0; i<len; i++)
		paired[i] = (double *)calloc(len, sizeof(double));

	for (i=0; i<len; i++)
		for (j=0; j<len; j++)
			paired[i][j] = INFINITE_D;

	for (j=1; j<len; j++)
		for (i=0; i<len-j; i++)
			paired[i][i+j] = (paired_tree[i][i+j] + beta * paired_seq[i][i+j]) / 2;

	return paired;
}


/*
 * modifies the desired structure to partial constrained structure and
 * tests the accordance of desired structure to the gap columns (gap must be an unpaired position!!!)
 */
PartStruc *mod_setstruct(char *setstruc, int *gap_col, int len, int gapnr)
{
	PartStruc *pstruc;
	char *tmpstruc;
	int *bpstack;
	int i, offset=0, k=0;

	/* structure must have the same length as alignment length */
	if( strlen(setstruc)!=len )
	{
		fprintf(stderr, "Error <setstruc>: "
				"Structure length (%i) is different from alignment (%i).\n", (int) strlen(setstruc), len);
		exit(1);
	}

    /* substitute alternative symbols to standard dot-bracket notation */
    for (i=0; i<len; i++)
    {
        if (setstruc[i] == '<') setstruc[i] = '(';
        else if (setstruc[i] == '>') setstruc[i] = ')';
        else if (setstruc[i] == '-') setstruc[i] = '.';
        else if (setstruc[i] == '[') setstruc[i] = '(';
        else if (setstruc[i] == ']') setstruc[i] = ')';
    }

	/*
	 * remove alignment gap columns from the input structure
	 * if a gap column is not unpaired in the input structure then exit
	 */
	tmpstruc = (char *)malloc((len-gapnr) * sizeof(char));
	for (i=0; i<len; i++)
	{
		if( gap_col[offset] == i )
		{
			if( setstruc[gap_col[offset]] != '.')
			{
				fprintf(stderr, "Error <setstruc>: "
						"Position %i must be unpaired due to an alignment gap.\n", gap_col[offset]+1);
				exit(1);
			}
			offset++;
		}
		else
		{
			tmpstruc[i-offset] = setstruc[i];
		}
	}

	/* Initialization of PartStruc Object */
	pstruc = (PartStruc *)malloc(sizeof(PartStruc));
	pstruc->bp_left = (int *)malloc(len/2 * sizeof(int));
	pstruc->bp_right = (int *)malloc(len/2 * sizeof(int));
	pstruc->ss = (int *)malloc(len * sizeof(int));

	bpstack = (int *)malloc(len/2 * sizeof(int));

	pstruc->bp_nr = 0;
	pstruc->ss_nr = 0;
	for (i=0; i<len-gapnr; i++)
	{
		switch(tmpstruc[i])
		{
			case '.':
				pstruc->ss[pstruc->ss_nr++] = i;
				break;
			case '(':
				bpstack[k++] = pstruc->bp_nr;
				pstruc->bp_left[pstruc->bp_nr++] = i;
				break;
			case ')':
				pstruc->bp_right[bpstack[--k]] = i; /* take last open base from stack */
				break;
			default:
				printf("%i "
						"%c\n",i,tmpstruc[i]);
				fprintf(stderr, "Error <setstruc>: "
						"The only allowed characters are '(', ')', '<', '>', '[', ']', '.', and '-'.\n");
				exit(1);
		}
	}

	if( k )
	{
		fprintf(stderr, "Error <setstruc>: "
				"Number of left and right parenthesis are unequal!\n");
		exit(1);
	}

	/* Free excess space */
	pstruc->bp_left = (int *)realloc(pstruc->bp_left, pstruc->bp_nr * sizeof(int));
	pstruc->bp_right = (int *)realloc(pstruc->bp_right, pstruc->bp_nr * sizeof(int));
	pstruc->ss = (int *)realloc(pstruc->ss, pstruc->ss_nr * sizeof(int));

	free(bpstack);
	free(tmpstruc);

	return pstruc;
}


/*
 * get maximum expected accuracy (MEA) structure using a Nussinov-style recursion
 * replacing products in CYK algorithm by sums
 */
DynEntry ***nussinov(double **paired_pet, double *single_pet, PartStruc *evocon, int len)
{
	DynEntry ***T;
	int *consstr;
	double rel;
	int i, j, k, t;

	/* convert partial structure in list of integers */
	consstr = partstruc2list(evocon, len);

	/* initialize cube T which holds arrays
	 * column(i=1,..,len+1-j): base position;
	 * row(j=1,...,len): length of subsequence;
	 * depth(t): non terminals of SCFG: [0]=>'S', [1]=>'F', [2]=>'L'
	 * DynEntry: rel=>reliability,
	 *           lchild=>reference to left child rule or infinite for rule L->s,
	 *           rchild=>reference to right child rule or infinite
	*/
	T = (DynEntry***)malloc(len * sizeof(DynEntry**));
	for( i=0; i<len; i++ ) {
		T[i] = (DynEntry**)malloc(len * sizeof(DynEntry*));
		for( j=0; j<len; j++ ) {
			T[i][j] = (DynEntry*)malloc(3 * sizeof(DynEntry));
			for( t=0; t<3; t++ )
			{
				/* Initialization of DynEntry Object */
				T[i][j][t] = (DynEntry)malloc(sizeof(struct tagDynEntry));
				T[i][j][t]->rel = INFINITE_D;
				T[i][j][t]->lchild = NULL;
				T[i][j][t]->rchild = NULL;
			}
		}
	}

	/* write first row of cube T with reliabilities of unpaired bases in alignment */
	for( i=0; i<len; i++ )
	{
		if (consstr[i]!=INFINITE_I && consstr[i]!=-1)
			continue;

		/* L->s */
	    T[i][0][2]->rel = 2*alpha * single_pet[i]/len;

	    /* S->L */
	    T[i][0][0]->rel = T[i][0][2]->rel;
	    T[i][0][0]->lchild = T[i][0][2];
	}

	/* fill cube T */
	for( j=0; j<len; j++ )
	{
		for( i=0; i<len-j; i++ )
		{
			/* rules S->LS and F->LS (standard CYK algorithm) */
			for( k=1; k<j+1; k++ )
			{
				if (T[i][k-1][2]->rel==INFINITE_D || T[i+k][j-k][0]->rel==INFINITE_D)
					continue;

				/* S->LS: T[i,j][S] := max{ T[i,j][S], T[i,k-1][L] x T[i+k,j-k][S] } */
				rel = add(T[i][k-1][2]->rel, T[i+k][j-k][0]->rel);
				if (max(rel, T[i][j][0]->rel))
				{
					T[i][j][0]->rel = rel;
					T[i][j][0]->lchild = T[i][k-1][2];
					T[i][j][0]->rchild = T[i+k][j-k][0];
				}
				/* F->LS: T[i,j][F] := max{ T[i,j][F], T[i,k-1][L] x T[i+k,j-k][S] } */
				rel = add(T[i][k-1][2]->rel, T[i+k][j-k][0]->rel);
				if (max(rel, T[i][j][1]->rel))
				{
					T[i][j][1]->rel = rel;
					T[i][j][1]->lchild = T[i][k-1][2];
					T[i][j][1]->rchild = T[i+k][j-k][0];
				}
			}

			/*
			 * rule F->dFd (inner bond)
			 * loop consists of at least 3 unpaired bases
			 */
			if (j>3)
			{
				if ( (consstr[i]==INFINITE_I && consstr[i+j]==INFINITE_I) ||
					 (T[i+1][j-2][1]->rel!=INFINITE_D && consstr[i]==consstr[i+j] && consstr[i]!=-1) )
				{
					/* F->dFd: T[i,j][F] := max{ T[i,j][F], T[i+1,j-2][F] x Rel_paired[i,i+j] } */
					rel = add(T[i+1][j-2][1]->rel, 2 * paired_pet[i][i+j]/len);
					if (max(rel, T[i][j][1]->rel))
					{
						T[i][j][1]->rel = rel;
						T[i][j][1]->lchild = T[i+1][j-2][1];
						T[i][j][1]->rchild = NULL;
					}
				}
			}

			/*
			 * rule L->dFd (open bond)
			 * loop consists of at least 3 unpaired bases
			 */
			if (j>3)
			{
				if ( (consstr[i]==INFINITE_I && consstr[i+j]==INFINITE_I) ||
					 (T[i+1][j-2][1]->rel!=INFINITE_D && consstr[i]==consstr[i+j] && consstr[i]!=-1) )
				{
					/* L->dFd: T[i,j][L] := max{ T[i,j][L], T[i+1][j-2][F] x Rel_paired[i,i+j] } */
					rel = add(T[i+1][j-2][1]->rel, 2 * paired_pet[i][i+j]/len);
					if (max(rel, T[i][j][2]->rel))
					{
						T[i][j][2]->rel = rel;
						T[i][j][2]->lchild = T[i+1][j-2][1];
						T[i][j][2]->rchild = NULL;
					}
				}
			}

			/* S->L: T[i,j][S] := max{ T[i,j][S], T[i,j][L] } */
			if (T[i][j][2]->rel!=INFINITE_D)
			{
				rel = T[i][j][2]->rel;
				if (max(rel, T[i][j][0]->rel))
				{
					T[i][j][0]->rel = rel;
					T[i][j][0]->lchild = T[i][j][2];
					T[i][j][0]->rchild = NULL;
				}
			}
		}
	}

	free(consstr);

	if (T[0][len-1][0]->rel==INFINITE_D)
	{
		fprintf(stderr, "The consensus RNA structure could not be predicted!\n");
		exit(1);
	}

	return T;
}


/*
 * convert PartStruc Object in a list of integers:
 * sign paired bases by incremental index,
 * unpaired bases as '-1' and unconstrained as INFINITIVE
 */
int *partstruc2list(PartStruc *evocon, int len)
{
	int *consstr;
	int i;

	/* Initialization */
	consstr = (int *)malloc(len * sizeof(int));
	for (i=0; i<len; i++)
		consstr[i] = INFINITE_I;

	if (setevo_flag) {
		for (i=0; i<evocon->ss_nr; i++)
			consstr[evocon->ss[i]] = -1;

		for (i=0; i<evocon->bp_nr; i++)
		{
			consstr[evocon->bp_left[i]] = i;
			consstr[evocon->bp_right[i]] = i;
		}
	}

	return consstr;
}


/*
 * convert PartStruc Object in a list of integers:
 * sign paired bases by index of their pairing partner,
 * unpaired bases as '-1' and unconstrained as INFINITIVE
 */
int *partstruc2str(PartStruc *evocon, int len) {
	int *consstr;
	int i;

	/* Initialization */
	consstr = (int *)malloc(len * sizeof(int));
	for (i=0; i<len; i++)
		consstr[i] = INFINITE_I;

	if (setevo_flag) {
		for (i=0; i<evocon->ss_nr; i++)
			consstr[evocon->ss[i]] = -1;

		for (i=0; i<evocon->bp_nr; i++)
		{
			consstr[evocon->bp_left[i]] = evocon->bp_right[i];
			consstr[evocon->bp_right[i]] = evocon->bp_left[i];
		}
	}

	return consstr;
}


/*
 * convert PartStruc Object in a string
 */
char *partstruc2string(PartStruc *evocon, int len)
{
	char *dotb;
	int i;

	/* Initialization */
	dotb = (char *)calloc(len+1, sizeof(char));
	for (i=0; i<len; i++)
		dotb[i] = '.';
	dotb[i] = '\0';

  if (setevo_flag) {
  	for (i=0; i<evocon->bp_nr; i++)
		{
			dotb[evocon->bp_left[i]] = '(';
			dotb[evocon->bp_right[i]] = ')';
		}

		for (i=0; i<evocon->ss_nr; i++)
			dotb[evocon->ss[i]] = 'x';
	}

	return dotb;
}


/*
 * convert a list of integers in PartStruc Object:
 * sign paired bases by index of their pairing partner,
 * unpaired bases as '-1' and unconstrained as INFINITIVE
 */
PartStruc *str2partstruc(int *struc, int len) {
	PartStruc *pstruc;
	int i;

	/* Initialization of PartStruc Object */
	pstruc = (PartStruc *)malloc(sizeof(PartStruc));
	pstruc->bp_left = (int *)malloc(len/2 * sizeof(int));
	pstruc->bp_right = (int *)malloc(len/2 * sizeof(int));
	pstruc->ss = (int *)malloc(len * sizeof(int));

	pstruc->bp_nr = 0;
	pstruc->ss_nr = 0;

	for (i=0; i<len; i++) {
		if( struc[i] == -1 )
			pstruc->ss[pstruc->ss_nr++] = i;
		else
			if( struc[i] != -1 && struc[i] != INFINITE_I && i < struc[i] ) {
				pstruc->bp_left[pstruc->bp_nr] = i;
				pstruc->bp_right[pstruc->bp_nr++] = struc[i];
			}
	}

	/* Free excess space */
	pstruc->bp_left = (int *)realloc(pstruc->bp_left, pstruc->bp_nr * sizeof(int));
	pstruc->bp_right = (int *)realloc(pstruc->bp_right, pstruc->bp_nr * sizeof(int));
	pstruc->ss = (int *)realloc(pstruc->ss, pstruc->ss_nr * sizeof(int));

	return pstruc;
}


/*
 * add two values considering infinite values
 */
double add(double x, double y)
{
	if (x!=INFINITE_D && y!=INFINITE_D)
		return x + y;
	else
		return INFINITE_D;
}


/*
 * return maximum of two values considering infinite values
 * true if first value is maximum, else false
 */
_Bool max(double x, double y)
{
	if (x!=INFINITE_D && y!=INFINITE_D)
		if (x>=y)
			return 1;
		else
			return 0;
	else
		if (y==INFINITE_D)
			return 1;
		else
			return 0;
}


/*
 * backtrack through the CYK-cube from T[0,n-1][S] to T[0..n-1,0][L] where n=alignment length
 * get the consensus secondary structure of the alignment
 */
double backtracking(DynEntry ***T, int len, char *struc)
{
	int *rid;   /* stack of non-terminals: S->0, F->1, L->2 */
	int *depth;   /* stack counting number of building rules for each non-terminals */
	//DynEntry *c1;   /* stack of left children of non-terminals */
	//DynEntry *c2;   /* stack of right children of non-terminals */
	DynEntry **c1;   /* stack of left children of non-terminals */
	DynEntry **c2;   /* stack of right children of non-terminals */
	int s_len, r_len, d_len, c1_len, c2_len;
	int i, offset;

	/* Initialization */
	rid = (int *)malloc(3*len * sizeof(int));
	int tmp_len = 2;
	if (len/5 > tmp_len) {
		tmp_len = len/5;
	}
	depth = (int *)malloc(tmp_len * sizeof(int));
	//c1 = (DynEntry *)malloc(len * sizeof(DynEntry));
	//c2 = (DynEntry *)malloc(len * sizeof(DynEntry));
	c1 = (DynEntry **)calloc(2*len, sizeof(DynEntry *));
	c2 = (DynEntry **)calloc(2*len, sizeof(DynEntry *));

	s_len = r_len = d_len = c1_len = c2_len = 0;
	struc[s_len++] = 'S';
	struc[len] = '\0';
	rid[r_len++] = 0;
	depth[d_len++] = 0;

	/* get children of T[0,n-1][S] */
	c1[c1_len++] = &T[0][len-1][0]->lchild;
	c2[c2_len++] = &T[0][len-1][0]->rchild;

	while(1)
	{
		//printf("%s %i %i %i %i %i\n", struc, s_len, r_len, c1_len, c2_len, d_len);
		//for(i=0; i<d_len; i++) printf("%i ", depth[i]);printf("\n");
		//for(i=0; i<r_len; i++) printf("%i ", rid[i]);printf("\n");

		switch( rid[r_len-1] )
		{
			/* actual non-terminal is S */
			case 0:
				/* rule S->LS */
				//if (c2[c2_len-1]!=NULL)
				if (*c2[c2_len-1]!=NULL)
				{
					s_len++;
					offset = strchr(struc, 'S') - struc;
					for( i=s_len; i>offset+1; i-- )
						struc[i] = struc[i-1];
					struc[offset] = 'L';
					struc[offset+1] = 'S';

					/* S rule */
					rid[r_len++] = 0;
					//c1[c1_len++] = c2[c2_len-1]->lchild;
					//c2[c2_len] = c2[c2_len-1]->rchild;
					c1[c1_len++] = &(*c2[c2_len-1])->lchild;
					c2[c2_len] = &(*c2[c2_len-1])->rchild;
					c2_len++;

					/* L rule */
					rid[r_len++] = 2;
					//c2[c2_len++] = c1[c1_len-2]->rchild;
					//c1[c1_len] = c1[c1_len-2]->lchild;
					c2[c2_len++] = &(*c1[c1_len-2])->rchild;
					c1[c1_len] = &(*c1[c1_len-2])->lchild;
					c1_len++;
					depth[d_len-1]++;
					depth[d_len++] = 1;
				}
				/* rule S->L */
				else
				{
					offset = strchr(struc, 'S') - struc;
					struc[offset] = 'L';

					rid[r_len++] = 2;
					//c2[c2_len++] = c1[c1_len-1]->rchild;
					//c1[c1_len] = c1[c1_len-1]->lchild;
					c2[c2_len++] = &(*c1[c1_len-1])->rchild;
					c1[c1_len] = &(*c1[c1_len-1])->lchild;
					c1_len++;
					depth[d_len-1]++;
				}
				break;
			/* actual non-terminal is F */
			case 1:
				/* rule F->LS */
				//if (c2[c2_len-1]!=NULL)
				if (*c2[c2_len-1]!=NULL)
				{
					s_len++;
					offset = strchr(struc, 'F') - struc;
					for( i=s_len; i>offset+1; i--)
						struc[i] = struc[i-1];
					struc[offset] = 'L';
					struc[offset+1] = 'S';

					/* S rule */
					rid[r_len++] = 0;
					//c1[c1_len++] = c2[c2_len-1]->lchild;
					//c2[c2_len] = c2[c2_len-1]->rchild;
					c1[c1_len++] = &(*c2[c2_len-1])->lchild;
					c2[c2_len] = &(*c2[c2_len-1])->rchild;
					c2_len++;

					/* L rule */
					rid[r_len++] = 2;
					//c2[c2_len++] = c1[c1_len-2]->rchild;
					//c1[c1_len] = c1[c1_len-2]->lchild;
					c2[c2_len++] = &(*c1[c1_len-2])->rchild;
					c1[c1_len] = &(*c1[c1_len-2])->lchild;
					c1_len++;
					depth[d_len-1]++;
					depth[d_len++] = 1;
				}
				/* rule F->dFd */
				else
				{
					s_len += 2;
					offset = strchr(struc, 'F') - struc;
					for( i=s_len; i>offset+2; i--)
						struc[i] = struc[i-2];
					struc[offset] = '(';
					struc[offset+1] = 'F';
					struc[offset+2] = ')';

					rid[r_len++] = 1;
					//c2[c2_len++] = c1[c1_len-1]->rchild;
					//c1[c1_len] = c1[c1_len-1]->lchild;
					c2[c2_len++] = &(*c1[c1_len-1])->rchild;
					c1[c1_len] = &(*c1[c1_len-1])->lchild;
					c1_len++;
					depth[d_len-1]++;
				}
				break;
			/* actual non-terminal is L */
			case 2:
				/* rule L->dFd */
				//if (c1[c1_len-1]!=NULL)
				if (*c1[c1_len-1]!=NULL)
				{
					s_len += 2;
					offset = strchr(struc, 'L') - struc;
					for( i=s_len; i>offset+2; i--)
						struc[i] = struc[i-2];
					struc[offset] = '(';
					struc[offset+1] = 'F';
					struc[offset+2] = ')';

					rid[r_len++] = 1;
					//c2[c2_len++] = c1[c1_len-1]->rchild;
					//c1[c1_len] = c1[c1_len-1]->lchild;
					c2[c2_len++] = &(*c1[c1_len-1])->rchild;
					c1[c1_len] = &(*c1[c1_len-1])->lchild;
					c1_len++;
					depth[d_len-1]++;
				}
				/* rule L->s */
				else
				{
					offset = strchr(struc, 'L') - struc;
					struc[offset] = '.';

					for( i=0; i<depth[d_len-1]; i++ )
					{
						c1_len--;
						c2_len--;
						r_len--;
					}
					d_len--;
				}
				break;
			default: printf("%i %i\n", rid[r_len-1], r_len); abort();
		}

		/* leave the loop if all non-terminals are replaced by terminals */
		if (r_len==1) break;
	}

	free(rid);
	free(depth);
	free(c1);
	free(c2);

	return T[0][len-1][0]->rel;
}


/*
 * add gaps in a string
 */
char *include_gap_columns(char *db, int *gap_col, int len, int gapnr, char *newdb)
{
	int i, j;

	strcpy(newdb, db);

	for( i=0; i<gapnr; i++ ) {
		for( j=len+i+1; j>=gap_col[i]; j-- )
			newdb[j] = newdb[j-1];
		newdb[gap_col[i]] = '-';
	}

	return newdb;
}


/*
 * print the reliabilities of the input alignment in a file
 * that can be printed as dotplot using 'drawplot' from Bjarne Knudsen
 * gap-columns are filled with 0-reliabilities
 */
int create_ppfile(double **paired, double *single, int *gap_col, int len, int gapnr, char *ppfile)
{
	FILE *fp;
	int i, j, k, l;

	if (ppfile != NULL)
	{
	    if ((fp = fopen(ppfile, "w")) == NULL)
	    {
	    	fprintf(stderr, "ppfile: Error in opening file '%s'\n", ppfile);
	    	return 1;
	    }

		fprintf(fp, "%i\n", len+gapnr);

		/* base pair reliabilities */
		k=0;
		for( i=0; i<len; i++ )
		{
			/* add gap column */
			while (1) {
				if (k < gapnr && gap_col[k] == i+k) {
					for( j=0; j<len+gapnr; j++ )
						fprintf(fp, "%10.9e ", .0);
					fprintf(fp, "\n");
					k++;
				}
				else
					break;
			}
			l=0;
			for( j=0; j<len; j++ )
			{
				/* add gap column */
				while (1) {
					if (l < gapnr && gap_col[l] == j+l) {
						fprintf(fp, "%10.9e ", .0);
						l++;
					}
					else
						break;
				}
				if (paired[i][j]==INFINITE_D)
					if (paired[j][i]==INFINITE_D)
						fprintf(fp, "%10.9e ", .0);
					else
						fprintf(fp, "%10.9e ", paired[j][i]);
				else
					fprintf(fp, "%10.9e ", paired[i][j]);
			}
			/* add gap column at the end of the line */
			while (l < gapnr) {
				fprintf(fp, "%10.9e ", .0);
				l++;
			}
			fprintf(fp, "\n");
		}
		/* add gap column at the end of the file */
		while (k < gapnr) {
			for( j=0; j<len+gapnr; j++ )
				fprintf(fp, "%10.9e ", .0);
			fprintf(fp, "\n");
			k++;
		}
		fprintf(fp, "\n");

		/* unpaired reliabilities */
		k=0;
		for( i=0; i<len; i++ ) {
			/* add gap column */
			while (1) {
				if (k < gapnr && gap_col[k] == i+k) {
					fprintf(fp, "%10.9e ", .0);
					k++;
				}
				else
					break;
			}
			fprintf(fp, "%10.9e ", single[i]);
		}
		/* add gap column at the end of the line */
		while (k < gapnr) {
			fprintf(fp, "%10.9e ", .0);
			k++;
		}
		fprintf(fp, "\n");

	    if (fclose(fp) != 0)
	    {
	    	fprintf(stderr, "ppfile: Error in closing file '%s'\n", ppfile);
	    	return 1;
	    }
	}

	return 0;
}


/*
 * ensemble diveristy
 */
double ensemble_diversity(double **paired, int len)
{
	int i, j, c=0;

	double ensdiv = .0;
	for( i=0; i<len; i++ )
	{
		for( j=i+4; j<len; j++ )
		{
			ensdiv += paired[i][j] * (1 - paired[i][j]);
			c++;
		}
	}

	//fprintf(stderr, "len = %i\tc = %i\n", len, c);
	return ensdiv/c;
}


/*
 * fill Nussinov-style data structure with all rule probabilities for all subregions
 */
double ***nussinov_subopt(double **paired_pet, double *single_pet, int len)
{
	double ***P;
	int i, j, k, t;

	/* initialize cube P which holds arrays
	 * column(i=1,..,len+1-j): base position;
	 * row(j=1,...,len): length of subsequence;
	 * depth(t): non terminals of SCFG: [0]=>'S', [1]=>'F', [2]=>'L'
	 * each field stores sum of all possible structure reliabilities in sub region
	 * -> kind of partition function
	*/
	P = (double***)malloc(len * sizeof(double**));
	for( i=0; i<len; i++ ) {
		P[i] = (double**)malloc(len * sizeof(double*));
		for( j=0; j<len; j++ ) {
			P[i][j] = (double*)malloc(3 * sizeof(double));
			/* Initialization of sum reliabilities */
			for( t=0; t<3; t++ )
				P[i][j][t] = INFINITE_D;
		}
	}

	/* write first row of cube T with reliabilities of unpaired bases in alignment */
	for( i=0; i<len; i++ )
	{
		/* L->s */
	    P[i][0][2] = 2*alpha * single_pet[i]/len;

	    /* S->L */
	    P[i][0][0] = P[i][0][2];
	}

	/* fill cube T */
	for( j=0; j<len; j++ )
	{
		for( i=0; i<len-j; i++ )
		{
			/* rules S->LS and F->LS (standard CYK algorithm) */
			for( k=1; k<j+1; k++ )
			{
				if (P[i][k-1][2]==INFINITE_D || P[i+k][j-k][0]==INFINITE_D)
					continue;

				/* S->LS: P[i,j][S] := P[i,j][S] + P[i,k-1][L] + P[i+k,j-k][S] */
				if( P[i][j][0]==INFINITE_D )
					P[i][j][0] = 0.;
				P[i][j][0] = add(P[i][j][0], add(P[i][k-1][2], P[i+k][j-k][0]));

				/* F->LS: P[i,j][F] := P[i,j][F] + P[i,k-1][L] + P[i+k,j-k][S] */
				if( P[i][j][1]==INFINITE_D )
					P[i][j][1] = 0.;
				P[i][j][1] = add(P[i][j][1], add(P[i][k-1][2], P[i+k][j-k][0]));
			}

			/*
			 * rule F->dFd (inner bond)
			 * loop consists of at least 3 unpaired bases
			 */
			if (j>3)
			{
				/* F->dFd: P[i,j][F] := P[i,j][F] + P[i+1,j-2][F] x Rel_paired[i,i+j] */
				if( P[i][j][1]==INFINITE_D )
					P[i][j][1] = 0;
				P[i][j][1] = add(P[i][j][1], add(P[i+1][j-2][1], 2 * paired_pet[i][i+j]/len));
			}

			/*
			 * rule L->dFd (open bond)
			 * loop consists of at least 3 unpaired bases
			 */
			if (j>3)
			{
				/* L->dFd: P[i,j][L] := P[i,j][L] + P[i+1][j-2][F] x Rel_paired[i,i+j] */
				if( P[i][j][2]==INFINITE_D )
					P[i][j][2] = 0;
				P[i][j][2] = add(P[i][j][2], add(P[i+1][j-2][1], 2 * paired_pet[i][i+j]/len));
			}

			/* S->L: P[i,j][S] := P[i,j][S] + P[i,j][L] */
			if (P[i][j][2]!=INFINITE_D)
			{
				if( P[i][j][0]==INFINITE_D )
					P[i][j][0] = 0;
				P[i][j][0] = add(P[i][j][0], P[i][j][2]);
			}
		}
	}

	if (P[0][len-1][0]==INFINITE_D)
	{
		fprintf(stderr, "The consensus RNA structure could not be predicted!\n");
		exit(1);
	}

	return P;
}


/*
 * backtracking through the full Nussinov-style data structure by structure sampling
 */
double backtracking_subopt(double ***P, double **paired_pet, double *single_pet, int len, char *struc)
{
	int rid[3*len];   /* stack of non-terminals: S->0, F->1, L->2 */
	int depth[len/5];   /* stack counting number of building rules for each non-terminals */
    int istack[3*len];   /* stack of i indices (start position) of subsequences */
    int jstack[3*len];   /* stack of j indices (length) of subsequences */
    unsigned long sam_prob[len];   /* stochastic probabilities of all rules that can occur */
    int klist[len];   /* list of indices k which separate sequence in two subsequences */
	int s_len, r_len, d_len, ij_len, sam_len;
	int i, j, k, offset, sample;
	double score;
	int CONST = 1000000;
	struct tagSCFG scfgProbs;

	s_len = r_len = d_len = ij_len = sam_len = 0;
	struc[s_len++] = 'S';
	struc[len] = '\0';
	rid[r_len++] = 0;
	depth[d_len++] = 0;
	istack[ij_len] = 0;
	jstack[ij_len++] = len-1;
	score = 0;

	/* stochastic context-free grammar (scfg) rule probabilities taken from Knudsen et al. (2003) */
	scfgProbs.S_LS = 0.868534;
	scfgProbs.S_L = 0.131466;
	scfgProbs.F_dFd = 0.787640;
	scfgProbs.F_LS = 0.212360;
	scfgProbs.L_s = 0.894603;
	scfgProbs.L_dFd = 0.105397;

	while(1)
	{
		i = istack[ij_len-1];
		j = jstack[ij_len-1];

		//printf("%s %i %i %i %i %i %i %i %i %8.2f\n", struc, s_len, r_len, d_len, ij_len, i, j, rid[r_len-1], depth[d_len-1], score);

		switch( rid[r_len-1] )
		{
			/* actual non-terminal is S */
			case 0:
				/* P[i,j][S] =: sample { P[i,j][L] if j=0                           : S->L->s,
				 *                       P[i+1][j-2][F] + Rel_paired[i,i+j]         : S->L->dFd,
				 *                       P[i,k-1][L] + P[i+k,j-k][S] forall 1<=k<=j : S->LS }
				 */
				if (j==0)
				{
					offset = strchr(struc, 'S') - struc;
					struc[offset] = '.';
					score = add(score, 2*alpha * single_pet[i]/len);

					for( k=0; k<depth[d_len-1]; k++ )
					{
						ij_len--;
						r_len--;
					}
					d_len--;
				}
				else
				{
					/* sample from rules S->LS and S->L->dFd */
					sam_len = 0;

					for( k=1; k<j+1; k++ )
					{
						if( P[i][k-1][2]!=INFINITE_D && P[i+k][j-k][0]!=INFINITE_D )
						{
							klist[sam_len] = k;
							sam_prob[sam_len++] = round( CONST * scfgProbs.S_LS * add( P[i][k-1][2]/pow(2,k), P[i+k][j-k][0]/pow(2,j-k+1) ) );
							//sam_prob[sam_len++] = round( CONST * add( P[i][k-1][2]/pow(2,k), P[i+k][j-k][0]/pow(2,j-k+1) ) );
							//printf("%8.6f %8.6f %i %li\n", P[i][k-1][2]/pow(2,k), P[i+k][j-k][0]/pow(2,j-k+1), k, sam_prob[sam_len-1]);
						}
					}
					change_list_sum_to_max(sam_prob, 2, sam_len);
					//for( k=1; k<sam_len+1; k++ ) printf("%i %li\n", k, sam_prob[k-1]);

					if( j>3 && P[i+1][j-2][1]!=INFINITE_D)
					{
						sam_prob[sam_len++] = round( CONST * scfgProbs.S_L * scfgProbs.L_dFd * add( P[i+1][j-2][1]/pow(2,j-1), 2 * paired_pet[i][i+j]/len ) );
						//sam_prob[sam_len++] = round( CONST * add( P[i+1][j-2][1]/pow(2,j-1), 2 * paired_pet[i][i+j]/len ) );
						//printf("%8.6f %8.6f bp %li\n", P[i+1][j-2][1]/pow(2,j-1), paired_pet[i][i+j], sam_prob[sam_len-1]);
					}
					else
					{
						sam_prob[sam_len++] = 0;
						//printf("bp0 %li\n", sam_prob[sam_len-1]);
					}

					sample = statistical_sampling(sam_prob, sam_len);
					//printf("sample: %i %i %i\n", sample, sam_len, klist[sample]);

					if( sample<sam_len-1 )
					{
						/* S->LS */
						s_len++;
						offset = strchr(struc, 'S') - struc;
						for( k=s_len; k>offset+1; k-- )
							struc[k] = struc[k-1];
						struc[offset] = 'L';
						struc[offset+1] = 'S';

						/* S rule */
						rid[r_len++] = 0;
						istack[ij_len] = i+klist[sample];
						jstack[ij_len++] = j-klist[sample];
						depth[d_len-1]++;

						/* L rule */
						rid[r_len++] = 2;
						istack[ij_len] = i;
						jstack[ij_len++] = klist[sample]-1;
						depth[d_len++] = 1;
					}
					else
					{
						//if (j<=3 || P[i+1][j-2][1]==INFINITE_D) break;

						/* S->L->dFd */
						s_len += 2;
						offset = strchr(struc, 'S') - struc;
						for( k=s_len; k>offset+2; k--)
							struc[k] = struc[k-2];
						struc[offset] = '(';
						struc[offset+1] = 'F';
						struc[offset+2] = ')';
						score = add(score, 2 * paired_pet[i][i+j]/len);

						rid[r_len++] = 1;
						istack[ij_len] = i+1;
						jstack[ij_len++] = j-2;
						depth[d_len-1]++;
					}
				}

				break;
			/* actual non-terminal is F */
			case 1:
				/* P[i,j][F] =: sample { P[i+1,j-2][F] + Rel_paired[i,i+j]          : F->dFd,
				 *                       P[i,k-1][L] + P[i+k,j-k][S] forall 1<=k<=j : F->LS }
				 */
				/* sample from rules F->LS and F->dFd */
				sam_len = 0;
				for( k=1; k<j+1; k++ )
				{
					if( P[i][k-1][2]!=INFINITE_D && P[i+k][j-k][0]!=INFINITE_D )
					{
						klist[sam_len] = k;
						sam_prob[sam_len++] = round( CONST * scfgProbs.F_LS * add( P[i][k-1][2]/pow(2,k), P[i+k][j-k][0]/pow(2,j-k+1) ) );
						//sam_prob[sam_len++] = round( CONST * add( P[i][k-1][2]/pow(2,k), P[i+k][j-k][0]/pow(2,j-k+1) ) );
						//printf("%i %li\n", k, sam_prob[sam_len-1]);
					}
				}
				change_list_sum_to_max(sam_prob, 2, sam_len);
				//for( k=1; k<sam_len+1; k++ ) printf("%i %li\n", k, sam_prob[k-1]);

				if( j>3 && P[i+1][j-2][1]!=INFINITE_D)
				{
					sam_prob[sam_len++] = round( CONST * scfgProbs.F_dFd * add( P[i+1][j-2][1]/pow(2,j-1), 2 * paired_pet[i][i+j]/len ) );
					//sam_prob[sam_len++] = round( CONST * add( P[i+1][j-2][1]/pow(2,j-1), 2 * paired_pet[i][i+j]/len ) );
					//printf("bp %li\n", sam_prob[sam_len-1]);

				}
				else
					sam_prob[sam_len++] = 0;

				sample = statistical_sampling(sam_prob, sam_len);
				//printf("sample: %i %i %i\n", sample, sam_len, klist[sample]);

				if( sample<sam_len-1 )
				{
					/* S->LS */
					s_len++;
					offset = strchr(struc, 'F') - struc;
					for( k=s_len; k>offset+1; k-- )
						struc[k] = struc[k-1];
					struc[offset] = 'L';
					struc[offset+1] = 'S';

					/* S rule */
					rid[r_len++] = 0;
					istack[ij_len] = i+klist[sample];
					jstack[ij_len++] = j-klist[sample];
					depth[d_len-1]++;

					/* L rule */
					rid[r_len++] = 2;
					istack[ij_len] = i;
					jstack[ij_len++] = klist[sample]-1;
					depth[d_len++] = 1;
				}
				else
				{
					//if (j<=3 || P[i+1][j-2][1]==INFINITE_D) break;

					/* S->L->dFd */
					s_len += 2;
					offset = strchr(struc, 'F') - struc;
					for( k=s_len; k>offset+2; k--)
						struc[k] = struc[k-2];
					struc[offset] = '(';
					struc[offset+1] = 'F';
					struc[offset+2] = ')';
					score = add(score, 2 * paired_pet[i][i+j]/len);

					rid[r_len++] = 1;
					istack[ij_len] = i+1;
					jstack[ij_len++] = j-2;
					depth[d_len-1]++;
				}

				break;
			/* actual non-terminal is L */
			case 2:
				/* P[i,j][L] =: sample { P[i,j][L] if j=0                           : L->s,
				 *                       P[i+1][j-2][F] + Rel_paired[i,i+j]         : L->dFd }
				 */
				if (j==0)
				{
					offset = strchr(struc, 'L') - struc;
					struc[offset] = '.';
					score = add(score, 2*alpha * single_pet[i]/len);

					for( k=0; k<depth[d_len-1]; k++ )
					{
						ij_len--;
						r_len--;
					}
					d_len--;
				}
				else
				{
					if (j>3)
					{
						s_len += 2;
						offset = strchr(struc, 'L') - struc;
						for( k=s_len; k>offset+2; k--)
							struc[k] = struc[k-2];
						struc[offset] = '(';
						struc[offset+1] = 'F';
						struc[offset+2] = ')';
						score = add(score, 2 * paired_pet[i][i+j]/len);

						rid[r_len++] = 1;
						istack[ij_len] = i+1;
						jstack[ij_len++] = j-2;
						depth[d_len-1]++;
					}
					else
					{
						fprintf(stderr, "Unpaired bases were not reached!\n");
						exit(1);
					}
				}

				break;
			default: printf("%i %i\n", rid[r_len-1], r_len); abort();
		}

		/* leave the loop if all non-terminals are replaced by terminals */
		if (r_len==1) break;
	}

	/*free(rid);
	free(depth);
	free(klist);*/

	return score;
}


/*
 * generate a random number from 0 to sum of entries in sam_prob
 */
int statistical_sampling(unsigned long *sam_prob, int sam_len)
{
	long sum = 0;
	int i, sample;

	for( i=0; i<sam_len; i++ )
		sum += sam_prob[i];

	sample = rand()%sum;
	//printf("Sample: %i\n",sample);

	sum = 0;
	for( i=0; i<sam_len; i++ )
	{
		sum += sam_prob[i];
		if (sample<=sum)
			return i;
	}

	return 0;
}


/*
 * find maximal value in the input list from index start to end
 * and the contribution (ratio) to the input sum of all values in the input list from index start to end
 * and adapts the list so that its output sum is the input maximal value by calculating each output value as its input ratio times the output sum
 */
void change_list_sum_to_max(unsigned long *sam_prob, int start, int end)
{
	int i;
	int sum = 0;
	int maxi = 0;

	for( i=start-1; i<end; i++ )
	{
		sum += sam_prob[i];
		if (max(sam_prob[i], maxi))
			maxi = sam_prob[i];
	}

	for( i=start-1; i<end; i++ )
		sam_prob[i] = maxi * sam_prob[i] / sum;
}


/*
 * Free Memory of structure Aln
 */
void FreeAln(Aln *align)
{
	int i, nr;

	nr = align->nr;

	for( i=0; i<nr; i++ )
		free(align->identifier[i]);
	free(align->identifier);

	for( i=0; i<nr; i++ )
		free(align->sequence[i]);
	free(align->sequence);

	free(align);
}


/*
 * Free Memory of structure SeqList
 */
void FreeSeqList(SeqList *seqlist)
{
	int i, nr;

	nr = seqlist->nr;

	free(seqlist->len);

	for( i=0; i<nr; i++ )
		free(seqlist->sequence[i]);
	free(seqlist->sequence);

	for( i=0; i<nr; i++ )
		free(seqlist->origid[i]);
	free(seqlist->origid);

	free(seqlist);
}

/*
 * Initialize structure tagPartStruc
 */
 void InitPartStruc(PartStruc *evocon)
 {
	 evocon->ss_nr = 0;
	 evocon->bp_nr = 0;
 }

/*
 * Free Memory of structure PartStruc
 */
void FreePartStruc(PartStruc *evocon)
{
	free(evocon->bp_left);
	free(evocon->bp_right);
	free(evocon->ss);

	free(evocon);
}


/*
 * Free Memory of structure DynEntry
 */
void FreeDynEntry(DynEntry ***T, int len)
{
	int i, j, t;

	for ( i=0; i<len; i++ ) {
		for ( j=0; j<len; j++ ) {
			for( t=0; t<3; t++ )  {
				free(T[i][j][t]);
			}
			free(T[i][j]);
		}
		free(T[i]);
	}
	free(T);
}


/*
 * FreeMemory of two dimensional array of length len
 */
void FreeArray(double **array, int len)
{
	int i;

	if(array) {
		for (i=0; i<len; i++)
			free(array[i]);
		free(array);
	}
}


/*
 * replace newline by string terminator
 */
void chomp(char *s)
{
	s[strcspn(s, "\n")] = '\0';
}


/*
 * help output if help_flag is set or no FASTA-file as input
 */
void usage()
{
	printf("PETfold v2.2\n"
	 "============\n"
	 "by Stefan E Seemann (seemann@rth.dk)\n"
	 "Reference: Seemann et al. Nucleic Acids Res. 36(20):6355-62, 2008\n"
	 "Web service: http://rth.dk/resources/petfold\n"
	 "\n"
	 "Usage:\n"
	 "  PETfold -f <file> [ options ] [ parameter settings ]\n"
	 "\n"
	 "  -f --fasta <file>          ... alignment in fasta format\n"
	 "Options:\n"
	 "  -s --setstruc <struc|file> ... calculates score for given structure in dot-bracket notation\n"
	 "  -t --settree <tree|file>   ... calculates score for given tree in Newick tree format\n"
	 "  --war                      ... fasta format output\n"
	 "  -r --ppfile <file>         ... writes PET reliabilities in file\n"
     "  -o --suboptimal <nr>       ... number of alternative structures found by sampling\n"
	 "  --ppfold                   ... reads ppfold calculated reliabilities from STDIN\n"
	 "  --verbose                  ... writes long output\n"
	 "  --help                     ... this output\n"
	 "Parameter settings:\n"
	 "  -p --setevocon_bp <reliab> ... reliab.threshold for conserved base pairs (default: 0.9)\n"
	 "  -u --setevocon_ss <reliab> ... reliab.threshold for conserved unpaired bases (default: 1)\n"
	 "  -a --setalpha <nr>         ... weighting factor for unpaired reliabilities (default: 0.2)\n"
	 "  -b --setbeta <nr>          ... weighting factor for thermodynamic overlap (default: 1)\n"
	 "  -g --setgap <nr>           ... max. percent of gaps in alignment column (default:0.25)\n"
	 "\n"
	 "You may set the environment variable PETFOLDBIN to path of files scfg.rate and article.grm:\n"
	 "export PETFOLDBIN=<PATH>\n"
	 "\n"
	 "Parse sampling of suboptimal structures:\n"
         "PETfold -f <file> -o 10000 | grep Subopt | sort -k 4 -ur | awk '{print $3,$4}' | less\n\n"
	);

	exit(1);
}
