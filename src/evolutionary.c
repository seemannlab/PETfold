/*
 * evolutionary.c
 *
 *  Created on: 06.02.2011
 *      Author: Stefan Seemann
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "evolutionary.h"
#include "petfoldlibs.h"
#include "file.h"
#include "phyl.h"
#include "grammar.h"
#include "col.h"
#include "align.h"
#include "newcolprob.h"
#include "inout.h"

int *get_prob_tree(Aln * align, double **paired_tree, double *single_tree, char *settree)
{
	PartStruc *dummy = malloc(1 * sizeof(PartStruc));
	InitPartStruc(dummy);
	int *struc = __get_pfold(align, paired_tree, single_tree, dummy, settree);
	free(dummy);

	return struc;
}

int *get_prob_tree_constraint(Aln* align, double **paired_tree, double *single_tree, PartStruc *partstruc, char *settree)
{
	return __get_pfold(align, paired_tree, single_tree, partstruc, settree);
}

int *__get_pfold(Aln *align_gf, double **paired_tree, double *single_tree, PartStruc *partstruc, char *settree)
{
	FILE *grammarfp;
	Grammar *grammar;
	Entry **entry_list;
	Entry *entry;
	Align *align;
	Dist *dist;
	Phyl *tmpphyl, *phyl;
	double *collikelhd;	/* Array of column probabilities to be created by the tree */
	double oldnll, nll;
	//char *file;
	char file[200];
	Alignstr *alignstr;
	int wlen;
	Colprob *cp;
	Inside *e;
	Outside *f;
	Pairprob *pp;
	PpCyk *ppcyk;
	double logp_limit;
	//Edouble totprob; /* column likelihood */
    char field[MAXCOLW];
    int align_bp_col, alignpos_col, seq_bp_col, seqpos_col;
    double robust;
	int i, j, p, w;
	int *struc = NULL;

	/* Pfold's default settings */
	logp_limit = -3;
	robust = 0.01;
	//totprob = Dbl2Edbl(0);

	if( setevo_flag )
	{
		/* open rate file
		 * read grammar */
		if( getenv("PETFOLDBIN") )
			sprintf(file, "%s/scfg.rate", getenv("PETFOLDBIN"));
		else
			sprintf(file, "scfg.rate");
		if ( (grammarfp = fopen(file, "r")) == NULL) {
			fprintf(stderr, "findphyl: Error in opening file '%s'\n", file);
		    exit(1);
		}
	    grammar = ReadGrammar(grammarfp);
	    if (fclose(grammarfp) != 0) {
	        fprintf(stderr, "scfg: Error in closing grammar file\n");
	        exit(1);
	    }

	    /* copy Aln object in Entry object */
	    //entry_list = (Entry **)malloc((align_gf->nr+2)*sizeof(Entry *));
			entry_list = (Entry **)malloc((align_gf->nr+1)*sizeof(Entry *));
	    Aln2EntryList(align_gf, entry_list);
		  //for (i = 0; entry_list[i] != NULL; i++) PrintEntry(stdout, entry_list[i]);

	    /* read the sequences */
	    align = Col2Align(entry_list, (int (*)(char, void *))FindSym,  (void *)grammar);

	    /* phylogenetical tree is not given and is estimated by neighbour joining approach
		 * adapted from $Pfold:findphyl.c
		 */
		if( ! strcmp("", settree) )
		{
		    /* find phylogeny of sequences
		     * using the neighbour joining approach */
		    dist = AlignDist(align, grammar);
		    tmpphyl = Neighbour(dist);
		    FixPhyl(tmpphyl);
		    entry = PhylEntry(tmpphyl, "tree");
				//PrintEntry(stdout, entry);
		    phyl = ReadColPhyl(entry);

				/* free memory */
				FreeDist(dist);
		}
		/* END adapted from Pfold:findphyl.c */
		else
		{
			char *tmptree = (char *) malloc(BUFFERSIZE * sizeof(char));
			//char* tmptree = (char *)malloc((strlen(settree)+3) * sizeof(char));
			strcpy(tmptree,settree);
			strcat(tmptree, ":0");
			tmpphyl = ReadPhyl(tmptree, 1);
			free(tmptree);
		  entry = PhylEntry(tmpphyl, "tree");
		  phyl = ReadColPhyl(entry);
		}
		/* free memory */
		FreePhyl(tmpphyl);

		if (verbose_flag)
		{
			printf("PFOLD:\n NJ tree = ");
			PrintPhyl(stdout, phyl);
		}

		/* branch lengths are not given and are estimated by maximum likelihood
		 * adapted from Pfold:mltree.c
		 */;
		if( settree == 0 || ! strcmp("", settree) || strstr(settree, ":") == NULL )
		{
			initphyl(phyl->root, grammar, align);
			collikelhd = (double *)calloc(align->len, sizeof(double));
			oldnll = -1;
			nll = -1;
			while (oldnll == -1 || oldnll-nll > 0.1) {
				oldnll = nll;
				setmat(phyl->root, grammar); /* Setting matrices */
				downcalc(phyl->root, align); /* Down calculations */
				upcalc(phyl->root, align); /* Up calculations */
				if (verbose_flag)
					fprintf(stderr, " Optimize tree:");
				nll = optimize(phyl->root, align, grammar, phyl, collikelhd); /* Optimize */
				if (verbose_flag)
					fprintf(stderr, " NLL = %f\n", nll);
			}
			free(collikelhd);
			FixPhyl(phyl);
			UpdatePhylEntry(phyl, entry);
			phyl = ReadColPhyl(entry);

			/* store tree in Newick format if not given */
			//Phyl2Char(&settree, phyl);
			Phyl2Char(settree, phyl);
		}
		/* END adapted from Pfold:mltree.c */

		if (verbose_flag)
		{
			printf(" ML tree = ");
			PrintPhyl(stdout, phyl);
			//printf(" ML tree = %s\n", settree);
		}

		/* calculation of structure reliabilities using the tree, the evolutionary model and the grammar (SCFG)
		 * adapted from Pfold:scfg.c
		 */

		/* open grammar file
		 * read grammar */
		if( getenv("PETFOLDBIN") )
			sprintf(file, "%s/article.grm", getenv("PETFOLDBIN"));
		else
			sprintf(file, "article.grm");
		if ((grammarfp = fopen(file, "r")) == NULL) {
		    fprintf(stderr, "scfg: Error in opening file '%s'\n", file);
		    exit(1);
		}
		grammar = ReadGrammar(grammarfp);
		grammar->mindist = 4;
	    if (fclose(grammarfp) != 0) {
	        fprintf(stderr, "scfg: Error in closing grammar file\n");
	        exit(1);
	    }
	    RobustGrammar(grammar, robust);

	    alignstr = (Alignstr *)malloc(sizeof(Alignstr));
	    alignstr->align = align;

	    /* initialize Pfold predicted RNA secondary structure */
	    alignstr->str = (int *)malloc(align_gf->len * sizeof(int));

	    if( partstruc->ss_nr == 0 && partstruc->bp_nr == 0 ) {
				align_bp_col = ReadColno(entry_list[0], "align_bp");
				alignpos_col = ReadColno(entry_list[0], "alignpos");
				seq_bp_col = ReadColno(entry_list[0], "seq_bp");
				seqpos_col = ReadColno(entry_list[0], "seqpos");
				if (align_bp_col == 0 && seq_bp_col == 0)
			  	for (i = 0; i < align_gf->len; i++)
				  	alignstr->str[i] = -1;
				else
			  	for (i = 0; i < align_gf->len; i++) {
				  	if (align_bp_col != 0)
					  	GetField(field, entry_list[0], i+1, align_bp_col);
				  	else
					  	GetField(field, entry_list[0], i+1, seq_bp_col);
				  	if (StrCmp(field, "s") == 0)
					  	alignstr->str[i] = SINGLE;
				  	else
					  	if (StrCmp(field, "d") == 0)
						 		alignstr->str[i] = DOUBLE;
					  	else {
						 		alignstr->str[i] = FindPair(entry_list[0], i+1,
						  	align_bp_col, alignpos_col,
						  	seq_bp_col, seqpos_col)-1;
					  	}
				  }
	    }
	    else {
	    	int *struc = partstruc2str(partstruc, align_gf->len);
	    	for (i=0; i < align_gf->len; i++) {
	    		alignstr->str[i] = -1;
	    		if (struc[i] != INFINITE_I) {
	    			if (struc[i] == -1)
	    				alignstr->str[i] = SINGLE;
	    			else
	    				if (struc[i] >= i+grammar->mindist || struc[i]+grammar->mindist <= i)
	    					alignstr->str[i] = struc[i];
	    		}
	    	}
	    	free(struc);
	    }

		alignstr->map = (int *)malloc((align_gf->len+1) * sizeof(int));
	    for (i = 0; i < align_gf->len; i++)
	    	alignstr->map[i] = i+1;
	    alignstr->map[i] = align_gf->len;
	    wlen = alignstr->align->len;

	    /* allocate memory */
	    e = MakeInside(grammar, wlen);
	    f = MakeOutside(grammar, wlen);
	    cp = MakeColprob(grammar, wlen);
	    pp = MakePairprob(wlen);
	    ppcyk = MakePpCyk(wlen);

	    /* inside-outside algorithm to find most occurring structure (MEA) */
	    InitCol(phyl, grammar, wlen, alignstr, logp_limit);
	    InitColprob(cp, alignstr);
	    InitInside(e, cp);
	    InitOutside(e, f, cp);
	    //totprob = InitPairprob(pp, e, f, cp);
	    InitPairprob(pp, e, f, cp);
	    InitPpCyk(ppcyk, pp);
	    alignstr->str = PpCykStr(ppcyk);
	    //for(i=0;i<align_gf->len;i++) printf("%i %i\t",alignstr->map[i],alignstr->str[i]);

	    /* add structure to list of Entries */
	    //for (i = 0; entry_list[i] != NULL; i++)
	    	//AddStr(entry_list[i], alignstr, pp);
	    //for (i = 0; entry_list[i] != NULL; i++) PrintEntry(stdout, entry_list[i]);

	    /* write reliability matrix */
	    //PrintPairprob(pairprobfp, pp, pp, alignstr, totprob);
	    for (i=0; i<align_gf->len; i++)
		{

			for (j=0; j<align_gf->len; j++)
			{
				w = (i > j) ? i - j : j - i;
				p = (i > j) ? j : i;
				paired_tree[i][j] = Edbl2Dbl(pp->dbl[p][w]);
			}

			single_tree[i] = Edbl2Dbl(pp->sgl[i]);
		}

		/* free memory */
		for( i=0; i<align_gf->nr; i++ )
			FreeEntry(entry_list[i]);
		free(entry_list[i]);
		free(entry_list);

	  struc = (int *)malloc(align_gf->len * sizeof(int));
		memcpy(struc, alignstr->str, align_gf->len * sizeof(int));

		free(alignstr->map);
		free(alignstr->str);
		free(alignstr);

		FreeAlign(align);
	}
	/* evolutionary model is switched off */
	else
	{
		/* write reliabilities to interface */
		for (i=0; i<align_gf->len; i++)
		{
			for (j=0; j<align_gf->len; j++)
				paired_tree[i][j] = .0;

			single_tree[i] = .0;
		}
		return 0;
	}

	/* free dynamic memory */
	FreePhyl(phyl);
	FreePpCyk(ppcyk);
	FreePairprob(pp);

	return struc;
}


int *get_ppfold(int len, double **paired_tree, double *single_tree)
{
	int nr1 = 16*len*(len-1)/2+16*len+1;
	int nr2 = 8*len+1;
	char *string;
	string = (char *)malloc(nr1 * sizeof(char));
	char a[17];
	char *string2;
	string2 = (char *)malloc(nr2 * sizeof(char));
	char b[9];
	int *str;

	/* read STDIN: time java -jar ../ppfold/PPfold_forPETfold.jar myRRE.fasta --onlyCT -t treefilename.newick | ./PETfold -f myRRE.fasta --ppfold */
	if( fgets(string, nr1, stdin) == NULL ) {
		fprintf(stderr, "PPfold input stream is not readable!\n");
		exit(1);
	}
	if( fgets(string2, nr2, stdin) == NULL ) {
		fprintf(stderr, "PPfold input stream is not readable!\n");
		exit(1);
	}

	int i,j;
	int count = 0;

	/* read ppfold reliabilities */
	for(i=0; i<len; i++) {
		for(j=1; j<len-i; j++) {
			sprintf(a, "%.16s", &string[count]);
			paired_tree[i][j+i] = (double) atof(a);
			paired_tree[j+i][i] = (double) atof(a);
			//fprintf(stderr,"%16.8e",paired_tree[i][j+i]);
			count += 16;
		}
	}

	for(i=0; i<len; i++) {
		sprintf(a, "%.16s", &string[count]);
		single_tree[i] = (double) atof(a);
		count +=16;
	}

	/* read ppfold structure */
	str = (int *)malloc(len * sizeof(int));
	count = 0;
	for(i=0; i<len; i++) {
		sprintf(b, "%.8s", &string2[count]);
		//fprintf(stderr, "%s %i\n", b, atoi(b));
		str[i] = (int) atoi(b);
		count += 8;
	}

	free(string);
	free(string2);

	return str;
}


/*
 * Converts Aln object in a list of Entry objects.
 */
int Aln2EntryList(Aln *align, Entry **entry_list)
{
  int i, j, k;
  Entry *entry;
  char *type = "RNA";

  for (i=0; i<align->nr; i++) {
	  if ((entry = MyNewEntry(type, align->identifier[i], align->len)) == NULL)
		  return 1;

	  k = 0;
	  for (j=0; j<align->len; j++)
	  {
		  if (align->sequence[i][j] != '-')
		  {
			  k++;
			  sprintf(entry->line[j], "N %c %i %i\n", align->sequence[i][j], k, j+1);
		  }
		  else {
			  sprintf(entry->line[j], "N %c . %i\n", align->sequence[i][j], j+1);
		  }
	  }

	  entry_list[i] = entry;
  }
  entry_list[i] = NULL;

  return 0;
}


/*
 * Makes an entry, ready for filling out.
 * Adapted from Pfold:col.c
 */
Entry *MyNewEntry(char *type, char *name, int len)
{
  Entry *entry;
  int i;

  entry = MakeEntry();

  entry->textlen = 7;
  entry->text = (char **)malloc(entry->textlen * sizeof(char *));

  entry->text[0] = (char *)malloc((22+strlen(type)) * sizeof(char *));
  sprintf(entry->text[0], "; TYPE              %s\n", type);

  entry->text[1] = (char *)malloc((22+strlen(name)) * sizeof(char *));
  sprintf(entry->text[1], "; COL 1             label\n");

  entry->text[2] = (char *)malloc((22+strlen(type)) * sizeof(char *));
  sprintf(entry->text[2], "; COL 2             residue\n");

  entry->text[3] = (char *)malloc((22+strlen(name)) * sizeof(char *));
  sprintf(entry->text[3], "; COL 3             seqpos\n");

  entry->text[4] = (char *)malloc((22+strlen(type)) * sizeof(char *));
  sprintf(entry->text[4], "; COL 4             alignpos\n");

  entry->text[5] = (char *)malloc((22+strlen(name)) * sizeof(char *));
  sprintf(entry->text[5], "; ENTRY             %s\n", name);

  entry->text[6] = (char *)malloc(14 * sizeof(char *));
  sprintf(entry->text[6], "; ----------\n");

  entry->len = len;

  entry->line = (char **)malloc(len * sizeof(char *));
  for (i = 0; i < len; i++) {
    entry->line[i] = (char *)malloc(2 * sizeof(char *));
    strcpy(entry->line[i], "\n");
  }

  entry->endline = (char *)malloc(14 * sizeof(char *));
  strcpy(entry->endline, "; **********\n");

  return entry;
}


/*
 * Sets alignment names, and numbers.
 * Copied from Pfold:mltree.c
 */
void initphyl(PhylNode *pnode, Grammar *grammar, Align *align)
{
  PhylNode *child;
  NodeInfo *pn_info;
  char *name;
  int numterm;
  int len;
  int i;

  numterm = grammar->term->count;

  name = (char *)pnode->elm;

  pnode->elm = (NodeInfo *)malloc(sizeof(NodeInfo));

  if (pnode->uplen < 0.001)
    pnode->uplen = 0.001;

  pn_info = (NodeInfo *)pnode->elm;
  len = align->len;

  pn_info->seqnum = SeqNumber(align, name);
  pn_info->grammar = grammar;
  pn_info->pmatrix = NULL;

  pn_info->p_bot_down = (Matrix **)malloc(len * sizeof(Matrix *));
  for (i = 0; i < len; i ++)
    pn_info->p_bot_down[i] = MakeMatrix(1, numterm);

  pn_info->p_top_down = (Matrix **)malloc(len * sizeof(Matrix *));
  for (i = 0; i < len; i ++)
    pn_info->p_top_down[i] = NULL;

  pn_info->p_bot_up = (Matrix **)malloc(len * sizeof(Matrix *));
  for (i = 0; i < len; i ++)
    pn_info->p_bot_up[i] = NULL;

  pn_info->p_top_up = (Matrix **)malloc(len * sizeof(Matrix *));
  for (i = 0; i < len; i ++)
    pn_info->p_top_up[i] = MakeMatrix(1, numterm);

  for (child = pnode->child; child != NULL; child = child->brother)
    initphyl(child, grammar, align);
}


/*
 * Calculates evolutionary matrices.
 * Copied from Pfold:mltree.c
 */
void setmat(PhylNode *pnode, Grammar *grammar)
{
  Sgrp *sgrp;
  LListCounter *lcount;
  PhylNode *child;
  NodeInfo *pn_info;
  Matrix *temp;

  //int number = pnode->number;
  //int seqnum = pn_info->seqnum;
  //int numterm = pn_info->grammar->term->count;
  pn_info = ((NodeInfo *)pnode->elm);

  lcount = MakeCounter(grammar->sgrp, FIRST);
  sgrp = Next(lcount);
  free(lcount);

  FreeMatrix(pn_info->pmatrix);
  pn_info->pmatrix =
    TransposeMatrix(temp = ExpMatrix(pnode->uplen, sgrp->eigen,
				     sgrp->diag, sgrp->inveigen));
  FreeMatrix(temp);

  for (child = pnode->child; child != NULL; child = child->brother)
    setmat(child, grammar);
}


/* Calculate probability distributions downwards in the tree
 * Copied from Pfold:mltree.c
 */
void downcalc(PhylNode *pnode, Align *align)
{
  PhylNode *child;
  NodeInfo *pn_info;
  int pos;
  int i;
  int numterm;
  int symno;

  pn_info = (NodeInfo *)(pnode->elm);

  numterm = pn_info->grammar->term->count;

  /* Do calculations for children */
  for (child = Child(pnode); child != NULL; child = Brother(child))
    downcalc(child, align);

  for (pos = 0; pos < align->len; pos++) {
    /* Set nucleotide probability vector */
    if (pn_info->seqnum == -1) { /* Internal node */
      for (i = 0; i < numterm; i++)
	pn_info->p_bot_down[pos]->entry[0][i] = Dbl2Edbl(1.);
    }
    else { /* Sequence here */
      symno = align->seq[pn_info->seqnum][pos];
      for (i = 0; i < numterm; i++)
	pn_info->p_bot_down[pos]->entry[0][i]
	  = pn_info->grammar->quickdist[symno][i];
    }

    /* Calculate evolution */
    for (child = Child(pnode); child != NULL; child = Brother(child))
      for (i = 0; i < numterm; i++)
	MulEdouble(&pn_info->p_bot_down[pos]->entry[0][i],
		   ((NodeInfo *)(child->elm))->p_top_down[pos]->entry[0][i]);

    FreeMatrix(pn_info->p_top_down[pos]);

    pn_info->p_top_down[pos]
      = MulMatrix(pn_info->p_bot_down[pos], pn_info->pmatrix);
  }
}


/* Calculate probability distributions upwards in the tree
 * Copied from Pfold:mltree.c
 */
void upcalc(PhylNode *pnode, Align *align)
{
  PhylNode *child;
  NodeInfo *pn_info;
  int pos;
  int i;
  int numterm;

  pn_info = (NodeInfo *)(pnode->elm);

  numterm = pn_info->grammar->term->count;

  for (pos = 0; pos < align->len; pos++) {

    /* Set nucleotide probability vector */
    for (i = 0; i < numterm; i++)
      pn_info->p_top_up[pos]->entry[0][i] = Dbl2Edbl(1.);

    /* Calculate evolution */
    if (pnode->parent != NULL) {
      for (i = 0; i < numterm; i++)
	MulEdouble(&pn_info->p_top_up[pos]->entry[0][i],
	     ((NodeInfo *)pnode->parent->elm)->p_bot_up[pos]->entry[0][i]);

      for (child = Child(pnode->parent); child != NULL; child = Brother(child))
	if (child != pnode) {
	  for (i = 0; i < numterm; i++)
	    MulEdouble(&pn_info->p_top_up[pos]->entry[0][i],
		    ((NodeInfo *)(child->elm))->p_top_down[pos]->entry[0][i]);
	}
    }

    FreeMatrix(pn_info->p_bot_up[pos]);

    pn_info->p_bot_up[pos]
      = MulMatrix(pn_info->p_top_up[pos], pn_info->pmatrix);
  }

  /* Do calculations for children */
  for (child = Child(pnode); child != NULL; child = Brother(child))
	  upcalc(child, align);
}


/* Optimize branch lengths
 * Copied from Pfold:mltree.c
 */
double optimize(PhylNode *pnode, Align *align, Grammar *grammar, Phyl *phyl, double *collikelhd)
{
  PhylNode *child;
  NodeInfo *pn_info;
  LListCounter *lcount;
  Sgrp *sgrp;
  Edouble e, sum;
  Matrix *temp, *p_temp, *old, *mat, *new;
  int pos;
  int i;
  int numterm;
  double nll;
  double l, x;
  int error;
  int finish;

  /* Do calculations for children */
  for (child = Child(pnode); child != NULL; child = Brother(child))
    nll = optimize(child, align, grammar, phyl, collikelhd);

  if (pnode->parent == NULL)
    return nll;

  lcount = MakeCounter(grammar->sgrp, FIRST);
  sgrp = Next(lcount);
  free(lcount);

  pn_info = (NodeInfo *)(pnode->elm);

  numterm = pn_info->grammar->term->count;

  p_temp = NULL;

  x = log(pnode->uplen+0.001);
  l = exp(x);
  finish = 0;
  InitMinimize(x-0.2, x+0.2, 0.001);
  while ((error = Minimize(&x, &nll)) == 0) {

    if (x < -10) {
      x = -10;
      finish = 10;
    }
    else if (x > 6) {
      x = 6;
      finish++;
    }
    if (finish > 4)
      break;
    l = exp(x);
    nll = 0;

    p_temp =
      TransposeMatrix(temp = ExpMatrix(l, sgrp->eigen,
				       sgrp->diag, sgrp->inveigen));
    FreeMatrix(temp);

    for (pos = 0; pos < align->len; pos++) {
      temp = MulMatrix(pn_info->p_bot_down[pos], p_temp);
      sum = Dbl2Edbl(0.);
      for (i = 0; i < numterm; i++) {
	e = ProdEdouble(sgrp->freq->entry[0][i],
			ProdEdouble(temp->entry[0][i],
				    pn_info->p_top_up[pos]->entry[0][i]));
	AddEdouble(&sum, e);
      }
      FreeMatrix(temp);

      nll -= Edbl2Dbl(LogEdouble(sum));

      //fprintf(stderr,"%f %20.18f\n",Edbl2Dbl(LogEdouble(sum)),Edbl2Dbl(sum));
      collikelhd[pos] = Edbl2Dbl(sum);
    }
    FreeMatrix(p_temp);
  }

  if (error == 1) {
    fprintf(stderr, "problems\n");
    exit(1);
  }

  p_temp = TransposeMatrix(temp = ExpMatrix(l, sgrp->eigen,
					    sgrp->diag, sgrp->inveigen));
  FreeMatrix(temp);

  pnode->uplen = exp(x);
  /* Adjust probabilities for relevant neighbours */

  if (pnode->parent != NULL) {
    for (pos = 0; pos < align->len; pos++) {
      old = pn_info->p_top_down[pos];
      new = MulMatrix(pn_info->p_bot_down[pos], p_temp);
      mat = ((NodeInfo *)pnode->parent->elm)->p_bot_down[pos];
      for (i = 0; i < numterm; i++) {
	if (Edbl2Dbl(old->entry[0][i]) != 0)
	  DivEdouble(&mat->entry[0][i], old->entry[0][i]);
	MulEdouble(&mat->entry[0][i], new->entry[0][i]);
      }

      for (child = Child(pnode->parent); child != NULL; child = Brother(child))
	if (child != pnode) {
	  mat = ((NodeInfo *)child->elm)->p_top_up[pos];
	  for (i = 0; i < numterm; i++) {
	    if (Edbl2Dbl(old->entry[0][i]) != 0)
	      DivEdouble(&mat->entry[0][i], old->entry[0][i]);
	    MulEdouble(&mat->entry[0][i], new->entry[0][i]);
	  }
	}
      FreeMatrix(new);
    }
  }

  FreeMatrix(p_temp);

  return nll;
}


/*
 * calculate probability of constrained (partial) structure in the evolutionary model
 */
double get_fold_evol_partial_prob(double evolconstraint, double evolensemble)
{
	return exp(evolconstraint-evolensemble);
}


/*
 * calculate probability of constrained ensemble (partial structure)
 */
double get_constraint_evolprob(Aln * align, char *tree, PartStruc *evocon)
{
	return __get_evolprob(align, tree, evocon);
}


/*
 * calculate probability of ensemble
 */
double get_ensemble_evolprob(Aln* align, char *tree)
{
	PartStruc *dummy = calloc(1, sizeof(PartStruc));

	double prob = __get_evolprob(align, tree, dummy);
	free(dummy);

	return prob;
}

/*
 * function used by get_constraint_evolprob and get_ensemble_evolprob
 */
double __get_evolprob(Aln* align, char *tree, PartStruc *evocon)
{
	FILE *grammarfp;
	Grammar *grammar;
	Entry **entry_list;
	Entry *entry;
	Align *alignment;
	Phyl *phyl;
	char file[200];
	Alignstr *alignstr;
	int wlen;
	Colprob *cp;
	Inside *e;
	Outside *f;
	Pairprob *pp;
	double logp_limit;
	Edouble totprob; /* column likelihood */
    char field[MAXCOLW];
    int align_bp_col, alignpos_col, seq_bp_col, seqpos_col;
    double robust;
	int i;

	/* Pfold's default settings */
	logp_limit = -3;
	robust = 0.01;
	totprob = Dbl2Edbl(0);

	/* open rate file
	 * read grammar */
	if( getenv("PETFOLDBIN") )
		sprintf(file, "%s/scfg.rate", getenv("PETFOLDBIN"));
	else
		sprintf(file, "scfg.rate");
	if ( (grammarfp = fopen(file, "r")) == NULL) {
		fprintf(stderr, "findphyl: Error in opening file '%s'\n", file);
	    exit(1);
	}
	grammar = ReadGrammar(grammarfp);
    if (fclose(grammarfp) != 0) {
        fprintf(stderr, "scfg: Error in closing grammar file\n");
        exit(1);
    }

    /* copy Aln object in Entry object */
    entry_list = (Entry **)malloc((align->nr+2)*sizeof(Entry *));
    Aln2EntryList(align, entry_list);
	//for (i = 0; entry_list[i] != NULL; i++) PrintEntry(stdout, entry_list[i]);

	/* read the sequences */
    alignment = Col2Align(entry_list, (int (*)(char, void *))FindSym,  (void *)grammar);

  /* initialize phylogeny by given tree */
	char* tmptree = (char *)malloc((strlen(tree)+3) * sizeof(char *));
	strcpy(tmptree, tree);
	strcat(tmptree, ":0");
	phyl = ReadPhyl(tmptree, 1);
	free(tmptree);
    entry = PhylEntry(phyl, "tree");
    phyl = ReadColPhyl(entry);
	/* END adapted from Pfold:mltree.c */

	/* calculation of structure reliabilities using the tree, the evolutionary model and the grammar (SCFG)
	 * adapted from Pfold:scfg.c
	 */

	/* open grammar file
	 * read grammar */
	if( getenv("PETFOLDBIN") )
		sprintf(file, "%s/article.grm", getenv("PETFOLDBIN"));
	else
		sprintf(file, "article.grm");
	if ((grammarfp = fopen(file, "r")) == NULL) {
	    fprintf(stderr, "scfg: Error in opening file '%s'\n", file);
	    exit(1);
	}
	grammar = ReadGrammar(grammarfp);
	grammar->mindist = 4;
    if (fclose(grammarfp) != 0) {
    	fprintf(stderr, "scfg: Error in closing grammar file\n");
	    exit(1);
	}
	RobustGrammar(grammar, robust);

	alignstr = (Alignstr *)malloc(sizeof(Alignstr));
	alignstr->align = alignment;

	/* initialize Pfold predicted RNA secondary structure */
	alignstr->str = (int *)malloc(align->len * sizeof(int));

	if( evocon->ss_nr == 0 && evocon->bp_nr == 0 ) {
		align_bp_col = ReadColno(entry_list[0], "align_bp");
		alignpos_col = ReadColno(entry_list[0], "alignpos");
		seq_bp_col = ReadColno(entry_list[0], "seq_bp");
		seqpos_col = ReadColno(entry_list[0], "seqpos");
		if (align_bp_col == 0 && seq_bp_col == 0)
			for (i = 0; i < align->len; i++)
				alignstr->str[i] = -1;
			else
				for (i = 0; i < align->len; i++) {
					if (align_bp_col != 0)
						GetField(field, entry_list[0], i+1, align_bp_col);
					else
						GetField(field, entry_list[0], i+1, seq_bp_col);
					if (StrCmp(field, "s") == 0)
						alignstr->str[i] = SINGLE;
					else
						if (StrCmp(field, "d") == 0)
							alignstr->str[i] = DOUBLE;
						else {
							alignstr->str[i] = FindPair(entry_list[0], i+1, align_bp_col, alignpos_col, seq_bp_col, seqpos_col)-1;
						}
				}
	}
	else {
		int *struc = partstruc2str(evocon, align->len);
		for (i=0; i < align->len; i++) {
			alignstr->str[i] = -1;
			if (struc[i] != INFINITE_I)
				if (struc[i] >= i+grammar->mindist || struc[i]+grammar->mindist <= i)
					alignstr->str[i] = struc[i];
		}
		free(struc);
	}

	alignstr->map = (int *)malloc((align->len+1) * sizeof(int));
	for (i = 0; i < align->len; i++)
		alignstr->map[i] = i+1;
    alignstr->map[i] = align->len;
    wlen = alignstr->align->len;

	/* allocate memory */
	e = MakeInside(grammar, wlen);
	f = MakeOutside(grammar, wlen);
	cp = MakeColprob(grammar, wlen);
	pp = MakePairprob(wlen);

	/* inside-outside algorithm to find most occurring structure (MEA) */
	InitCol(phyl, grammar, wlen, alignstr, logp_limit);
	InitColprob(cp, alignstr);
	InitInside(e, cp);
	InitOutside(e, f, cp);
	totprob = InitPairprob(pp, e, f, cp);

	/* free memory */
	for( i=0; i<align->nr+1; i++ )
		free(entry_list[i]);
	free(entry_list);

	free(alignstr->map);
	free(alignstr->str);
	free(alignstr);

	return Edbl2Dbl(LogEdouble(totprob));
}

/* Write phylogeny in Newick format to char */
/*void Phyl2Char(char** tree, Phyl *phyl)
{
	int offset = 0;
	char *ltree = *tree;

	offset = __Phyl2Char(ltree, offset, phyl->root);
	offset += snprintf(ltree+offset, BUFFERSIZE-offset-1, ";");
	ltree = (char *)realloc(ltree, offset+2 * sizeof(char *));
	*tree = ltree;
}*/
void Phyl2Char(char* tree, Phyl *phyl)
{
	int offset = 0;

	offset = __Phyl2Char(tree, offset, phyl->root);
	offset += snprintf(tree+offset, BUFFERSIZE-offset-1, ";");
}

/* Recursive function used by Phyl2Char */
int __Phyl2Char(char *tree, int offset, PhylNode *pnode)
{
	PhylNode *pnode2;

	if (pnode->child == NULL) { /* We have a leaf */
		offset += snprintf(tree+offset, BUFFERSIZE-offset-1, "%s", (char *)pnode->elm);
	}
	else {
		offset = offset + snprintf(tree+offset, BUFFERSIZE-offset-1, "(");
		if (pnode->elm != NULL) { /* We have an element */
			offset += snprintf(tree+offset, BUFFERSIZE-offset-1, "%s,", (char *)pnode->elm);
		}
		for (pnode2 = pnode->child; pnode2 != NULL; pnode2 = pnode2->brother) {
			offset = __Phyl2Char(tree, offset, pnode2);
			if (pnode2->brother != NULL)
				offset += snprintf(tree+offset, BUFFERSIZE-offset-1, ",");
		}
		offset += snprintf(tree+offset, BUFFERSIZE-offset-1, ")");
	}
	if (pnode->uplen != 0)
		offset += snprintf(tree+offset, BUFFERSIZE-offset-1, ":%g", pnode->uplen);

	return offset;
}
