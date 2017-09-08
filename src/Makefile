BINPATH = ../bin
HPATH = /opt/ViennaRNA/2.0.7/include/ViennaRNA/
LPATH = /opt/ViennaRNA/2.0.7/lib/libRNA.a

# Targets
all:       $(BINPATH)/PETfold $(BINPATH)/PETcofold $(BINPATH)/drawdot $(BINPATH)/drawphyl
petfold:   $(BINPATH)/PETfold
petcofold: $(BINPATH)/PETcofold
pfold: 	   $(BINPATH)/drawdot $(BINPATH)/drawphyl

$(BINPATH)/PETfold:	petfold.c file.o llist.o col.o grammar.o edouble.o matrix.o align.o phyl.o optimize.o inout.o newcolprob.o search.o
			#gcc -Wall -g -O1 -o $(BINPATH)/PETfold petfold.c petfoldlibs.o thermodynamic.o evolutionary.o file.o llist.o col.o grammar.o edouble.o matrix.o align.o phyl.o optimize.o inout.o newcolprob.o search.o $(LPATH) -lm
			gcc -Wall -o $(BINPATH)/PETfold petfold.c petfoldlibs.c thermodynamic.c evolutionary.c file.o llist.o col.o grammar.o edouble.o matrix.o align.o phyl.o optimize.o inout.o newcolprob.o search.o -I$(HPATH) $(LPATH) -lm -fopenmp

$(BINPATH)/PETcofold: 	petcofold.c petfoldlibs.o thermodynamic.o evolutionary.o file.o llist.o col.o grammar.o edouble.o matrix.o align.o phyl.o optimize.o inout.o newcolprob.o search.o
			gcc -Wall -fopenmp -o $(BINPATH)/PETcofold petcofold.c petcofoldlibs.c petfoldlibs.o thermodynamic.o evolutionary.o file.o llist.o col.o grammar.o edouble.o matrix.o align.o phyl.o optimize.o inout.o newcolprob.o search.o -I$(HPATH) $(LPATH) -lm

$(BINPATH)/drawphyl:	drawphyl.c file.o llist.o phyl.o grammar.o optimize.o col.o align.o matrix.o edouble.o
			gcc -O3 -g -o $(BINPATH)/drawphyl drawphyl.c file.o llist.o phyl.o grammar.o optimize.o col.o align.o matrix.o edouble.o -lm

$(BINPATH)/drawdot:     drawdot.c file.o llist.o
	                gcc -O3 -o $(BINPATH)/drawdot drawdot.c file.o llist.o

petfoldlibs.o:	petfoldlibs.c petfoldlibs.h
		gcc -O3 -g -c petfoldlibs.c

thermodynamic.o:	thermodynamic.c thermodynamic.h
			gcc -O3 -g -c thermodynamic.c -I$(HPATH) -lm

evolutionary.o:	evolutionary.c evolutionary.h
		gcc -O3 -g -c evolutionary.c

align.o:	align.c align.h file.o
		gcc -O3 -g -c align.c

col.o:		col.c col.h file.o
		gcc -O3 -g -c col.c

edouble.o:	edouble.c edouble.h file.o
		gcc -O3 -g -c edouble.c

file.o: 	file.c file.h
		gcc -O3 -g -c file.c

grammar.o:	grammar.c grammar.h
		gcc -O3 -g -c grammar.c

inout.o:	inout.c inout.h
		gcc -O3 -g -c inout.c
		 
llist.o: 	llist.c llist.h
		gcc -O3 -g -c llist.c

matrix.o:	matrix.c matrix.h file.o
		gcc -O3 -g -c matrix.c

newcolprob.o:	newcolprob.c newcolprob.h phyl.o llist.o col.o matrix.o edouble.o search.o
		gcc -O3 -g -c newcolprob.c

optimize.o:	optimize.c optimize.h
		gcc -O3 -g -c optimize.c

phyl.o:		phyl.c phyl.h
		gcc -O3 -g -c phyl.c

search.o:	search.c search.h
		gcc -O3 -g -c search.c

clean:
		rm -f *.o $(BINPATH)/PETfold $(BINPATH)/PETcofold $(BINPATH)/drawdot $(BINPATH)/drawphyl
