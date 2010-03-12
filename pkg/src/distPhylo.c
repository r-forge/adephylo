/*
  Coded by Thibaut Jombart (tjombart@imperial.ac.uk), March 2010.
  Distributed with the adephylo package for the R software.
  Licence: GPL >=2.

   Notes:
   These functions implement several different phylogenetic distances between all pairs of tips in a phylogeny.
   These functions require sptips.c and adesub.c.
*/


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Utils.h>
#include "adesub.h"
#include "sptips.h"


/* 
   =====================
   UTILITARY FUNCTIONS
   =====================
*/


/*
   === FIND THE LENGTH OF AN EDGE ===
   == for internal use only ==
  - the edge is identified by the descending node
  - ances, desc, and brlength must be created using vecintalloc
  - N is the number of edges to represent the tree
*/
double findedgelength(int *ances, int *desc, int *brlength, int N, int myNode){
	int posi=0;


	/* find the edge */
	posi = intAinB(myNode, desc, N);

	if(posi==0){
		printf("\n Likely error in findedgelength: edge not found");
		return(0.0);
	}

	/* return corresponding length */
	return(brlength[posi]);
} /* end findedgelength */





/*
   === FIND THE NUMBER OF DIRECT DESCENDENTS (DD) OF A NODE ===
   == for internal use only ==
   - ances, desc, and brlength must be created using vecintalloc
   - N is the number of edges to represent the tree
*/
int findNbDD(int *ances, int *desc, int N, int myNode){
	int i, nbDD=0;


	/* browse the ances vector */
	for(i=1; i<=N; i++){
		if(ances[i] == myNode) {
			nbDD++;
		}
	}

	if(nbDD==0){
		printf("\n Likely error in findNbDD: no direct descendent found.\n");
	}

	/* return corresponding length */
	return(nbDD);
} /* end findedgelength */







/*
  === DISTANCE(s) BETWEEN TWO TIPS ===
  == for internal use only ==
  - ances, desc, and brlength must be created using vecintalloc
  - N is the number of edges to represent the tree
  - 'method' indicates the type of distance: 1) patristic 2) nNodes 3) Abouheif 4) sumDD
  - edges are identified by their descending node
  - for patristic distances, the set of edge used is: {output of sp2tips} U {tipA, tipB} \ {output of mrca2tips}
  - for all others: {output of sp2tips}
*/
int dist2tips(int *ances, int *desc, int *brlength, int N, int tipA, int tipB, int *res, int *resSize, int method){
	/* declarations */
	int *path, *lengthPath, *myMrca;
	int i, res;


	/* allocate memory */
	vecintalloc(&path, N);
	lengthPath = (int *) calloc(1, sizeof(int));
	myMrca = (int *) calloc(1, sizeof(int));


	/* find the shortest path between the two tips */
	 sp2tips(ances, desc, N, tipA, tipB, path, lengthPath);


	 /* compute the distance */
	 switch( method )
	 {
	 case 1: /* patristic distance */
		 /* find the mrca */
		 *myMrca = 0;
		 *myMrca = mrca2tips(ances, desc, tipA, tipB, N);

		 /* remove mrca from the path */
		 intANotInB(path, myMrca, lengthPath, 1, path, lengthPath);

		 /* add tips to the path */
		 *lengthPath = *lengthPath + 1;
		 path[*lengthPath] = tipA;
		 *lengthPath = *lengthPath + 1;
		 path[*lengthPath] = tipB;

		 /* compute length */
		 res=0.0;
		 for(i=1; i<=*lengthPath; i++){
			 res += findedgelength(ances, desc, brlength, N, path[i]);
		 }
		 break;

	 case 2: /* number of nodes */
		 res = double(*lengthPath);
		 break;

	 case 3: /* prod DD (Abouheif) */
		 res=1.0;
		 for(i=1; i<=*lengthPath; i++){
			 res *= findNbDD(ances, desc, N, path[i]);
		 }
		 break;

	 case 4: /* sum DD */
		 res=0.0;
		 for(i=1; i<=*lengthPath; i++){
			 res += findNbDD(ances, desc, N, path[i]);
		 }
		 break;

	 default :
		 res=0.0;
		 printf("\n\n Likely error in dist2tips: unknown method (%d):", method);
		 break;
	 }

	/* free memory */
	freeintvec(path);
	free(lengthPath);
	free(myMrca);

	return(res)
} /* end dist2tips */







/*
  === ... BETWEEN ALL PAIRS OF TIPS ===
  == for internal/external uses ==
  - all arguments are passed from R
  - N is the number of edges to represent the tree
  - nTips is the total number of tips in the tree
  - resSize is the total size of the output vector; it can't be known in advance, so a fake value has to be passed
  - resId indicates how the final result should be cut
*/
void distalltips(int *ances, int *desc, int *N, int *nTips, int *res, int *resId, int *resSize){
	/* declarations */
	int i, j, k, m, idPair;
	int *ancesLoc, *descLoc, *tempRes, *tempResSize; /* must use dynamic allocation */

	/* allocate memory for local variables */
	vecintalloc(&ancesLoc, *N);
	vecintalloc(&descLoc, *N);
	vecintalloc(&tempRes, *N);
	tempResSize = (int *) calloc(1, sizeof(int));


	/* create local vectors for ancestors and descendents */
	ancesLoc[0] = *N;
	descLoc[0] = *N;
	for(i=0; i< *N; i++){
		ancesLoc[i+1] = ances[i];
		descLoc[i+1] = desc[i];
	}

	
	/* perform computations for all pairs of tips (indexed 'i,j') */
	*tempResSize = 0;
	*resSize = 0;
	m = 0; /* used to browse 'res' and 'resId' */
	idPair = 0;
	
	/* printf("\ngot to 1"); */
	/* debugging*/
/*	printf("\nancesLoc:\n");
	for(i=1; i<= *N;i++){
		printf(" %d", ancesLoc[i]);
	}

	printf("\ndesc:\n");
	for(i=1; i<= *N;i++){
		printf(" %d", descLoc[i]);
	}

	printf("\nN: %d", *N);
*/
	for(i=0; i<=(*nTips-2); i++){
		for(j=(i+1); j<=(*nTips-1); j++){
			/* temp results are save in tempRes and tempResSize */
			idPair++;
			sp2tips(ancesLoc, descLoc, *N, i+1, j+1, tempRes, tempResSize); /* i+1 and j+1 are tips id */

			/* copy temp results to returned results */
			*resSize = *resSize + *tempResSize;
			for(k=1; k <= *tempResSize; k++){
				res[m] = tempRes[k];
				resId[m] = idPair;
				m++;
			}
		}
	}
	/* printf("\ngot to 4"); */

	/* free memory */
	freeintvec(ancesLoc);
	freeintvec(descLoc);
	freeintvec(tempRes);
	free(tempResSize);

} /* end sptips */




/* TESTING */
/*

library(adephylo)
tre=rtree(10)
plot(tre)
nodelabels()
tiplabels()

res <- resId <- integer(1e5)
resSize=as.integer(1e5)

# void spalltips(int *ances, int *desc, int *N, int *nTips, int *res, int *resId, int *resSize){

toto <- .C("spalltips", as.integer(tre$edge[,1]), as.integer(tre$edge[,2]), nrow(tre$edge), as.integer(nTips(tre)), res, resId, resSize)
toto[[5]] <- toto[[5]][1:toto[[7]]]
toto[[6]] <- toto[[6]][1:toto[[7]]]

res <- split(toto[[5]], toto[[6]])
res

*/
