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
   =============================
   UTILITARY (INTERNAL) FUNCTIONS
   =============================
*/


/*
   === FIND THE LENGTH OF AN EDGE ===
   == for internal use only ==
  - the edge is identified by the descending node
  - ances, desc, and brlength must be created using vecintalloc
  - N is the number of edges to represent the tree
*/
double findedgelength(int *ances, int *desc, double *brlength, int N, int myNode){
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
int dist2tips(int *ances, int *desc, double *brlength, int N, int tipA, int tipB, int method){
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
		 intANotInB(path, myMrca, *lengthPath, 1, path, lengthPath);

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
		 res = static_cast<double> (*lengthPath);
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
   ==========================
     MAIN  (EXTERNAL) FUNCTION
   ==========================
*/



/*
  === FIND DISTANCES BETWEEN ALL PAIRS OF TIPS ===
  == for internal/external uses ==
  - all arguments are passed from R
  - N is the number of edges to represent the tree
  - nTips is the total number of tips in the tree
  - 'method' indicates the type of distance: 1) patristic 2) nNodes 3) Abouheif 4) sumDD
*/
void distalltips(int *ances, int *desc, double *brlength, int *N, int *nTips, double *res, int *resSize, int *method){
	/* declarations */
	int i, j, k, temp;
	int *ancesLoc, *descLoc, *brlengthLoc; /* must use dynamic allocation */


	/* check resSize */
	temp = (*nTips)*(nTips-1) / 2;
	if(*resSize !=  temp) {
		printf("\n Likely error in distalltips: resSize is %d, and should be %d.\n", *resSize, temp);
		return;
	}


	/* allocate memory for local variables */
	vecintalloc(&ancesLoc, *N);
	vecintalloc(&descLoc, *N);
	vecalloc(&brlengthLoc, *N);


	/* create local vectors for ancestors, descendents and branch lengths */
	ancesLoc[0] = *N;
	descLoc[0] = *N;
	brlengthLoc[0] = static_cast<double>(*N) ; /* conversion (casting) int->double*/
	for(i=0; i< *N; i++){
		ancesLoc[i+1] = ances[i];
		descLoc[i+1] = desc[i];
		brlengthLoc[i+1] = brlength[i];
	}


	/* perform computations for all pairs of tips (indexed 'i,j') */
	k = 0; /* used to browse 'res' and 'resId' */

	for(i=1; i<=(*nTips-1); i++){
		for(j=(i+1); j<=(*nTips); j++){
			res[k] = dist2tips(ancesLoc, descLoc, brlengthLoc, *N, i, j, *method);
			k++;
		}
	}

	/* free memory */
	freeintvec(ancesLoc);
	freeintvec(descLoc);
	freevec(brlengthLoc);

} /* end distalltips */




/* TESTING */
/*

library(adephylo)
tre=rtree(10)
plot(tre)
nodelabels()
tiplabels()

n <- as.integer(nTips(tre))
resSize=as.integer(n*(n-1)/2)
res <- integer(resSize)


# void distalltips(int *ances, int *desc, double *brlength, int *N, int *nTips, double *res, int *resSize, int *method){

## nb nodes
toto <- .C("distalltips", as.integer(tre$edge[,1]), as.integer(tre$edge[,2]), as.double(tre$edge.length), nrow(tre$edge), n, res, length(res), as.integer(2))
res <- toto[[6]]
res

## patristic
toto <- .C("distalltips", as.integer(tre$edge[,1]), as.integer(tre$edge[,2]), as.double(tre$edge.length), nrow(tre$edge), n, res, length(res), as.integer(1))
res <- toto[[6]]
res

## Abou
toto <- .C("distalltips", as.integer(tre$edge[,1]), as.integer(tre$edge[,2]), as.double(tre$edge.length), nrow(tre$edge), n, res, length(res), as.integer(3))
res <- toto[[6]]
res

## sumDD
toto <- .C("distalltips", as.integer(tre$edge[,1]), as.integer(tre$edge[,2]), as.double(tre$edge.length), nrow(tre$edge), n, res, length(res), as.integer(4))
res <- toto[[6]]
res

*/
