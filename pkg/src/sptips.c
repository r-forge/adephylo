/* Notes:
   these functions are used to find the shortest path between specified pairs of tips.
   The algorithm proceeds as follows:
   1) find the paths (pathA, pathB) to the root
   2) find the MRCA, defined as the first term of pathA in pathB (same as the converse)
   3) from A, go back to MRCA, adding crossed nodes to the result, not including the MRCA
   4) from B, go back to MRCA, adding crossed nodes to the result, not including the MRCA
   5) add the MRCA to the output
   6) return the output
*/


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Utils.h>
#include "adesub.h"



/* 
   =====================
   UTILITARY FUNCTIONS
   =====================
*/


/* 
   === REPLICATE %IN% / MATCH OPERATOR FOR INTEGERS === 
   == for internal use only ==
   - *b has to be a vector created using vecintalloc
   - returns 0 if there are no matches, and the index of the first match otherwise
*/
int intAinB(int a, int *b, int lengthB){
	if(lengthB == 0) return(0); /* avoid comparison with empty vector */

	int i=1;

	while(i < lengthB){
		if(b[i]==a) {
			return(i);
		} else {
			i++;
		}
	}
	return(0);
} /* intAinB */





/* 
   === REPLICATE SETDIFF MATCH OPERATOR FOR INTEGERS === 
   == for internal use only ==
   - *b has to be a vector created using vecintalloc
   - finds (possibly duplicated) elements of a not in b
*/
void intANotInB(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize){
	int i;

	/* a few checks */
	if(lengthA==0) return;
	if(lengthB==0){
		*resSize = 0; 
		return;
	}


	for(i=1; i<=lengthA; i++){
		if(intAinB(a[i], b, lengthB)==0){
			*resSize++;
			res[*resSize] = a[i];
		}
	}

} /* intANotInB */






/*
  === UNION OF TWO INTEGER VECTORS ===
  == for internal use only ==
  - a, b, and res have to be created by vecintalloc
  - returns unique(c(a,b))
*/
void unionInt(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize){
	if(lengthA==0 && lengthB && 0) {
		*res = 0;
		*resSize = 0;
		return;
	}

	int i, idx;

	res[1] = a[1]; /* initialization of temp results */
	*resSize = 1;

        /* For a */
	for(i=1;i<=lengthA;i++){
		idx = intAinB(a[i], res, *resSize); /* check if element is in res */
		if(idx==0) {
			*resSize++;
			res[*resSize] = a[i];
		}
	}

        /* For b */
	for(i=1;i<=lengthB;i++){
		idx = intAinB(b[i], res, *resSize); /* check if element is in res */
		if(idx==0) {
			*resSize++;
			res[*resSize] = b[i];
		}
	}
} /* unionInt */






/* 
   === INTERSECTION OF TWO INTEGER VECTORS ===
   == for internal use only ==
   - a, b, and res have to be created by vecintalloc
*/
void intersectInt(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize){
	if((lengthA * lengthB) ==0) {
		*res = 0;
		*resSize = 0;
		return;
	}
	int i, idx;

	*resSize = 0;

        /* store elements of a present in b */
	for(i=1;i<=lengthA;i++){
		idx = intAinB(a[i], b, lengthB) * intAinB(a[i], res, *resSize); /* a in b and not already in res */
		if(idx != 0) {
			*resSize++;
			res[*resSize] = a[i];
		}
	}
} /* intersectInt */





/* 
   === FIND THE PATH FROM A TIP TO THE ROOT === 
   == for internal use only ==
   - ances, desc and path must have been created using vecintalloc; their indices start at 1.
   - N is the number of edges in the tree, i.e. number of rows of $edge
*/
void pathTipToRoot(int tip, int *ances, int *desc, int N, int *res, int *resSize){
	int i, curNode=0, keepOn=1, nextNodeId;

	curNode = tip;

	while(keepOn==1){
		nextNodeId = intAinB(curNode, desc, N);
		if(nextNodeId != 0){
			*resSize++;
			res[*resSize] = ances[nextNodeId];
			curNode = ances[nextNodeId];
		} else {
			keepOn = 0;
		}
	}
} /* pathTipToRoot */





/*
  === FIND THE MRCA BETWEEN TWO TIPS ===
  == for internal use only ==
  - a and b are two tips
  - ances and desc must be created using vecintalloc
  - N is the number of edges
*/
int mrca2tips(int *ances, int*desc, int a, int b, int N){
	int *pathAroot, *pathBroot, *lengthPathA, *lengthPathB, i, res, idx;

	/* allocate memory */
	vecintalloc(&pathAroot, N);
	vecintalloc(&pathBroot, N);

	/* find paths to the root */
	pathTipToRoot(a, ances, desc, N, pathAroot, lengthPathA);
	pathTipToRoot(b, ances, desc, N, pathBroot, lengthPathB);

	/* initialization */
	i=1;
	idx = intAinB(pathAroot[1], pathBroot, *lengthPathA);
	
	while(idx==0){
		i++;
		idx = intAinB(pathAroot[i], pathBroot, *lengthPathA);
		if(i > *lengthPathA){ /* that would indicate an error */
			idx=0;
			printf("\n Likely error: no MRCA found between specified tips.");
				}
	}

	/* free memory */
	freeintvec(pathAroot);
	freeintvec(pathBroot);

	return(pathAroot[idx]);
} /* end mrca */






/*
  === FIND SHORTEST PATH BETWEEN TWO TIPS ===
  == for internal use only ==
  - ances and desc must be created using vecintalloc
  - N is the number of edges to represent the tree
*/
void sp2tips(int *ances, int *desc, int N, int tipA, int tipB, int *res, int *resSize){
	/* declarations */
	int k, *pathAroot, *pathBroot, *lengthPathA, *lengthPathB, myMrca;

	/* allocate memory */
	vecintalloc(&pathAroot, N);
	vecintalloc(&pathBroot, N);

	/* find paths to the root */
	pathTipToRoot(tipA, ances, desc, N, pathAroot, lengthPathA);
	pathTipToRoot(tipB, ances, desc, N, pathBroot, lengthPathB);

	/* find the MRCA between both tips */
	myMrca = mrca2tips(ances, desc, tipA, tipB, N);

	/* go back the paths and stop at MRCA (exclude MRCA) */
	/* for A */
	k = 1;
	*resSize = 0;
	while(pathAroot[k] != myMrca){ /* ### reprendre ici.*/
		*resSize++;
		res[*resSize] = pathAroot[k];
		k++;
	}

	/* for B */
	k = 1;
	while(pathBroot[k] != myMrca){ /* ### reprendre ici.*/
		*resSize++;
		res[*resSize] = pathBroot[k];
		k++;
	}

	/* add the MRCA */
	*resSize++;
	res[*resSize] = myMrca;


	/* free memory */
	freeintvec(pathAroot);
	freeintvec(pathBroot);

} /* end sp2tips */







/*
  === FIND SHORTEST PATH BETWEEN ALL PAIRS OF TIPS ===
  == for internal/external uses ==
  - all arguments are passed from R
  - N is the number of edges to represent the tree
  - nTips is the total number of tips in the tree
  - resSize is the total size of the output vector; it can't be known in advance, so a fake value has to be passed
  - resId indicates how the final result should be cut
*/
void spalltips(int *ances, int *desc, int N, int nTips, int *res, int *resId, int *resSize){
	/* declarations */
	int i, j, k, finalResId, *ancesLoc, *descLoc, *tempRes, *tempResSize, idPair;

	/* allocate memory for local variables */
	vecintalloc(&ancesLoc, N);
	vecintalloc(&descLoc, N);
	vecintalloc(&tempRes, N);


	/* create local vectors for ancestors and descendents */
	ancesLoc[0] = N;
	descLoc[0] = N;
	for(i=0; i<N; i++){
		ancesLoc[i+1] = ances[i];
		descLoc[i+1] = desc[i];
	}


	/* perform computations for all pairs of tips (indexed 'i,j') */
	*tempResSize=0;
	*resSize=0;
	finalResId=0; /* used to browse 'res' and 'resId' */
	idPair=0;
	
	for(i=0; i<=(N-2); i++){
		for(j=(i+1); j<=(N-1); j++){
			/* temp results are save in tempRes and tempResSize */
			idPair++;
			sp2tips(ancesLoc, descLoc, N, i+1, j+1, tempRes, tempResSize); /* i+1 and j+1 are tips id */

			/* copy temp results to returned results */
			*resSize += *tempResSize;
			for(k=1; k <= *tempResSize; k++){
				res[finalResId] = tempRes[k];
				resId[finalResId] = idPair;
				finalResId++;
			}
		}
	}

	/* free memory */
	freeintvec(ancesLoc);
	freeintvec(descLoc);
	freeintvec(tempRes);

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

.C("spalltips", as.integer(tre$edge[,1]), as.integer(tre$edge[,2]), nrow(tre$edge), as.integer(nTips(tre)), res, resId, resSize)
*/
