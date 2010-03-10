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
   === replicate %in% / match operator for integers === 
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
   === replicate setdiff match operator for integers === 
   == for internal use only ==
   - *b has to be a vector created using vecintalloc
   - finds (possibly duplicated) elements of a not in b
*/
void intANotInB(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize){
	int i;

	/* a few checks */
	if(lengthA==0) return;
	if(lengthB==0){
		*resSize = 0; /* ### have to check that*/
		return;
	}


	for(i=1; i<=lengthA; i++){
		if(intAinB(a[i], b, lengthB)==0){
			*resSize++;
			res[resSize] = a[i];
		}
	}

} /* intANotInB */






/*
  === union of two integer vectors ===
  == for internal use only ==
  - a, b, and res have to be created by vecintalloc
  - returns unique(c(a,b))
*/
void unionInt(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize){
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
			res[resSize] = b[i];
		}
	}
} /* unionInt */






/* 
   === intersection of two integer vectors ===
   == for internal use only ==
   - a, b, and res have to be created by vecintalloc
*/
void intersectInt(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize){
	int resSize, i, idx;

	resSize = 0;

        /* store elements of a present in b */
	for(i=1;i<=lengthA;i++){
		idx = intAinB(a[i], b, lengthB) * intAinB(a[i], temp, *resSize); /* a in b and not already in temp */
		if(idx != 0) {
			*resSize++;
			res[*resSize] = a[i];
		}
	}
} /* intersectInt */





/* 
   === find the path from a tip to the root === 
   == for internal use only ==
   - ances, desc and path must have been created using vecintalloc; their indices start at 1.
   - N is the number of edges in the tree, i.e. number of rows of $edge
*/
void pathTipToRoot(int tip, int *ances, int *desc, int N, int *res, int *resSize){
	int i, curNode=0, keepOn=1, pathSize=0;

	curNode = tip;

	while(keepOn==1){
		nextNodeId = intAinB(curNode, desc, N);
		if(nextNodeId != 0){
			*resSize++;
			res[resSize] = ances[nextNodeId];
			curNode = ances[nextNodeId];
		} else {
			keepOn = 0;
		}
	}
} /* pathTipToRoot */





/*
  === find the MRCA between two tips ===
  == for internal use only ==
  - a and b are two tips
  - ances and desc must be created using vecintalloc
  - N is the number of edges
*/
int mrca2tips(int *ances, int*desc, int a, int b, int N){
	int *pathAroot, *pathBroot, *lengthPathA, *lengthPathB, i, res;

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
		if(i > lengthPathA){ /* that would indicate an error */
			idx=0;
			printf("\n Likely error: no MRCA found between specified tips.")
				}
	}

	/* free memory */
	freeintvec(pathAroot);
	freeintvec(pathBroot);

	return(pathAroot[idx]);
} /* end mrca */






/*
  === Find shortest path between two tips ===
  == for internal use only ==
  - ances and desc must be created using vecintalloc
  - N is the number of edges to represent the tree
*/
void sp2tips(int *ances, int *desc, int N, int tipA, int tipB, int N, int **res, inr *resSize){
	/* declarations */
	int k, *pathAroot, *pathBroot, *lengthPathA, *lengthPathB, myMrca;

	/* allocate memory */
	vecintalloc(&pathAroot, N);
	vecintalloc(&pathBroot, N);

	/* find paths to the root */
	pathTipToRoot(tipA, ances, desc, N, pathAroot, lengthPathA);
	pathTipToRoot(tipB, ances, desc, N, pathBroot, lengthPathB);

	/* find the MRCA between both tips */
	myMrca = mrca2tips(ances, desc, tipsA, tipsB, N);

	/* go back the paths and stop at MRCA (exclude MRCA) */
	/* for A */
	k = 1;
	*resSize = 0;
	while(pathAroot[k] != myMrca){ /* ### reprendre ici.*/
		resSize++;
		res[resSize] = pathAroot[k];
		k++;
	}

	/* for B */
	k = 1;
	while(pathBroot[k] != myMrca){ /* ### reprendre ici.*/
		resSize++;
		res[resSize] = pathBroot[k];
		k++;
	}

	/* add the MRCA */
	resSize++;
	res[resSize] = myMrca;


	/* free memory */
	freeintvec(pathAroot);
	freeintvec(pathBroot);

} /* end sp2tips */







/*
  === Find shortest path between two tips ===
  == for internal/external uses ==
  - ances, desc, tipsA, tipsB are passed from R
  - N is the number of edges to represent the tree
*/
void sptips(int *ances, int *desc, int N, int tipA, int tipB, int N, int **res, inr *resSize){
	/* declarations */
	int i, k, *pathAroot, *pathBroot, *lengthPathA, *lengthPathB, myMrca;

	/* allocate memory for local variables */
	vecintalloc(&ancesLoc, N);
	vecintalloc(&descLoc, N);


	/* create local vectors for ancestors and descendents */
	ancesLoc[0] = N;
	descLoc[0] = N;
	for(i=0; i<N; i++){
		ancesLoc[i+1] = ances[i];
		descLoc[i+1] = desc[i];
	}


	/* find paths to the root */
	pathTipToRoot(tipA, ances, desc, N, pathAroot, lengthPathA);
	pathTipToRoot(tipB, ances, desc, N, pathBroot, lengthPathB);


	/* find the MRCA between both tips */
	myMrca = mrca(ances, desc, tipsA, tipsB, N);


	/* go back the paths and stop at MRCA (exclude MRCA) */
	/* for A */
	k = 1;
	*resSize = 0;
	while(pathAroot[k] != myMrca){ /* ### reprendre ici.*/
		resSize++;
		res[resSize] = pathAroot[k];
		k++;
	}


	/* for B */
	k = 1;
	while(pathBroot[k] != myMrca){ /* ### reprendre ici.*/
		resSize++;
		res[resSize] = pathBroot[k];
		k++;
	}

	/* add the MRCA */
	resSize++;
	res[resSize] = myMrca;


	/* free memory */
	freeintvec(pathAroot);
	freeintvec(pathBroot);

} /* end sptips */


























void sharedAll(int *matAll, int *nRow, int *nCol, double *resVec)
{
/* Declare local C variables */
	int i, i1, i2, j, k, n, p, nbAll, **mat, temp;
	n = *nRow;
	p = *nCol;

	int nLoc=p/2;

/* Memory allocation for local C variables */

	tabintalloc(&mat, n, p); /* function from ade4 */

/* Local reconstruction of the matrix of alleles
   ! beware: this matrix is to be used from 1 to n and 1 to p included,
   and not from zero to n and p excluded as it is common in C */
	k = 0;
	for (j=1; j<=p; j++) {
		for (i=1; i<=n; i++) {
			mat[i][j] = matAll[k];
			k = k + 1;
		}
	}

/* == Main Computations: ==
   - i1, i2: indices of genotypes
   - j: index of allele
   - n: number of genotypes
   - p number of columns in mat (i.e. twice the number of loci)
   - each term in mat codes an allele (NAs are coded by 0)
*/

	k=0; /* counter used to browse resVec */
	for(i1=1; i1<=(n-1); i1++){
		for(i2=(i1+1); i2<=n; i2++){
			/* Used for debugging
			   printf("\n\n debug: ## %d-%d ##",i1,i2);
			*/

			resVec[k] = 0.0; /* Initialization of the result */
			nbAll = 0; /* counts the number of types alleles */
			for(j=1; j<=nLoc; j++){
				/* Used for debugging
				   printf("\n debug: j=%d",j);
				   printf("\n debug: mat[i1,j]=%d",mat[i1][j]);
				   printf("\n debug: mat[i1,j]=%d",mat[i1][j+nLoc]);
				   printf("\n debug: mat[i2,j]=%d",mat[i2][j]);
				   printf("\n debug: mat[i2,j+nLoc]=%d",mat[i2][j+nLoc]);
				*/
				if(mat[i1][j] != 0 && mat[i1][j+nLoc] != 0 && 
				   mat[i2][j] != 0 && mat[i2][j+nLoc] != 0){
					/* Used for debugging 
					   printf("\n debug: alleles are typed");
					*/
					nbAll+=2;
					/* Used for debugging
					   printf("\n debug: nbAll=%d", nbAll);
					*/
					/* Compare alleles: 
					   -> either both alleles are in common, 
					   -> or no allele are common, 
					   -> or there is one common allele */
					/* both alleles common */
					if((mat[i1][j] == mat[i2][j] 
					    && mat[i1][j+nLoc] == mat[i2][j+nLoc])
					   || (mat[i1][j] == mat[i2][j+nLoc]
					       && mat[i1][j+nLoc] == mat[i2][j])){
						resVec[k] += 2.0;
					} else if(!( /* if not 'all alleles differe' */
							  mat[i1][j] != mat[i2][j] 
							  && mat[i1][j] != mat[i2][j+nLoc]
							  && mat[i1][j+nLoc] != mat[i2][j] 
							  && mat[i1][j+nLoc] != mat[i2][j+nLoc])
						) resVec[k]++;

				} /* end if */
			} /* end for j in 1:(nLoc) */

			/* Divide the number of shared alleles by the number of typed alleles */
			if(nbAll > 0) resVec[k] = resVec[k]/nbAll;
			/*printf("\n debug: resVec[i1,i2]/nbAll (%d,%d)=# %f #", i1,i2,resVec[k]);*/
			k++;

		} /* end for i2 */
	} /* end for i1*/

	/* Free allocated memory */
	freeinttab(mat);

} /* end sharedAll */
