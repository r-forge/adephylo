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


/* === intersection of two integer vectors === */

/* 
=== replicate %in% / match operator for integers === 
- *b has to be a vector created using vecintalloc
- returns 0 if there are no matches, and the index of the first match otherwise
*/
int intAinB(int a, int *b, int lengthB){
	int i=1;

	while(i < lengthB){
		if(b[i]==a) {
			return(i);
		} else {
			i++;
		}
	}
	return(0);
}



/* === union of two integer vectors === */
void unionInt(int *a, int *b, int lengthA, int lengthB, int *res){
	int *temp, resSize, i, idx;
	vecintalloc(&temp, lengthA+lengthB); /* Initial size: as large as it could get */
	
	temp[1] = a[1]; /* initialization */
	resSize = 1;

        /* For a */
	for(i=1;i<=lengthA;i++){
		idx = intAinB(a[i], temp)
			if(idx != 0) {
				resSize++;
				temp[resSize] = a[i];
			}
	}

        /* For b */
	for(i=1;i<=lengthB;i++){
		idx = intAinB(b[i], temp)
			if(idx != 0) {
				resSize++;
				temp[resSize] = b[i];
			}
	}

	/* create final result */
	vecintalloc(&res, resSize);
	for(i=1; i<=resSize; i++){
		res[i] = temp[i];
	}

}



/* 
=== find the path from a tip to the root === 
- ances, desc and path must have been created using vecintalloc; their indices start at 1.
*/
void pathTipToRoot(int tip, int *ances, int *desc, int *path, int N){
	int i, curNode=0, keepOn=1, pathSize=0;

	vecintalloc(&path, N); /* Initial size: as large as it could get */

	curNode = tip;

	while(keepOn==1){
		nextNodeId = intAinB(curNode, desc);
		if(nextNodeId != 0){
			pathSize++;
			path[pathSize] = ances[nextNodeId];
			curNode = ances[nextNodeId];
		} else {
			keepOn = 0;
		}
	}
	
	/* resize path vector to minimal required size */ 
	int *temp=path;
	freeintvec(path);
	vecintalloc(&path, pathSize);
	for(i=1;i<=pathSize;i++){
		path[i] = temp[i];
	}
	path[0] = pathSize;

	freeintvec(temp);
}


/* === MRCA === */



/* === diff */

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
