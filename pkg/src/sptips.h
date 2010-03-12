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


int intAinB(int a, int *b, int lengthB);
void intANotInB(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize);
void unionInt(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize);
void intersectInt(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize);
void pathTipToRoot(int tip, int *ances, int *desc, int N, int *res, int *resSize);
int mrca2tips(int *ances, int*desc, int a, int b, int N);
void sp2tips(int *ances, int *desc, int N, int tipA, int tipB, int *res, int *resSize);
void spalltips(int *ances, int *desc, int *N, int *nTips, int *res, int *resId, int *resSize);
