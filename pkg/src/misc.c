#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>


/* interpolate values in a density */
/* x contais n values for which we want y=f(x) */
void predict_density(double *densx, double *densy, int *densn, double *x, double *y, int *n){
	int i, idx;

	for(i=0;i<*n;i++){
		idx=0;
		while(idx < *densn && x[i]>densx[idx]){
			idx++;
		}
		if(idx==0 || idx == *densn){
			y[i] = densy[idx]/2;
		} else {
			y[i] = (densy[idx] + densy[idx+1])/2;
		}
	}
}
