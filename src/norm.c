#include <math.h>

double norma(int dim, double* x) {
	int sum = 0;
	int i;
	for (i = 0; i < dim; ++i) {
		sum += fabs(x[i]);	
	}
	return sum;
}
