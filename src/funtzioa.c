#include <math.h>

void f(int dim, double* x, double* fx) {
	fx[0] = x[0] + x[1] - 10;
	fx[1] = pow(x[0], 2) - x[1] - 4;
}

void jakobiarra(int dim, double* x, double* jx) {
	jx[0] = 1;
	jx[1] = 1;
	jx[2] = 2*x[0];
	jx[3] = -1;
}
