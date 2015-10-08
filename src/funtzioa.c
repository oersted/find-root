#include <math.h>

void f(int dim, double* x, double* fx) {
	fx[0] = -1 + pow(x[0], 2) + pow(x[1], 2) - cos(x[0]) * sin(x[1]);
	fx[1] = x[0] * x[1] - cos(x[1]) - sin(x[0]);
}

void jakobiarra(int dim, double* x, double* jx) {
	jx[0] = 2*x[0] + sin(x[0]) * sin(x[1]);
	jx[1] = 2*x[1] - cos(x[0]) * cos(x[1]);
	jx[2] = x[1] - cos(x[0]);
	jx[3] = x[0] + sin(x[1]);
}
