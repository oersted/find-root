void f(int dim, double* x, double* fx) {
	fx[0] = x[0]*x[0] - 1;
}

void jakobiarra(int dim, double* x, double* fx) {
	fx[0] = 2*x[0];
}

double norma(int dim, double* x) {
	return 0;
}
