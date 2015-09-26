#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "datuak_lortu.h"
#include "funtzioak.h"

#define TRY(c, r) if (r == 0) {ERR = c; goto EXCEPT;}

int ERR;
size_t SIZE;
char HAVE_NORM = 0;

int norm(int dim, double* x, double* n) {
	if (HAVE_NORM) {
		norma(dim, x, norm);
	} else {
		gsl_vector_view x_gsl = gsl_vector_view_array(x, dim);
		*n = gsl_blas_dnrm2(&x_gsl.vector);
	}

	return 0;
}

int dist(int dim, double* x0, double* x1, double* d) {
	gsl_vector_view x0_gsl = gsl_vector_view_array(x0, dim);
	gsl_vector_view x1_gsl = gsl_vector_view_array(x1, dim);

	TRY(-1, gsl_vector_sub(&x0_gsl.vector, &x1_gsl.vector) == 0)
	TRY(-2, norm(dim, x0, d) == 0)

	return 0;

	EXCEPT:
	{
		switch(ERR) {
			case -1:
				fprintf(stderr, "[x] Distantzia kalkulatzeko bektoreen \
                                 arteko kenketa egitean errore kritiko \
                                 bat egon da.");
				break;
			case -2:
				fprintf(stderr, "[x] Distantzia kalkulatzeko norma kalku\
                                 latzean errore kritiko bat egon da.");
				break;
		}

		return ERR;
	}
}

void ebatzi(int dim, double* x0, double* fx, double* jx, double* x) {
	int s;
	
	gsl_vector_view x_gsl = gsl_vector_view_array(x, dim);
	gsl_vector_view x0_gsl = gsl_vector_view_array(x0, dim);
	gsl_vector_view fx_gsl = gsl_vector_view_array(fx, dim);
	gsl_matrix_view jx_gsl = gsl_matrix_view_array(fx, dim, dim);
	gsl_permutation* p = gsl_permutation_alloc(4);

	gsl_linalg_LU_decomp(&jx_gsl.matrix, p, &s);
	gsl_linalg_LU_solve(&jx_gsl.matrix, p, &fx_gsl.vector, &x_gsl.vector);
	
	gsl_vector_sub(&x_gsl.vector, &x0_gsl.vector);
	gsl_vector_scale(&x_gsl.vector, -1);

	return 0;
}

int findroot(int dim, double* x0, double* x, double tol) {
	double errorea;
	double* x1;
	double* x2;
	double* fx;
	double* jx;
	double* aux;

	gsl_set_error_handler_off();

	TRY(-1, (x1 = (double*) malloc(SIZE)) != NULL)
	TRY(-2, (x2 = (double*) malloc(SIZE)) != NULL)
	TRY(-3, (fx = (double*) malloc(SIZE)) != NULL)
	TRY(-4, (jx = (double*) malloc(SIZE)) != NULL)

	memcpy(x1, x0, SIZE);
	TRY(-5, norm(dim, x1, &errorea) == 0)

	while (err > tol) {
		f(dim, x1, fx);
		jakobitarra(dim, x1, jx);

		TRY(-6, ebatzi(dim, x1, fx, jx, x2) == 0)

		TRY(-7, dist(dim, x2, x1, &errorea) == 0)

		aux = x1;
		x1 = x2;
		x2 = aux;
	}

	return 0;

	EXCEPT:
	{
		switch(ERR) {
			case -1:
			case -2:
			case -3:
			case -4:
				fprintf(stderr, "[x] Ezin izan da memoria nahikoa erreserbatu.");
				break;
			case -5:
				fprintf(stderr, "[x] Norma kalkulatzean errore kritiko bat egon da.");
				break;
			case -7:
				fprintf(stderr, "[x] Distantzia kalkulatzean errore kritiko bat egon da.");
				break;
		}

		free(x1);
		free(x2);
		free(fx);
		free(jx);

		return ERR;
	}
}

int main(int argc, char** args) {
	int dim, i;	
	char* path_ptr;	
	double tol;
	double* x;
	double* x0;

	TRY(-1, argc == 2)
	path = args[1];

	TRY(-2, datuak_lortu(path, &dim, &x0, &tol) > 0)
	SIZE = dim * sizeof(double);

	x = (double*) malloc(SIZE);	
	TRY(-3, findroot(dim, x0, x, tol) == 0)

	printf("Emaitza: (");
	printf("%f", x[0]);
	for (i=1; i<dim; ++i)
		printf(", %lf", x[i]);
	printf(")");

	return 0;

	EXCEPT:
	{
		switch(ERR) {
			case -1:
				printf("Usage: findroot path\n");
				break;
			case -2:
				fprintf(stderr, "[x] Ezin izan da konfigurazio fitxategia ondo irakurri.\n");
				break;			
		}

		free(x);
		free(x0);

		return ERR;
	 }
}
