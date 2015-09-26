#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "datuak_lortu.h"
#include "funtzioak.h"

// Exception system configuration
#define ERR_T int
#define ERR_V ERRNO
#define OK 0
#include "exceptions.h"

// Globals

size_t SIZE;

// User options

char HAVE_NORM = 0;

ERR_T norm(int dim, double* x, double* n) {
	if (HAVE_NORM) {
		norma(dim, x, n);
	} else {
		gsl_vector_view x_gsl = gsl_vector_view_array(x, dim);
		TRY(-1, (*n = gsl_blas_dnrm2(&x_gsl.vector)) < 0)
	}

	return OK;

	EXCEPT: {
		switch(ERR_V) {
			case -1:
				fprintf(stderr, "[x] Ezin izan da norma kalkulatu.");
				break;
		}

		return ERR_V;
	}
}

ERR_T dist(int dim, double* x0, double* x1, double* d) {
	gsl_vector_view x0_gsl = gsl_vector_view_array(x0, dim);
	gsl_vector_view x1_gsl = gsl_vector_view_array(x1, dim);

	TRY(-1, gsl_vector_sub(&x0_gsl.vector, &x1_gsl.vector))
	TRY(-2, norm(dim, x0, d))

	return OK;

	EXCEPT: {
		switch (ERR_V) {
			case -1:
				fprintf(stderr,
						"[x] Distantzia kalkulatzeko bektoreen arteko kenketa \
						egitean errore kritiko bat egon da.");
				break;
			case -2:
				fprintf(stderr,
						"[x] Distantzia kalkulatzeko norma kalkulatzean errore \
						kritiko bat egon da.");
				break;
		}

		return ERR_V;
	}
}

ERR_T solve(int dim, double* x0, double* fx, double* jx, double* x) {
	int s;

	gsl_vector_view x_gsl = gsl_vector_view_array(x, dim);
	gsl_vector_view x0_gsl = gsl_vector_view_array(x0, dim);
	gsl_vector_view fx_gsl = gsl_vector_view_array(fx, dim);
	gsl_matrix_view jx_gsl = gsl_matrix_view_array(fx, dim, dim);
	gsl_permutation* p = gsl_permutation_alloc(4);

	TRY(-1, gsl_linalg_LU_decomp(&jx_gsl.matrix, p, &s))
	TRY(-2, gsl_linalg_LU_solve(&jx_gsl.matrix, p,&fx_gsl.vector,&x_gsl.vector))

	TRY(-3, gsl_vector_sub(&x_gsl.vector, &x0_gsl.vector))
	TRY(-4, gsl_vector_scale(&x_gsl.vector, -1))

	return OK;

	EXCEPT: {
		switch(ERR_V) {
			case -1:
				fprintf(stderr, "[x] Ezin izan da LU deskonposaketa egin.");
				break;
			case -2:
				fprintf(stderr, "[x] Ezin izan da sistema lineala ebatzi.");
				break;
			case -3:
				fprintf(stderr, "[x] Ezin izan da bektore kenketa egin.");
				break;
			case -4:
				fprintf(stderr, "[x] Ezin izan da bektore eskalaketa egin.");
				break;
		}

		return ERR_V;
	}
}

ERR_T findroot(int dim, double* x0, double* x, double tol) {
	double errorea;
	double* x1;
	double* x2;
	double* fx;
	double* jx;
	double* aux;

	// Inicialize pointers with to free safely
	x1 = NULL;
	x2 = NULL;
	fx = NULL;
	jx = NULL;

	TRY(-1, (x1 = (double*) malloc(SIZE)) != NULL)
	TRY(-2, (x2 = (double*) malloc(SIZE)) != NULL)
	TRY(-3, (fx = (double*) malloc(SIZE)) != NULL)
	TRY(-4, (jx = (double*) malloc(SIZE)) != NULL)

	gsl_set_error_handler_off();

	memcpy(x1, x0, SIZE);
	TRY(-5, norm(dim, x1, &errorea))

	while (errorea > tol) {
		f(dim, x1, fx);
		jakobiarra(dim, x1, jx);

		TRY(-6, solve(dim, x1, fx, jx, x2))

		TRY(-7, dist(dim, x2, x1, &errorea))

		aux = x1;
		x1 = x2;
		x2 = aux;
	}

	return OK;

	EXCEPT: {
		switch(ERR_V) {
			case -1:
			case -2:
			case -3:
			case -4:
				fprintf(stderr,
						"[x] Ezin izan da memoria nahikoa erreserbatu.");
				break;
			case -5:
				fprintf(stderr,
						"[x] Norma kalkulatzean errore kritiko bat egon da.");
				break;
			case -7:
				fprintf(stderr,
						"[x] Distantzia kalkulatzean errore kritiko bat egon \
						da.");
				break;
		}

		free(x1);
		free(x2);
		free(fx);
		free(jx);

		return ERR_V;
	}
}

int main(int argc, char** args) {
	int dim, i;
	char* path;
	double tol;
	double* x;
	double* x0;

	// Inicialize pointers with to free safely
	x = NULL;
	x0 = NULL;

	// Save the conf file path
	TRY(-1, argc == 2)
	path = args[1];

	// Get conf file data
	TRY(-2, datuak_lortu(path, &dim, &x0, &tol) > 0)
	SIZE = dim * sizeof(double);

	// Find root
	TRY(-3, (x = (double*) malloc(SIZE)) != NULL)
	TRY(-4, findroot(dim, x0, x, tol) == 0)

	// Output result
	printf("Emaitza: (");
	printf("%f", x[0]);
	for (i = 1; i < dim; ++i)
		printf(", %lf", x[i]);
	printf(")");

	return 0;

	EXCEPT: {
		switch(ERR_V) {
			case -1:
				printf("Usage: findroot path\n");
				break;
			case -2:
				fprintf(stderr,
						"[x] Ezin izan da konfigurazio fitxategia ondo \
						irakurri.\n");
				break;
			case -3:
				fprintf(stderr,
						"[x] Ezin izan da memoria nahikoa erreserbatu.");
				break;
			case -4:
				fprintf(stderr, "[x] Ezin izan da emaitza kalkulatu.");
				break;
		}

		free(x);
		free(x0);

		return ERR_V;
	}
}
