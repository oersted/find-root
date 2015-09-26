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
		TRY(-1, (*n = gsl_blas_dnrm2(&x_gsl.vector)) >= 0)
	}

	return OK;

	EXCEPT: {
		switch(ERR_V) {
			case -1:
				fprintf(stderr, "[x] Ezin izan da norma kalkulatu.\n");
				break;
		}

		return ERR_V;
	}
}

ERR_T dist(int dim, double* x0, double* x1, double* d) {
	gsl_vector_view x0_gsl = gsl_vector_view_array(x0, dim);
	gsl_vector_view x1_gsl = gsl_vector_view_array(x1, dim);

	TRY(-1, gsl_vector_sub(&x0_gsl.vector, &x1_gsl.vector) == OK)
	TRY(-2, norm(dim, x0, d) == OK)

	return OK;

	EXCEPT: {
		switch (ERR_V) {
			case -1:
				fprintf(stderr,
						"[x] Distantzia kalkulatzeko bektoreen arteko kenketa "
						"egitean errore kritiko bat egon da.\n");
				break;
			case -2:
				fprintf(stderr,
						"[x] Distantzia kalkulatzeko norma kalkulatzean errore "
						"kritiko bat egon da.\n");
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

	TRY(-1, gsl_linalg_LU_decomp(&jx_gsl.matrix, p, &s) == OK)
	TRY(-2,
		gsl_linalg_LU_solve(&jx_gsl.matrix,p,&fx_gsl.vector,&x_gsl.vector)== OK)

	TRY(-3, gsl_vector_sub(&x_gsl.vector, &x0_gsl.vector) == OK)
	TRY(-4, gsl_vector_scale(&x_gsl.vector, -1) == OK)

	return OK;

	EXCEPT: {
		switch(ERR_V) {
			case -1:
				fprintf(stderr, "[x] Ezin izan da LU deskonposaketa egin.\n");
				break;
			case -2:
				fprintf(stderr, "[x] Ezin izan da sistema lineala ebatzi.\n");
				break;
			case -3:
				fprintf(stderr, "[x] Ezin izan da bektore kenketa egin.\n");
				break;
			case -4:
				fprintf(stderr, "[x] Ezin izan da bektore eskalaketa egin.\n");
				break;
		}

		return ERR_V;
	}
}

ERR_T solve_1dim(int dim, double* x0, double* fx, double* jx, double* x) {

	TRY(-1, jx != 0);
	*x = *x0 - (*fx)/(*jx);

	return OK;

	EXCEPT:
		return ERR_V;
}

ERR_T findroot(int dim, double* x0, double* x, double tol) {
	double errorea;
	double* x1;
	double* fx;
	double* jx;
	double* aux;

	// Inicialize pointers with to free safely
	x1 = NULL;
	fx = NULL;
	jx = NULL;

	TRY(-1, (x1 = (double*) malloc(SIZE)) != NULL)
	TRY(-2, (fx = (double*) malloc(SIZE)) != NULL)
	TRY(-3, (jx = (double*) malloc(SIZE)) != NULL)

	gsl_set_error_handler_off();

	memcpy(x1, x0, SIZE);
	TRY(-4, norm(dim, x1, &errorea) == OK)

	while (errorea > tol) {
		f(dim, x1, fx);
		jakobiarra(dim, x1, jx);

		TRY(-5, jx != 0);
		errorea = (*fx)/(*jx);
		*x = *x1 - errorea;

		aux = x1;
		x1 = x;
		x = aux;
	}

	return OK;

	EXCEPT: {
		switch(ERR_V) {
			case -1:
			case -2:
			case -3:
				fprintf(stderr,
						"[x] Ezin izan da memoria nahikoa erreserbatu.\n");
				break;
			case -4:
				fprintf(stderr,
						"[x] Norma kalkulatzean errore kritiko bat egon da.\n");
				break;
			case -5:
				fprintf(stderr,
						"[x] Distantzia kalkulatzean errore kritiko bat egon "
						"da.\n");
				break;
		}

		free(x1);
		free(fx);
		free(jx);

		return ERR_V;
	}
}

int main(int argc, char** argv) {
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
	path = argv[1];

	// Get conf file data
	TRY(-2, datuak_lortu(path, &dim, &x0, &tol) > 0)
	SIZE = dim * sizeof(double);

	// Find root
	TRY(-3, (x = (double*) malloc(SIZE)) != NULL)
	TRY(-4, findroot(dim, x0, x, tol) == OK)

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
						"[x] Ezin izan da konfigurazio fitxategia ondo "
						"irakurri.\n");
				break;
			case -3:
				fprintf(stderr,
						"[x] Ezin izan da memoria nahikoa erreserbatu.\n");
				break;
			case -4:
				fprintf(stderr, "[x] Ezin izan da emaitza kalkulatu.\n");
				break;
		}

		free(x);
		free(x0);

		return ERR_V;
	}
}
