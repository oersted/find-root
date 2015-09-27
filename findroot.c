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

/*
 * The result will be stored in x and x0 will contain the difference between x
 * and the result of the previous iteration.
 */
ERR_T findroot(int dim, double tol, double* x0, double* x) {
	int s;
	double errorea;
	double* fx;
	double* jx;
	double* aux;

	// INITIALIZATION

	// Initialize pointers with NULL to free safely
	fx = NULL;
	jx = NULL;

	TRY(-1, (fx = (double*) malloc(SIZE)) != NULL)
	TRY(-2, (jx = (double*) malloc(SIZE)) != NULL)

	gsl_set_error_handler_off();

	gsl_vector_view x_gsl = gsl_vector_view_array(x, dim);
	gsl_vector_view x0_gsl = gsl_vector_view_array(x0, dim);
	gsl_vector_view fx_gsl = gsl_vector_view_array(fx, dim);
	gsl_matrix_view jx_gsl = gsl_matrix_view_array(fx, dim, dim);
	gsl_permutation* p = gsl_permutation_alloc(4);

	// NEWTON-RAPHSON LOOP

	/*
	 * a, b and c are vectors of dim dimensions
	 * FX and JX are square matrixes of dim dimensions
	 *
	 * b = a - c
	 * WHERE JX * c = FX
	 *
	 * a beign the result of the previous iteration (or the initial point)
	 * b beign the result of the current iteration (or the final result)
	 */

	// x0 == a

	x = x0;
	TRY(-3, norm(dim, x, &errorea) == OK)

	// x == a
	// x0 == a

	while (errorea > tol) {
		f(dim, x, fx);
		jakobiarra(dim, x, jx);

		TRY(-4, gsl_linalg_LU_decomp(&jx_gsl.matrix, p, &s) == OK)
		TRY(-5,
			gsl_linalg_LU_solve(
				&jx_gsl.matrix, p, &fx_gsl.vector, &x0_gsl.vector) == OK
		)

		// x == a
		// x0 == c

		TRY(-6, norm(dim, x0, &errorea) == OK)

		TRY(-7, gsl_vector_sub(&x_gsl.vector, &x0_gsl.vector) == OK)

		// x == b
		// x0 == c
	}

	return OK;

	EXCEPT: {
		switch(ERR_V) {
			case -1:
			case -2:
				fprintf(stderr,
						"[x] Ezin izan da memoria nahikoa erreserbatu.\n");
				break;
			case -3:
				fprintf(stderr,
						"[x] Hasierako puntuaren norma kalkulatzean errore "
						"kritiko bat egon da.\n");
				break;
			case -4:
				fprintf(stderr,
						"[x] Ezin izan da JX * x = FX ekuazio sistema "
						"linealaren LU deskonposaketa egin.\n");
				break;
			case -5:
				fprintf(stderr,
						"[x] Ezin izan da JX * x = FX ekuazio sistema lineala "
						"ebatzi LU deskonposaketa erabiliz\n");
				break;
			case -6:
				fprintf(stderr,
						"[x] JX * x = FX ekuazio sistema linealaren emaitzaren "
						"norma kalkulatzean errore kritiko bat egon da.");
				break;
			case -7:
				fprintf(stderr,
						"[x] Bektoreen arteko kenketa egitean errore kritiko "
						"bat egon da");
				break;
		}

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

	// Initialize pointers with to free safely
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
	TRY(-4, findroot(dim, tol, x0, x) == OK)

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
