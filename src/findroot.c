#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include "datuak_lortu.h"
#include "funtzioa.h"
#include "norma.h"

// Exception system configuration
#define ERR_T int
#define ERR_V ERRNO
#define OK 0
#include "exceptions.h"

// User options
#define DEFAULT_TOLERANCE 1.0e-12
#define DEFAULT_REL_TOL 1
#define DEFAULT_USER_NORM 0
#define DEFAULT_MAX_ITER 25
#define DEFAULT_MAX_DIV_ITER 10
#define DEFAULT_JX_REUSE 0

struct options {
	double tolerance;
	char rel_tol;
	char user_norm;
	unsigned int max_iter;
	unsigned int max_div_iter;
	unsigned int jx_reuse;
};

// Result
struct additional_data {
	double max_error;
	double* fx;
	unsigned int iter_count;
	double delta_t;
};

void handler(const char* reason, const char* file, int line, int gsl_errno) {
	fprintf(stderr, "[x] ");

	switch(gsl_errno) {
		case GSL_EDOM:
			fprintf(stderr, "GSL DOMAIN ERROR: ");
			break;
		case GSL_ERANGE:
			fprintf(stderr, "GSL RANGE ERROR: ");
			break;
		case GSL_ENOMEM:
			fprintf(stderr, "GSL NO MEMORY AVAILABLE: ");
			break;
		case GSL_EINVAL:
			fprintf(stderr, "GSL INVALID ARGUMENT: ");
			break;
		default:
			fprintf(stderr, "GSL ERROR: ");
			break;
	}

	fprintf(stderr, "%s\n", reason);
}

void output_vector(int dim, double* x) {
	int i;

	printf("(");
	printf("%.*g", DBL_DIG, x[0]);
	for (i = 1; i < dim; ++i)
		printf(", %.*g", DBL_DIG, x[i]);
	printf(")\n");
}

void output_result(int dim, double* result, struct additional_data* data) {
	printf("Erroa: ");
	output_vector(dim, result);

	printf("Errore maximoa: %.*g\n", DBL_DIG, data->max_error);

	printf("F(Erroa): ");
	output_vector(dim, data->fx);

	printf("Iterazio kopurua: %u\n", data->iter_count);

	printf("Denbora: %.*g seg\n", DBL_DIG, data->delta_t);
}

ERR_T norm(int dim, double* x, double* n, struct options* options) {
	if (options->user_norm) {
		norma(dim, x, n);
	} else {
		gsl_vector_view x_gsl = gsl_vector_view_array(x, dim);
		TRY(1, (*n = gsl_blas_dnrm2(&x_gsl.vector)) >= 0.0)
	}

	EXCEPT(
		case 1:
			fprintf(stderr, "[x] Ezin izan da norma kalkulatu, posible da "
					"bektorea handiegia izatea.\n");
			break;
	)

	FINALLY()
}

/*
 * The result will be stored in x and x0 will contain the difference between x
 * and the result of the previous iteration.
 */
ERR_T findroot(int dim, double* x0, double* x, struct options* options,
		struct additional_data* data) {
	int s;
	unsigned int iter_count, iter_div_count, jx_reuse_count;
	double max_err, max_err_prev, x_norm;
	double* fx;
	double* jx;
	clock_t begin, end;

	begin = clock();

	// INITIALIZATION

	// Initialize pointers with NULL to free safely
	fx = NULL;
	jx = NULL;

	TRY(1, (fx = (double*) malloc(dim * sizeof(double))) != NULL)
	TRY(2, (jx = (double*) malloc(dim * dim * sizeof(double))) != NULL)

	gsl_set_error_handler(&handler);

	gsl_vector_view x_gsl = gsl_vector_view_array(x, dim);
	gsl_vector_view x0_gsl = gsl_vector_view_array(x0, dim);
	gsl_vector_view fx_gsl = gsl_vector_view_array(fx, dim);
	gsl_matrix_view jx_gsl = gsl_matrix_view_array(jx, dim, dim);
	gsl_permutation* p = gsl_permutation_alloc(dim);

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

	memcpy(x, x0, dim * sizeof(double));

	// x == a
	// x0 == a

	iter_count = 0;
	iter_div_count = 0;
	max_err_prev = DBL_MAX;
	jx_reuse_count = options->jx_reuse;

	do {
		f(dim, x, fx);
		jakobiarra(dim, x, jx);

		// fx == FX
		// jx == JX

		// Reusage of the Jacobian matrix
		if (jx_reuse_count == options->jx_reuse) {
			TRY(3, gsl_linalg_LU_decomp(&jx_gsl.matrix, p, &s) == OK)
			jx_reuse_count = 0;
		} else {
			++jx_reuse_count;
		}

		TRY(4,
			gsl_linalg_LU_solve(
				&jx_gsl.matrix, p, &fx_gsl.vector, &x0_gsl.vector) == OK
		)

		// x == a
		// x0 == c

		TRY(5, norm(dim, x0, &max_err, options) == OK)

		TRY(6, gsl_vector_sub(&x_gsl.vector, &x0_gsl.vector) == OK)

		// Relative error
		if (options->rel_tol) {
			TRY(7, norm(dim, x, &x_norm, options) == OK)
			max_err /= x_norm;
		}

		// x == b
		// x0 == c

		// UPDATES

		++iter_count;

		// Track divergence iterations
		if (max_err > max_err_prev) {
			++iter_div_count;
		} else {
			iter_div_count = 0;
		}
		max_err_prev = max_err;

	} while (max_err > options->tolerance &&
			iter_count < options->max_iter &&
			iter_div_count < options->max_div_iter);

	end = clock();

	data->max_error = max_err;
	data->fx = fx;
	data->iter_count = iter_count;
	data->delta_t = (double) (end - begin) / (double) CLOCKS_PER_SEC;

	EXCEPT(
		case 1:
		case 2:
			fprintf(stderr,
					"[x] Ezin izan da memoria nahikoa erreserbatu.\n");
			break;
		case 3:
			fprintf(stderr,
					"[x] Ezin izan da JX * x = FX ekuazio sistema "
					"linealaren LU deskonposaketa egin.\n");
			break;
		case 4:
			fprintf(stderr,
					"[x] Ezin izan da JX * x = FX ekuazio sistema lineala "
					"ebatzi LU deskonposaketa erabiliz.\n");
			break;
		case 5:
			fprintf(stderr,
					"[x] JX * x = FX ekuazio sistema linealaren emaitzaren "
					"norma kalkulatzean errore kritiko bat egon da.\n");
			break;
		case 6:
			fprintf(stderr,
					"[x] Bektoreen arteko kenketa egitean errore kritiko "
					"bat egon da.");
			break;
		case 7:
			fprintf(stderr,
					"[x] x-ren norma kalkulatzean errore kritiko bat egon "
					"da.\n");
			break;
	)

	FINALLY(
		free(fx);
		free(jx);
		gsl_permutation_free(p);
	)
}

int main(int argc, char** argv) {
	int dim;
	char* path;
	double* x;
	double* x0;
	struct options options;
	struct additional_data additional_data;

	// Initialize pointers with to free safely
	x = NULL;
	x0 = NULL;

	// Save the conf file path
	TRY(1, argc == 2)
	path = argv[1];

	// Initialize options with default data
	options.tolerance = DEFAULT_TOLERANCE;
	options.rel_tol = DEFAULT_REL_TOL;
	options.user_norm = DEFAULT_USER_NORM;
	options.max_iter = DEFAULT_MAX_ITER;
	options.max_div_iter = DEFAULT_MAX_DIV_ITER;
	options.jx_reuse = DEFAULT_JX_REUSE;

	// Get conf file data
	TRY(2, datuak_lortu(path, &dim, &x0, &(options.tolerance)) > 0)

	// Find root
	TRY(3, (x = (double*) malloc(dim * sizeof(double))) != NULL)
	TRY(4, findroot(dim, x0, x, &options, &additional_data) == OK)

	output_result(dim, x, &additional_data);

	RETURN(0)

	EXCEPT(
		case 1:
			printf("Usage: findroot path\n");
			break;
		case 2:
			fprintf(stderr,
					"[x] Ezin izan da konfigurazio fitxategia ondo "
					"irakurri.\n");
			break;
		case 3:
			fprintf(stderr,
					"[x] Ezin izan da memoria nahikoa erreserbatu.\n");
			break;
		case 4:
			fprintf(stderr, "[x] Ezin izan da emaitza kalkulatu.\n");
			break;
	)

	FINALLY(
		free(x);
		free(x0);
	)
}
