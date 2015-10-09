#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include "funtzioa.h"
#include "norma.h"

// Exception system configuration
#define ERR_T int
#define ERR_V ERRNO
#define OK 0
#include "exceptions.h"

// Input
#define MAX_BUF 256
#define MATCH(line, str) !strncmp(line, str, strlen(str))

// User options
#define DEFAULT_TOLERANCE 1.0e-12
#define DEFAULT_REL_TOL 1
#define DEFAULT_MAX_ZERO_DIST 1.0e-12
#define DEFAULT_USER_NORM 0
#define DEFAULT_MAX_ITER 25
#define DEFAULT_MAX_DIV_ITER 10
#define DEFAULT_JX_REUSE 5

struct options {
	double tolerance;
	char rel_tol;
	double max_zero_dist;
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

// Function declarations

/*
 * The result will be stored in x and x0 will contain the difference between x
 * and the result of the previous iteration.
 */
ERR_T findroot(int dim, double* x0, double* x, struct options* options,
		struct additional_data* data);

ERR_T norm(int dim, double* x, double* n, struct options* options);

ERR_T input_data(char* path, int* dim, double** x0, struct options* options);

void output_result(int dim, double* result, struct additional_data* data);

void output_vector(int dim, double* x);

void handler(const char* reason, const char* file, int line, int gsl_errno);

// Program

int main(int argc, char** argv) {
	int dim;
	char* path;
	double* x = NULL;
	double* x0 = NULL;
	struct options options;
	struct additional_data additional_data;

	// Save the conf file path
	TRY(1, argc == 2)
	path = argv[1];

	// Initialize options with default data
	options.tolerance = DEFAULT_TOLERANCE;
	options.rel_tol = DEFAULT_REL_TOL;
	options.max_zero_dist = DEFAULT_MAX_ZERO_DIST;
	options.user_norm = DEFAULT_USER_NORM;
	options.max_iter = DEFAULT_MAX_ITER;
	options.max_div_iter = DEFAULT_MAX_DIV_ITER;
	options.jx_reuse = DEFAULT_JX_REUSE;

	// Get conf file data
	TRY(2, input_data(path, &dim, &x0, &options) == OK)

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
			printf("[?] MEMORIA KOPURUA: %lu\n", dim * sizeof(double));
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

ERR_T findroot(int dim, double* x0, double* x, struct options* options,
		struct additional_data* data) {
	int s;
	unsigned int iter_count, iter_div_count, jx_reuse_count;
	double max_err, max_err_prev, zero_dist;
	double* fx = NULL;
	double* jx = NULL;
	clock_t begin, end;

	// INITIALIZATION

	TRY(1, (fx = (double*) malloc(dim * sizeof(double))) != NULL)
	TRY(2, (jx = (double*) malloc(dim * dim * sizeof(double))) != NULL)

	gsl_set_error_handler(&handler);

	gsl_vector_view x_gsl = gsl_vector_view_array(x, dim);
	gsl_vector_view x0_gsl = gsl_vector_view_array(x0, dim);
	gsl_vector_view fx_gsl = gsl_vector_view_array(fx, dim);
	gsl_matrix_view jx_gsl = gsl_matrix_view_array(jx, dim, dim);
	gsl_permutation* p = gsl_permutation_alloc(dim);

	// NEWTON-RAPHSON LOOP
	begin = clock();

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
	max_err = DBL_MAX;
	max_err_prev = DBL_MAX;
	jx_reuse_count = options->jx_reuse;

	do {
		f(dim, x, fx);

		// Reusage of the Jacobian matrix
		if (jx_reuse_count == options->jx_reuse) {
			jakobiarra(dim, x, jx);
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
			// zero_dist is reused to save memory
			TRY(7, norm(dim, x, &zero_dist, options) == OK)
			max_err /= zero_dist;
		}

		// Zero dist
		TRY(8, norm(dim, fx, &zero_dist, options) == OK)

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

	} while ((max_err > options->tolerance ||
			zero_dist > options->max_zero_dist) &&
			iter_count < options->max_iter &&
			iter_div_count < options->max_div_iter);

	end = clock();

	if (iter_count == options->max_iter)
		printf("[!] Iterazio mugara iritsi da.\n");
	if (iter_div_count == options->max_div_iter)
		printf("[!] Iterazio dibergente mugara iritsi da.\n");

	data->max_error = max_err;
	data->fx = fx;
	data->iter_count = iter_count;
	data->delta_t = (double) (end - begin) / (double) CLOCKS_PER_SEC;

	EXCEPT(
		case 1:
			fprintf(stderr,
					"[x] Ezin izan da memoria nahikoa erreserbatu.\n");
			printf("[?] MEMORIA KOPURUA: %lu\n", dim * sizeof(double));
			break;
		case 2:
			fprintf(stderr,
					"[x] Ezin izan da memoria nahikoa erreserbatu.\n");
			printf("[?] MEMORIA KOPURUA: %lun", dim * dim * sizeof(double));
			break;
		case 3:
			fprintf(stderr,
					"[x] Ezin izan da JX * x = FX ekuazio sistema "
					"linealaren LU deskonposaketa egin.\n");
			printf("[?] HURBILPEN PARTZIALA: ");
			output_vector(dim, x);
			printf("[?] ERRORE MAXIMOA: %.*g\n", DBL_DIG, max_err);
			printf("[?] F(X): ");
			output_vector(dim, fx);
			printf("[?] ITERAZIO_KOPURUA: %u\n", iter_count);
			break;
		case 4:
			fprintf(stderr,
					"[x] Ezin izan da JX * x = FX ekuazio sistema lineala "
					"ebatzi LU deskonposaketa erabiliz.\n");
			printf("[?] HURBILPEN PARTZIALA: ");
			output_vector(dim, x);
			printf("[?] ERRORE MAXIMOA: %.*g\n", DBL_DIG, max_err);
			printf("[?] F(X): ");
			output_vector(dim, fx);
			printf("[?] ITERAZIO_KOPURUA: %u\n", iter_count);
			break;
		case 5:
			fprintf(stderr,
					"[x] JX * x = FX ekuazio sistema linealaren emaitzaren "
					"norma kalkulatzean errore kritiko bat egon da.\n");
			printf("[?] HURBILPEN PARTZIALA: ");
			output_vector(dim, x);
			printf("[?] F(X): ");
			output_vector(dim, fx);
			printf("[?] ITERAZIO_KOPURUA: %u\n", iter_count);
			break;
		case 6:
			fprintf(stderr,
					"[x] Bektoreen arteko kenketa egitean errore kritiko "
					"bat egon da.");
			printf("[?] F(X): ");
			output_vector(dim, fx);
			printf("[?] ITERAZIO_KOPURUA: %u\n", iter_count);
			break;
		case 7:
			fprintf(stderr,
					"[x] x-ren norma kalkulatzean errore kritiko bat egon "
					"da.\n");
			printf("[?] HURBILPEN PARTZIALA: ");
			output_vector(dim, x);
			printf("[?] F(X): ");
			output_vector(dim, fx);
			printf("[?] ITERAZIO_KOPURUA: %u\n", iter_count);
			break;
		case 8:
			fprintf(stderr,
					"[x] f(x)-ren norma kalkulatzean errore kritiko bat egon "
					"da.\n");
			printf("[?] HURBILPEN PARTZIALA: ");
			output_vector(dim, x);
			printf("[?] ERRORE MAXIMOA: %.*g\n", DBL_DIG, max_err);
			printf("[?] F(X): ");
			output_vector(dim, fx);
			printf("[?] ITERAZIO_KOPURUA: %u\n", iter_count);
			break;
	)

	FINALLY(
		free(fx);
		free(jx);
		gsl_permutation_free(p);
	)
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
			printf("[?] BEKTOREA: ");
			output_vector(dim, x);
			break;
	)

	FINALLY()
}

ERR_T input_data(char* path, int* dim, double** x0, struct options* options) {
	FILE* f;
	char line[MAX_BUF];
	int count = 0;

	TRY(1, (f = fopen(path, "r")) != NULL)

	// Read dimension
	*dim = 0;
	while (fgets(line, MAX_BUF, f) != NULL) {
		if (line[0] != '#' && (line[0] != '\0')) {
			if (MATCH(line, "dimentsioa")) {
				TRY(2, sscanf(line, "dimentsioa %d", dim) == 1)
				TRY(3, *dim > 0)
				break;
			}
		}
	}

	TRY(4, *dim != 0)

	TRY(5, (*x0 = (double*) malloc((*dim) * sizeof(double))) != NULL)

	// Read x0 and optional parameters
	while (fgets(line, MAX_BUF, f) != NULL) {
		if (line[0] != '#' && (line[0] != '\0')) {
			if (MATCH(line, "tolerantzia")) {
				TRY(6,
					sscanf(line, "tolerantzia %lf", &(options->tolerance)) == 1)
				TRY(7, options->tolerance > 0.0)
			} else if (MATCH(line, "erlatiboa")) {
				if (MATCH(line, "erlatiboa bai")) {
					options->rel_tol = 1;
				} else if (MATCH(line, "erlatiboa ez")) {
					options->rel_tol = 0;
				} else {
					TRY(8, 0);
				}
			} else if (MATCH(line, "max_zero_distantzia")) {
				TRY(9,
					sscanf(line, "max_zero_distantzia %lf",
							&(options->max_zero_dist)) == 1)
				TRY(10, options->max_zero_dist > 0.0)
			} else if (MATCH(line, "ordezko_norma")) {
				if (MATCH(line, "ordezko_norma bai")) {
					options->user_norm = 1;
				} else if (MATCH(line, "ordezko_norma ez")) {
					options->user_norm = 0;
				} else {
					TRY(11, 0);
				}
			} else if (MATCH(line, "iterazio_maximoa")) {
				TRY(12,	sscanf(line, "iterazio_maximoa %u",
							&(options->max_iter)) == 1)
				TRY(13, options->max_iter != 0)
			} else if (MATCH(line, "dibergentzia_iterazio_maximoa")) {
				TRY(14, sscanf(line, "dibergentzia_iterazio_maximoa %u",
						&(options->max_div_iter)) == 1)
			} else if (MATCH(line, "jakobiar_berrerabilpena")) {
				TRY(15, sscanf(line, "jakobiar_berrerabilpena %u",
							&(options->jx_reuse)) == 1)
			} else {
				TRY(16, count < *dim)
				TRY(17, sscanf(line, "%lf", (*x0) + count) == 1)
				++count;
			}
		}
	}

	TRY(18, count == *dim)


	EXCEPT(
		case 1:
			fprintf(stdout,
					"[x] Ezin izan da konfigurazio fitxategia irakurri\n");
			break;
		case 2:
			fprintf(stdout,
					"[x] Sintaxi desegokia dimentsioa emateko lerroan\n");
			printf("[?] LERROA: %s", line);
			break;
		case 3:
			fprintf(stdout,
					"[x] Dimentsioak zero baina handiagoa izan behar du\n");
			printf("[?] LERROA: %s", line);
			break;
		case 4:
			fprintf(stdout, "[x] Ez da aurkitu dimentsioa.");
			break;
		case 5:
			fprintf(stdout,
					"[x] Ezin izan da X0-rentzat memoria erreserbatu\n.");
			printf("[?] MEMORIA KOPURUA: %lu", (*dim) * sizeof(double));
			break;
		case 6:
			fprintf(stdout,
					"[x] Sintaxi desegokia tolerantzia emateko lerroan\n.");
			printf("[?] LERROA: %s", line);
			break;
		case 7:
			fprintf(stdout,
					"[x] Tolerantziak zero baina handiagoa izan behar du.\n");
			printf("[?] LERROA: %s", line);
			break;
		case 8:
			fprintf(stdout,
					"[x] Sintaxi desegokia tolerantzia erlatiboa aukeratzeko "
					"lerroan. Aukera honek 'bai' eta 'ez' balioak hartu "
					"ditzake soilik.\n");
			printf("[?] LERROA: %s", line);
			break;
		case 9:
			fprintf(stdout,
					"[x] Sintaxi okerra zerotik distantzia maximoa emateko "
					"lerroan.\n");
			printf("[?] LERROA: %s", line);
			break;
		case 10:
			fprintf(stdout,
					"[x] Zerotik distantzia maximoak zero baina handiagoa izan "
					"behar du.\n");
			printf("[?] LERROA: %s", line);
			break;
		case 11:
			fprintf(stdout,
					"[x] Sintaxi desegokia ordezko norma aukeratzeko "
					"lerroan. Aukera honek 'bai' eta 'ez' balioak hartu "
					"ditzake soilik.\n");
			printf("[?] LERROA: %s", line);
			break;
		case 12:
			fprintf(stdout,
					"[x] Sintaxi okerra iterazio kopuru maximoa emateko "
					"lerroan. Kopuruak ezin du negatiboa izan.\n");
			printf("[?] LERROA: %s", line);
			break;
		case 13:
			fprintf(stdout,
					"[x] Iterazio kopuru maximoak ezin du zero izan.\n");
			printf("[?] LERROA: %s", line);
			break;
		case 14:
			fprintf(stdout,
					"[x] Sintaxi okerra iterazio dibergente kopuru maximoa "
					"emateko lerroan. Kopuruak ezin du negatiboa izan.\n");
			printf("[?] LERROA: %s", line);
			break;
		case 15:
			fprintf(stdout,
					"[x] Sintaxi okerra jakobiarraren berrerabilpen ratioa "
					"emateko lerroan. Ratioak ezin du negatiboa izan.\n");
			printf("[?] LERROA: %s", line);
			break;
		case 16:
			fprintf(stdout,
					"[x] X0-ren elementu kopurura dimentsioa baina handiagoa "
					"da.\n");
			printf("[?] DIMENTSIOA: %i", *dim);
			break;
		case 17:
			fprintf(stdout,
					"[x] Sintaxi desegokia X0-ren elementu bat emateko lerro "
					"batean.\n");
			printf("[?] LERROA: %s", line);
			break;
		case 18:
			fprintf(stdout,
					"[x] Dimentsioa eta X0-ren tamaina ez datoz bat.\n");
			printf("[?] DIMENTSIOA: %i", *dim);
			printf("[?] TAMAINA: %i", count);
			break;
		case 19:
			fprintf(stdout,
					"[x] Konfigurazio fitxategia ixtean errore kritiko bat "
					"egon da.\n");
			break;
	)

	FINALLY(
		if (f != NULL) {
			TRY(19, fclose(f) == OK);

			// On fail, so that it doesn't try again and again
			f = NULL;
		}
	)
}

void output_result(int dim, double* result, struct additional_data* data) {
	printf("\n");

	printf("Erroa: ");
	output_vector(dim, result);

	printf("Errore maximoa: %.*g\n", DBL_DIG, data->max_error);

	printf("Erroaren irudia: ");
	output_vector(dim, data->fx);

	printf("Iterazio kopurua: %u\n", data->iter_count);

	printf("Denbora: %.*g seg\n", DBL_DIG, data->delta_t);

	printf("\n");
}

void output_vector(int dim, double* x) {
	int i;

	printf("(");
	printf("%.*g", DBL_DIG, x[0]);
	for (i = 1; i < dim; ++i)
		printf(", %.*g", DBL_DIG, x[i]);
	printf(")\n");
}

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
