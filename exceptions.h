/*
 * exceptions.h
 *
 *  Created on: 26 Sep 2015
 *      Author: oersted
 */

/*
 * ERR_T, ERR_N and OK can be defined by the user of this header before
 * including it.
 *
 * In fact, it is encouraged that they are defined explicitely in the src file.
 */

/*
 * The type of the global variable that will store the error values.
 */
#ifndef ERR_T
#define ERR_T int
#endif

/*
 * The name reserved to store the error values, this name shouln't be used in
 * the src file that includes this header.
 */
#ifndef ERR_N
#define ERR_V ERR_4578919
#endif

// Declaration of the error value container
ERR_T ERR_V;

/*
 * If a function returns this value, it will mean that it has been executed
 * without any issue.
 */
#ifndef OK
#define OK 0
#endif

/*
 * If the given expression returns a non OK value, the program will jump to the
 * EXCEPT label with the given error_value.
 */
#define TRY(error_value, expression) \
	if ((expression) != OK) { \
		ERR_V = error_value; \
		goto EXCEPT; \
	} \
