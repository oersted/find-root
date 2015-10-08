/*
 * exceptions.h
 *
 *  Created on: 26 Sep 2015
 *      Author: oersted
 */

/*
 * ERR_T, ERR_N can be defined by the user of this header before
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
#ifndef ERR_V
#define ERR_V ERR_4578919
#endif

// Declaration of the error value container
ERR_T ERR_V;

/*
 * If a function returns ERR_T, use this macro to return instead of the
 * standard C statement.
 *
 * It will ensure that the FINALLY section is executed before finishing the
 * execution of the function.
 *
 * To mantain abstraction, it is recommended to return only OK or ERR_V, but you
 * can return any ERR_T value.
 */
#define RETURN(value) \
	ERR_V = value; \
	goto FINALLY;

/*
 * If the given expression is false, the program will jump to the EXCEPT section
 * with the given error_value.
 */
#define TRY(error_value, expression) \
	if ((expression) == 0) { \
		ERR_V = error_value; \
		goto EXCEPT; \
	}

/*
 * If a function returns ERR_T it must have an EXCEPT section.
 *
 * Having and except section will ensure:
 * 	- The function will return OK if nothing goes wrong.
 * 	- The function will return ERR_V if something goes wrong.
 * 	- The FINALLY section will always be executed before the function has ended.
 *
 * The EXCEPT section should contain the content of a switch(ERR_V) block.
 */
#define EXCEPT(cases) \
	RETURN(OK) \
	\
	EXCEPT: \
		switch(ERR_V) { \
			cases \
		} \
		\
		RETURN(ERR_V) \

/*
 * If a function returns ERR_T it must have a FINALLY section, even if it's
 * empty.
 *
 * The code in the FINALLY section is ensured to be executed before the function
 * finishes.
 */
#define FINALLY(block) \
	FINALLY: \
		block \
		return ERR_V;
