#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#define MAX_BUF 256

/* hasierako datuak irakurtzen ditu. Lortzen badu, hau da, dimentsioak adierazi adina 
 datu irakurtzen baditu) irakurritako bektorearen dimentsioa itzuliko du
 bestela -1 edo 0 edo -(datukopurua)*/
int datuak_lortu(char* path, int *dptr, double **pptrptr, double *tptr) {
	FILE* f;
	char lerroa[MAX_BUF];
	int kont = 0;

	*dptr = 0;
	f = fopen(path, "r");
	if (f != NULL) {
		while (fgets(lerroa, MAX_BUF, f) != NULL) {
			if ((lerroa[0] != '#') && (lerroa[0] != '\0')) {
				if (strncmp(lerroa, "dimentsioa", 10) == 0) {
					*dptr = atoi(&(lerroa[11]));
					if (*dptr < 1) {
						printf("Dimentsioak 1 edo handiagoa izan behar du\n");
						return (-1);
					}
				}
			}
		}
		fclose(f);
		if (*dptr > 0)
			*pptrptr = (double*) malloc((*dptr) * sizeof(double));
		else {
			printf(
					"Dimentsioak agertu egin behar du eta 1 edo handiagoa izan behar du\n");
			return (*dptr);
		}
	} else {
		printf("ERROR: Ezin da fitxategia irakurri\n");
		return (-1);
	}
	//Berriz fitxategia ireki
	f = fopen(path, "r");
	if (f != NULL) {
		while (fgets(lerroa, MAX_BUF, f) != NULL) {
			if ((lerroa[0] != '#') && (lerroa[0] != '\0')
					&& (lerroa[0] != 'd')) {
				if (strncmp(lerroa, "tolerantzia", 11) == 0) {
					//toleratzia adierazi nahi du.
					*tptr = atof(&(lerroa[12]));
				} else {
					if (kont < *dptr) {
						(*pptrptr)[kont] = atof(lerroa);
						kont++;
					} else
						printf(
								"KONTUZ!: dimentsioa %d dela esan arren datu gehiago ageri dira fitxategian\n",
								*dptr);
				}
			}
		}
		// datu guztiak irakurri ditu eta bektorean jaso ditu
		if (kont < *dptr) {
			printf(
					"ERROREA!: dimentsia %d dela esan arren %d balio besterik ez dira ageri fitxategian\n",
					*dptr, kont);
			return (-kont);
		}
		fclose(f);
	} else {
		printf("ERROR: Ezin da fitxategia bigarren aldiz ireki...\n");
		return (-1);
	}
	return (*dptr);
}

