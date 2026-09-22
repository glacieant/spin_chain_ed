#include "randspin_utility.h"

/*Dependencies: randspin_utility.h*/

/*This is a program to find out the SZ spin states of a random Heisenberg spin chain*/

void randstate(int sitenum,int SZ,int *OneCount,int *LinLeast, int *LinMax, int *StateKets) {

	/*Parameters in the code*/

	int i, j, k, q, t, u;

	/***********************************************************************************************/

	/*Number of half site number, maximum attributable binary number for one
	half of the lattice*/

	int hfsitenum, halfmaxbinary;


	/***********************************************************************************************/

	/*Enumerating the states*/

	/*For StateKets[i], where i denotes the state number, the numerical value of StateKets[i]
	is equal to the binary number associated to the state.
	In the string of the binary number 1 stands for spin +half and 0 stands for spin -half*/

	hfsitenum = sitenum/2;


	halfmaxbinary = pow(2, hfsitenum) - 1;


	t = 0;
	k = 0;
	for (i = 0; i <= halfmaxbinary; i++) {
		if (OneCount[i] <= (hfsitenum + SZ)) {
			k = 0;
			for (j = 0; j <= halfmaxbinary; j++) {
				q = i + j*(halfmaxbinary + 1);
				if (OneCount[q] == (hfsitenum + SZ)) {
					LinLeast[i] = t;
					LinMax[j] = k;
					u = t + k;
					StateKets[u] = i + j*(halfmaxbinary + 1);
					k++;
				}
			}
			t = t + k;
		}
	}

}
