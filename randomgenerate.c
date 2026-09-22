#include "randspin_utility.h"
#include "mt19937ar.h"

/*Dependencies: randspin_utility.h, mt19937ar.h*/

/* Generates an array of random doubles between llim and ulim with uniform distribution,
or with the power-law distribution of DELTA. Uses mt19937ar.
RandomNumbers is an array of doubles(defined outside) which stores the random no,
nrand is the length of RandomNumbers */


void uniformdoublerand(double *RandomNumbers, int nrand, double llim, double ulim, unsigned long SEED) {

	int i;

	init_genrand(SEED);

	for (i = 0; i < nrand; i++) {
		RandomNumbers[i] = llim + (ulim - llim)* genrand_real1();
	}

}


void powerdoublerand(double *RandomNumbers, int nrand, double DELTA, unsigned long SEED) {

	int i;
	double gamma;

	init_genrand(SEED);

	for (i = 0; i < nrand; i++) {
		gamma = -DELTA*log(1.0-genrand_real2());
		RandomNumbers[i] = exp(-gamma);
	}

}



/*************************************************************************************************************************/

