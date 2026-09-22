#include "randspin_utility.h"

/*Dependencies: randspin_utility.h*/

/*This is a program to find out and diagonalize a hamiltonian of a random Heisenberg spin chain*/

void randomham(double *ExchangeCoupling,double *ZeemanField,int sitenum,int statenum, int *LinLeast, int *LinMax, int *StateKets, double *EigenValues,double *EigenVectors) {

	/*Parameters in the code*/

	int i, j, k, m, p;

	/*Hamiltonian matrix*/

	double *SpinHamiltonian;



	/***********************************************************************************************/


	/*Computing the Hamiltonian matrix*/

	SpinHamiltonian = allocate((size_t)statenum*statenum, sizeof(double));

	for (i = 0; i < statenum; i++) {
		k = i + statenum*i;
		for (m = 0; m < sitenum; m++) {
			if (m < (sitenum-1)) {
				SpinHamiltonian[k] = SpinHamiltonian[k] + ExchangeCoupling[m]* pointspin(StateKets[i], m)*pointspin(StateKets[i], (m + 1));
				j = spinflip(sitenum,StateKets[i], m, m+1, LinLeast, LinMax);
				if (j != -1) {
					p = j + statenum*i;
					SpinHamiltonian[p] = SpinHamiltonian[p] + 0.5*ExchangeCoupling[m];
				}
			}
			else {
				SpinHamiltonian[k] = SpinHamiltonian[k] + ExchangeCoupling[m]* pointspin(StateKets[i], m)*pointspin(StateKets[i], 0);
				j = spinflip(sitenum,StateKets[i], m, 0, LinLeast, LinMax);
				if (j != -1) {
					p = j + statenum*i;
					SpinHamiltonian[p] = SpinHamiltonian[p] + 0.5*ExchangeCoupling[m];
				}
			}
		}
		for (m = 0; m < sitenum; m++) {
			SpinHamiltonian[k] = SpinHamiltonian[k] + ZeemanField[m]*pointspin(StateKets[i], m);
		}
	}


	/***********************************************************************************************/

	/*Diagonalise hamiltonian with eigenvalue and eigenvector*/

	diagonaliseA(SpinHamiltonian, statenum, EigenValues, EigenVectors);

	/***********************************************************************************************/

	free(SpinHamiltonian);
}
