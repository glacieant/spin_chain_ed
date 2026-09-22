#include "randspin_utility.h"

/*This is a program to look at thermalization of heisenberg antiferromagnetic states*/

int main() {

	/***********************************************************************************************/

	FILE *fyl1, *fyl2;

	/***********************************************************************************************/

	/*Initial input parameters*/

	fyl1 = openfile("constants.txt", "r");

	/* Input order

	DISORDER *field or bond*
	DELTA *value*
	FIELD_DELTA *value*
	SZ *value*
	INIT_DISORDER_REALIZATION *value*
	MINSITE *value*
	MAXSITE *value*
	RATE *value*
	MEASURE *sz, bond or singlet*
	LOWNUM *value*
	HIGHNUM *value*
	PARTS *value*
	SBEGIN *value*
	SINGNUM *value*
	SEED *value*

	*/

	char DISORDER[1000], MEASURE[1000];

	readconstant(fyl1, "DISORDER", DISORDER);
	readconstant(fyl1, "MEASURE", MEASURE);

	/* Disorder parameter */

	const double DELTA = readdouble(fyl1, strcmp(DISORDER, "field") == 0 ? "FIELD_DELTA" : "DELTA");

	/*Number of lattice sites, half site number, overall SZ of the states,number of possible states, disorder realization
	for minimum number of lattice,minimum and maximum sitenum, rate in which disorder realizations would decrease, begining of
	S^2 sector(from lowest energy eigenstate) and number of singlet sector*/

	int sitenum, statenum, loweig, higheig;

	const int SZ = readint(fyl1, "SZ");
	const int INIT_DISORDER_REALIZATION = readint(fyl1, "INIT_DISORDER_REALIZATION");
	const int MINSITE = readint(fyl1, "MINSITE");
	const int MAXSITE = readint(fyl1, "MAXSITE");
	const int RATE = readint(fyl1, "RATE");
	const int LOWNUM = readint(fyl1, "LOWNUM");
	const int HIGHNUM = readint(fyl1, "HIGHNUM");
	const int PARTS = readint(fyl1, "PARTS");
	const int SBEGIN = readint(fyl1, "SBEGIN");
	const int SINGNUM = readint(fyl1, "SINGNUM");
	const unsigned long SEED = readlong(fyl1, "SEED");

	fclose(fyl1);


	/***********************************************************************************************/

	/*Setting clock to find execution time*/

	double total_time;
	clock_t start, end;

	start = clock();


	/*Array number for storing random numbers, Random number array, Random exchange couplings; Random Zeeman fields,
	state ket array, one-count array, two Lin list arrayes, eigenvalue array, eignevector array,
	temporary storage for the measure and errors computed and number of samples at each iteration, total number of samples  */

	int arraynum;
	double *RandomNumbers;
	double *ExchangeCoupling;
	double *ZeemanField;
	int *StateKets;
	int *OneCount;
	int *LinLeast;
	int *LinMax;
	double *EigenValues, *EigenVectors;
	double tempavg, merror;
	int sampnum;
	int totsamp;

	/*Number of disorder realization and realization loop, disorder average, standard error,
	number of lattice samples*/

	int disorder_realization, realization;
	double disorder_average, stderror_disorder;
	int latnum;
	int qpkk, bpcc ,stnm;


	/***********************************************************************************************/

	/*Computing the number of random numbers required to assign Zeeman fields or exchange couplings for all iterations*/

	latnum = (MAXSITE - MINSITE) / 2 + 1;
	arraynum = INIT_DISORDER_REALIZATION*
		((MINSITE - MAXSITE*1.0 / pow(RATE, latnum)) / (1.0 - 1.0 / RATE)
		+ 2.0*(1.0 / RATE - 1.0 / pow(RATE, latnum)) / pow(1 - 1.0 / RATE, 2));



	disorder_realization = INIT_DISORDER_REALIZATION;

	/***********************************************************************************************/

	/*Creating a big enough random number array to assign the Zeeman fields or exchange couplings for all the iterations*/

	RandomNumbers = allocate(arraynum, sizeof(double));
	if (strcmp(DISORDER, "field") == 0) {
		uniformdoublerand(&RandomNumbers[0], arraynum, -DELTA, DELTA, SEED);
	}
	else {
		powerdoublerand(&RandomNumbers[0], arraynum, DELTA, SEED);
	}


	/***********************************************************************************************/

	/*Creating an one-count table for binary integer strings encounterd in the problem */

	arraynum = pow(2,MAXSITE) - 1;
	OneCount = allocate((size_t)arraynum + 1, sizeof(int));
	for (qpkk = 0; qpkk <= arraynum ;qpkk++) {
		OneCount[qpkk] = fonecount(qpkk, MAXSITE);
	}

	/***********************************************************************************************/

	/*Evaluating the absolute spin difference average for different disorder realization*/

	bpcc = 0;


	for (sitenum = MINSITE; sitenum <= MAXSITE; sitenum += 2) {

		statenum = comb(sitenum, (sitenum / 2 - SZ));
		disorder_average = 0;
		stderror_disorder = 0;
		StateKets = allocate(statenum, sizeof(int));

		stnm = sitenum/2;
		stnm = pow(2, stnm) - 1;
		LinLeast = allocate(stnm + 1, sizeof(int));
		LinMax = allocate(stnm + 1, sizeof(int));

		randstate(sitenum,SZ,OneCount,LinLeast,LinMax,StateKets);

		loweig = statenum*LOWNUM/PARTS;
		higheig = statenum*HIGHNUM/PARTS;
		if (higheig > statenum - 1) {
			higheig = statenum - 1;
		}

		totsamp = 0;

		for (realization = 0; realization < disorder_realization; realization++) {

			ExchangeCoupling = allocate(sitenum, sizeof(double));
			ZeemanField = allocate(sitenum, sizeof(double));

			EigenValues = allocate(statenum, sizeof(double));
			EigenVectors = allocate((size_t)statenum*statenum, sizeof(double));


			/*Assign Zeeman fields; Assign Exchange couplings.*/

			for (qpkk = 0; qpkk<sitenum; qpkk++) {
				if (strcmp(DISORDER, "field") == 0) {
					ExchangeCoupling[qpkk] = 1.0;
					ZeemanField[qpkk] = RandomNumbers[bpcc];
				}
				else {
					ExchangeCoupling[qpkk] = RandomNumbers[bpcc];
				}
				bpcc++;
			}


			/*Creating and diagonalising the Hamiltonian of the problem */

			randomham(&ExchangeCoupling[0], &ZeemanField[0], sitenum, statenum, LinLeast, LinMax, StateKets, EigenValues, EigenVectors);

			/*Calculating absolute average spin diference*/

			if (strcmp(MEASURE, "sz") == 0) {
				tempavg = avgspingap(loweig, higheig, sitenum, statenum, EigenVectors, StateKets);
				merror = tempavg*tempavg;
				sampnum = 1;
			}
			else if (strcmp(MEASURE, "bond") == 0) {
				sampnum = localspindiff(loweig, higheig, sitenum, statenum, EigenVectors, StateKets, LinLeast, LinMax, &tempavg, &merror);
			}
			else {
				sampnum = localspindiff_singlet(loweig, higheig, SBEGIN, SINGNUM, sitenum, statenum, EigenVectors, StateKets, LinLeast, LinMax, &tempavg, &merror);
			}
			totsamp = totsamp + sampnum;

			disorder_average = disorder_average + tempavg;
			stderror_disorder = stderror_disorder + merror;
			free(ExchangeCoupling);
			free(ZeemanField);
			free(EigenValues);
			free(EigenVectors);
		}
		disorder_average = disorder_average / totsamp;
		stderror_disorder = stderror_disorder / totsamp;
		stderror_disorder = sqrt((stderror_disorder - disorder_average*disorder_average)/ totsamp);
		stderror_disorder = stderror_disorder / disorder_average;

		fyl2 = openfile("randoutput.txt", "a");
		end = clock();
		total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
		fprintf(fyl2, "%d\t%f\t%f\t%d\t%f\t%f\n", sitenum, log(disorder_average), stderror_disorder, disorder_realization ,DELTA ,total_time);
		fclose(fyl2);
		free(StateKets);
		free(LinLeast);
		free(LinMax);

		disorder_realization = disorder_realization / RATE;

	}



	/***********************************************************************************************/

	free(RandomNumbers);
	free(OneCount);




	end = clock();
	total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("Time taken: %f\n", total_time);

	return 0;

}
