#include "randspin_utility.h"

/*Dependencies: randspin_utility.h*/

/***********************************************************************************************/

/* Calculates C(x,y)*/


int comb(int x, int y) {

	int i , r;
	long double num;

	r = (int)(fmax(y,x-y));

	num = 1.0;
	for (i = 0; i < r; i++) {
		num = num*(x-i)/(r-i);
	}
	return (int)(floor(num+0.5));
}

/***********************************************************************************************/

/*Counting ones of a binary number x of size y */


int fonecount(int x, int y) {

	int onenum, p;

	onenum = 0;

	for (p = 0; p < y; p++)	{
		if ((x & (1 << p)) != 0) {
			onenum++;
		}
	}

	return onenum;
}

/***********************************************************************************************/


/*Subroutine to find out the spin at y'th lattice site when the binary number corresponding to the state is x*/



double pointspin(int x, int y) {

	int spinmask1;

	spinmask1 = pow(2, y);
	x = x & spinmask1;
	if (x == spinmask1) {
		return 0.5;
	}
	else {
		return -0.5;
	}
}


/***********************************************************************************************/

/* Function to flip any two arbitrary digit of a binary string */


int spinflip(int sitenum,int x, int j, int k, int *LinLeast,int *LinMax) {

	int r, s, bsitenum, binary1 , binary2, binary3;
	bsitenum = sitenum/2;
	binary1 = pow(2, bsitenum) - 1;
	binary2 = binary1 + 1;
	binary3 = pow(2, sitenum) - 1 - binary1;


	int y, spinmask3;
	spinmask3 = pow(2, j) + pow(2, k);
	y = x ^ spinmask3;
	if (pointspin(y,j) != pointspin(y,k)) {
		r = y&binary1;
		r = LinLeast[r];
		s = y&binary3;
		s = s / binary2;
		s = LinMax[s];
		return r + s;
	}
	else {
		return -1;
	}	
}



/***********************************************************************************************/


/*calling dsyev subroutine to diagonalize matrix SpinHamiltonian of dimension y. The array EigenValues[] would
contain the eigenvalues and the array EigenVectors[] would contain the eigenvectors*/


void diagonaliseA(double *SpinHamiltonian, int y, double *EigenValues, double *EigenVectors) {


	double *Work = NULL;
	int lwork, info;

	Work = allocate(1, sizeof(double));
	lwork = -1;
	dsyev_("V", "U", &y, SpinHamiltonian, &y, EigenValues, Work, &lwork, &info);
	lwork = Work[0];
	free(Work);
	Work = allocate(lwork, sizeof(double));
	dsyev_("V", "U", &y, SpinHamiltonian, &y, EigenValues, Work, &lwork, &info);

	int i;

	for (i = 0; i < y*y; i++) {
		EigenVectors[i] = SpinHamiltonian[i];
	}
	free(Work);


}

/***********************************************************************************************/

double spinsector(int sitenum, int statenum, int eigenum, double *EigenVectors, int *LinLeast, int *LinMax , int *StateKets) {

	int j, k, m, n, ti;
	double kappa;

	kappa = sitenum * 3.0 / 4;
	for (m = 0; m < statenum; m++) {			
		ti = m + statenum*eigenum;
		for (j = 0; j < sitenum; j++) {
			for (k = 0; k < j; k++) {
				kappa = kappa + 2.0 * EigenVectors[ti] * EigenVectors[ti] * pointspin(StateKets[m], j)* pointspin(StateKets[m], k);
				n = spinflip(sitenum, StateKets[m], j, k, LinLeast, LinMax);
				if (n != -1) {
					n = n + statenum*eigenum; 
					kappa = kappa + EigenVectors[n]*EigenVectors[ti];
				}										
			}
		}
	}

	return kappa;
}




/***********************************************************************************************/

/*   Calculates |<S_i.S_{i+1}>_{k+1}-<S_i.S_{i+1}>_k| for site i and state k. Sums over i and over k from loweig to higheig.
 Over bonds and neighbouring eigenstates: measure = sum of the differences, merror = sum of their squares.
 sitenum is no. of sites; statenum is no. of total states; StateKets is array which gives spin config (binary) corr. to state no.
 EigenVectors is the array of eigenvectors obtained from diagonalization. See randstate.c   */


/* This subroutine does not check whether the eigenstates are singlets of total spin and returns the number of samples in the average*/

int localspindiff(int loweig, int higheig, int sitenum, int statenum, double *EigenVectors, int *StateKets, int *LinLeast, int *LinMax , double *measure, double *merror) {

/*declaring integers and doubles for later use*/

	int i, j, k, l, m, n, p, q, bsitenum, norm;

	double agap;
	double bgap;
	double cgap;
	agap = 0.0;
	cgap = 0.0;
	norm = 0;
	bsitenum = sitenum -1;

	/*************************************************************************************************/

	/*The measure |<n|Si.Si+1|n>-<n+1|Si.Si+1|n+1>| is implemented. Here |n> is the  eigenkets and i is the
	site number. Summation over sitenumber is also done for averaging purposes.*/
	
	for (i = loweig; i < higheig; i++) {
		for (m = 0; m < bsitenum; m++) {
			bgap = 0.0;
			for (j = 0; j < statenum; j++) {
				k = j + statenum*i;
				l = j + statenum*(i+1);

				/* the Sz.Sz part is implemented */

				bgap = bgap +
					(EigenVectors[k] * EigenVectors[k] - EigenVectors[l] * EigenVectors[l])*
					(pointspin(StateKets[j], m)*pointspin(StateKets[j], m + 1));

				/* the (S+.S-+S-.S+)/2 part is implemented. Given by spin-flip*/

				n = spinflip(sitenum, StateKets[j], m, m + 1, LinLeast, LinMax);
				if (n != -1) {
					p = n + statenum*i;
					q = n + statenum*(i+1);
					bgap = bgap + 0.5*(EigenVectors[k] * EigenVectors[p] - EigenVectors[l] * EigenVectors[q]);
				}
			}
			agap = agap + fabs(bgap);
			cgap = cgap + fabs(bgap)*fabs(bgap);
			norm++;
		}
	}
	/*Same measure |<n|Si.Si+1|n>-<n+1|Si.Si+1|n+1>| implemented where i and i+1 correspond the i'th(final) lattice site
	and zeroeth lattice site*/

	for (i = loweig; i < higheig; i++) {
		bgap = 0.0;
		for (j = 0; j < statenum; j++) {
			k = j + statenum*i;
			l = j + statenum*(i+1);
			bgap = bgap +
				(EigenVectors[k] * EigenVectors[k] - EigenVectors[l] * EigenVectors[l])*
				(pointspin(StateKets[j], bsitenum)*pointspin(StateKets[j], 0));
			n = spinflip(sitenum, StateKets[j], bsitenum, 0, LinLeast, LinMax);
			if (n != -1) {
				p = n + statenum*i;
				q = n + statenum*(i+1);
				bgap = bgap + 0.5*(EigenVectors[k] * EigenVectors[p] - EigenVectors[l] * EigenVectors[q]);
			}
		}
		agap = agap + fabs(bgap);
		cgap = cgap + fabs(bgap)*fabs(bgap);
		norm++;
	}

	(*measure) = agap;
	(*merror) = cgap;
	
	return norm;
}

/***********************************************************************************************/


/*   Calculates |<S_i.S_{i+1}>_{k+1}-<S_i.S_{i+1}>_k| for site i and state k. Sums over i and over k from loweig to higheig.
 Over bonds and neighbouring eigenstates: measure = sum of the differences, merror = sum of their squares.
 sitenum is no. of sites; statenum is no. of total states; StateKets is array which gives spin config (binary) corr. to state no.
 EigenVectors is the array of eigenvectors obtained from diagonalization. See randstate.c   */


/* This subroutine only averages over a maximum of the lowest SINGNUM singlet eigenstates starting from the SBEGIN'th one in the
range and returns the number of samples in the average */


int localspindiff_singlet(int loweig, int higheig, int SBEGIN, int SINGNUM, int sitenum, int statenum, double *EigenVectors, int *StateKets, int *LinLeast, int *LinMax, double *measure, double *merror) {

	/*declaring integers and doubles for later use*/
	
	int i, j, k, l, m, n, p, q, p1, p2, bsitenum, singletcount, startcount, norm;

	double agap;
	double bgap;
	double cgap;
	agap = 0.0;
	cgap = 0.0;
	norm = 0;
	bsitenum = sitenum -1;

	/*************************************************************************************************/


	double ACC = 0.000001;

	p1 = -1;
	startcount = 0;
	for (i = loweig; i <= higheig; i++) {
		if (fabs(spinsector(sitenum, statenum, i, EigenVectors, LinLeast, LinMax, StateKets)) <= ACC) {
			startcount++;
			if (startcount == SBEGIN) {			
				p1 = i;
				break;
			}
		}
	}

	if(p1 < 0){

		/* No singlet sector found in the specified range */

		(*measure) = 0;
		(*merror) = 0;
	
		return norm;
	}

	p2=-1;

	singletcount = 1;

	for(i=p1+1;i<=higheig;i++) {

		if (fabs(spinsector(sitenum, statenum, i, EigenVectors, LinLeast, LinMax, StateKets)) <= ACC) {

			p2 = i;

		/*The measure |<n|Si.Si+1|n>-<n+1|Si.Si+1|n+1>| is implemented. Here |n> is the  eigenkets and i is the
		site number. Summation over sitenumber is also done for averaging purposes.*/

			for (m = 0; m < bsitenum; m++) {
				bgap = 0.0;
				for (j = 0; j < statenum; j++) {
					k = j + statenum*p1;
					l = j + statenum*p2;

					/* the Sz.Sz part is implemented */

					bgap = bgap +
						(EigenVectors[k] * EigenVectors[k] - EigenVectors[l] * EigenVectors[l])*
						(pointspin(StateKets[j], m)*pointspin(StateKets[j], m + 1));

					/* the (S+.S-+S-.S+)/2 part is implemented. Given by spin-flip*/
	
					n = spinflip(sitenum, StateKets[j], m, m + 1, LinLeast, LinMax);
					if (n != -1) {
						p = n + statenum*p1;
						q = n + statenum*p2;
						bgap = bgap + 0.5*(EigenVectors[k] * EigenVectors[p] - EigenVectors[l] * EigenVectors[q]);
					}
				}
				agap = agap + fabs(bgap);
				cgap = cgap + fabs(bgap)*fabs(bgap);
				norm++;
			}

			/*The measure |<n|Si.Si+1|n>-<n+1|Si.Si+1|n+1>| , implemented for the periodic boundary part.*/

			bgap = 0.0;
			for (j = 0; j < statenum; j++) {
				k = j + statenum*p1;
				l = j + statenum*p2;

				/* the Sz.Sz part is implemented */
		
				bgap = bgap +
					(EigenVectors[k] * EigenVectors[k] - EigenVectors[l] * EigenVectors[l])*
					(pointspin(StateKets[j], bsitenum)*pointspin(StateKets[j], 0));

				/* the (S+.S-+S-.S+)/2 part is implemented. Given by spin-flip*/
	
				n = spinflip(sitenum, StateKets[j], bsitenum, 0, LinLeast, LinMax);
				if (n != -1) {
					p = n + statenum*p1;
					q = n + statenum*p2;
					bgap = bgap + 0.5*(EigenVectors[k] * EigenVectors[p] - EigenVectors[l] * EigenVectors[q]);
				}
			}
			agap = agap + fabs(bgap);
			cgap = cgap + fabs(bgap)*fabs(bgap);
			norm++;
			singletcount++;
			p1 = p2;
		}

		if (singletcount == SINGNUM) {
			break;
		}
	}

	(*measure) = agap;
	(*merror) = cgap;

	return norm;
}

/***********************************************************************************************/

/*Calculating absolute average spin diference*/

double avgspingap(int loweig, int higheig, int sitenum, int statenum, double *EigenVectors, int *StateKets) {

	int i, j, k, l, m, norm;

	double agap;
	double bgap;
	agap = 0.0;
	norm = 0;

	for (m = 0; m < sitenum; m++) {
		for (i = loweig; i < higheig; i++) {
			bgap = 0.0;
			for (j = 0; j < statenum; j++) {
				k = j + statenum*i;
				l = j + statenum*(i + 1);
				bgap = bgap + (EigenVectors[k]*EigenVectors[k] - EigenVectors[l]*EigenVectors[l])*pointspin(StateKets[j], m);
			}
			agap=agap+fabs(bgap);
			norm++;
		}
	}
	agap = agap / norm;

	return agap;
}

/***********************************************************************************************/

void *allocate(size_t num, size_t size) {

	void *Block;

	Block = calloc(num, size);
	if (Block == NULL) {
		fprintf(stderr, "allocation error: %zu x %zu bytes\n", num, size);
		exit(EXIT_FAILURE);
	}
	return Block;
}

FILE *openfile(const char *name, const char *mode) {

	FILE *fyl;

	fyl = fopen(name, mode);
	if (fyl == NULL) {
		fprintf(stderr, "file error: %s\n", name);
		exit(EXIT_FAILURE);
	}
	return fyl;
}

void readconstant(FILE *fyl, const char *KEY, char *Value) {

	char Line[1000];
	char *eq, *start, *end;
	size_t keylength;

	keylength = strlen(KEY);
	rewind(fyl);
	while (fgets(Line, sizeof(Line), fyl) != NULL) {
		eq = strchr(Line, '=');
		if (eq != NULL && (size_t)(eq - Line) == keylength && strncmp(Line, KEY, keylength) == 0) {
			start = eq + 1;
			end = start + strcspn(start, "\r\n");
			if (*start == '"' && end > start + 1 && *(end - 1) == '"') {
				start++;
				end--;
			}
			memcpy(Value, start, end - start);
			Value[end - start] = '\0';
			return;
		}
	}
	fprintf(stderr, "constant error: %s\n", KEY);
	exit(EXIT_FAILURE);
}

int readint(FILE *fyl, const char *KEY) {

	char Value[1000];

	readconstant(fyl, KEY, Value);
	return atoi(Value);
}

double readdouble(FILE *fyl, const char *KEY) {

	char Value[1000];

	readconstant(fyl, KEY, Value);
	return strtod(Value, NULL);
}

unsigned long readlong(FILE *fyl, const char *KEY) {

	char Value[1000];

	readconstant(fyl, KEY, Value);
	return strtoul(Value, NULL, 10);
}

/***********************************************************************************************/
