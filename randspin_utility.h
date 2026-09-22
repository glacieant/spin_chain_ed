//
//  randspin_utility.h
//  

#ifndef LIBRARY_INCLUSION
#define LIBRARY_INCLUSION

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


#endif



#ifndef PROTOTYPE_DEFN
#define PROTOTYPE_DEFN

int fonecount(int x, int sitenum);
double pointspin(int x, int y);
int spinflip(int sitenum,int x, int j, int k, int *LinLeast,int *LinMax);
int comb(int x, int y);
void diagonaliseA(double *SpinHamiltonian, int dim, double *EigenValues, double *EigenVectors);
double spinsector(int sitenum, int statenum, int eigenum, double *EigenVectors, int *LinLeast, int *LinMax , int *StateKets);
int localspindiff(int loweig, int higheig, int sitenum, int statenum, double *EigenVectors, int *StateKets, int *LinLeast, int *LinMax , double *measure, double *merror);
int localspindiff_singlet(int loweig, int higheig, int SBEGIN, int SINGNUM, int sitenum, int statenum, double *EigenVectors, int *StateKets, int *LinLeast, int *LinMax, double *measure, double *merror);
double avgspingap(int loweig, int higheig, int sitenum, int statenum, double *EigenVectors, int *StateKets);
void *allocate(size_t num, size_t size);
FILE *openfile(const char *name, const char *mode);
void readconstant(FILE *fyl, const char *KEY, char *Value);
int readint(FILE *fyl, const char *KEY);
double readdouble(FILE *fyl, const char *KEY);
unsigned long readlong(FILE *fyl, const char *KEY);

extern void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
	double *w, double *Work, int *lwork, int *info);

void uniformdoublerand(double *RandomNumbers, int nrand, double llim, double ulim, unsigned long SEED);
void powerdoublerand(double *RandomNumbers, int nrand, double DELTA, unsigned long SEED);

void randomham(double *ExchangeCoupling,double *ZeemanField,int sitenum,int statenum, int *LinLeast, int *LinMax, int *StateKets, double *EigenValues,double *EigenVectors);
void randstate(int sitenum,int SZ,int *OneCount,int *LinLeast, int *LinMax, int *StateKets);

#endif
