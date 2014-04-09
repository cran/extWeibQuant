//Global Head
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>  
/* Rule of the names of variables:
1. All CAPIalf: RDAT RC, the one passed from/to R; RES only updated at the last
2. Capialf first: Dat, the original copy of the data. Do not change;
3. all small, local for computation;
*/
int ccenIndex(double *org, double tc, int tlen);
double csum(double *vec, int tlen);
double cmean(double *vec, int tlen);
double cvar(double *vec, int tlen);
double cvar_adj(double *vec, int tlen);
double csd(double *vec, int tlen);
void cvecPow(double *vec, double y, int tlen, double *nvec);
void cvecTimes(double *vec1, double *vec2, int tlen, double *tvec);
double cvecTimes_Sum(double *vec1, double *vec2, int tlen);
void cvecLog(double *vec, int tlen);
double dabs(double aa);
void CWbMle(double *Dat, double Cx, int n, int m, double conCr, int nIter, double *RES);
SEXP R2C_cenWeibullMLE(SEXP RDAT, SEXP RC, SEXP RN, SEXP RM, SEXP RConCr, SEXP RNIter);
double nllh_CWbMix(double *dat, double Cx, int n, int m, double *parm);
void eStep(double *dat, double Cx, int n, int m, double *parm, double *eVec1, double *eVec2);
void mStep(double *Dat, double *eVec, double Cx, int n, int m, double *Cmle, double conCr, int nIter);
void emWbMix(double *dat, double Cx, int n, int m, double *parm, double conCr, int nIter);
SEXP R2C_CWbMix(SEXP RDAT, SEXP RC, SEXP RN, SEXP RM, SEXP RIni, SEXP RConCr, SEXP RNIter);
double cdfWbMix(double x, double *parm, double qInt);
double quanWbMix(double *parm, double qInt);
void cbootstrap(double *odt, int n, int bn, double *ndt);
void BstpCWbMle(double *dat, double qInt, int *indSet, int n, int B, int K, double conCr, int nIter, double *qnTilVec);
SEXP R2C_boostrapCMLE(SEXP RDAT, SEXP RQInt, SEXP RIndSet, SEXP RB, SEXP RConCr, SEXP RNIter);
void BstpCWbMix(double *dat, double qInt, int *indSet, double *iniVec, 
				int n, int B, int K, double conCr, int nIter, 
				double *qnTilVec);
SEXP R2C_bstpWbMix(SEXP RDAT, SEXP RQInt, SEXP RIndSet, SEXP RIniVec, SEXP RB, SEXP RConCr, SEXP RNIter);
				





