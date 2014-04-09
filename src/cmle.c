#include "allhead.h"  
//Functions for CMLE
void CWbMle(double *Dat, double Cx, int n, int m, double conCr, int nIter, double *RES)
{
	/*The function to calculate censored Weibull MLE with T1 or T2*/
	/*Initialize the local variable for RES*/
	double conInd =1.0;
	double alf =1.0;
	double eta =1.0; 

	if (m ==0)
	{
		goto label;
	}
	else
	{
		/*First copy the data we need*/
		double lucd[m], pucd[m];
		for (int i=0; i<m; i++)
		{
			lucd[i] = pucd[i] = Dat[i];
		}
		cvecLog(lucd, m);
		
		//Variable do not change over the iterations
		double id1 = cmean(lucd, m);
		
		//Variables changes with time
		double np1, np2, d12, d1, pucdsum;
		
		//Temporary variables
		double oal, ial;
		double flag = 100.0;
		int icount = 0;
		
		/* the real calculate begins here*/
		alf = cmean(Dat, m)/csd(Dat, m);
		while(flag > conCr && icount< nIter)
		{
			oal = alf;
			cvecPow(Dat, alf, m, pucd);
			np1 = cvecTimes_Sum(pucd, lucd, m);
			pucdsum = csum(pucd, m);
			d12 = (double)(n-m) * pow(Cx, alf);
			np2 = d12 * log(Cx);
			
			d1 =  pucdsum + d12;
			ial = (np1 + np2)/d1 - id1;
			
			if (ial < conCr || isnan(ial))
			{
				conInd =2.0;
				goto label;
			}
			else
			{
				alf = 1.0/ial;
				flag = dabs((ial - 1.0/oal));
				//flag = dabs((alf - oal));
			}
			icount++;
		}
		if (icount< nIter)
		{	
			conInd = 0;
			eta=pow(d1/(double) m, ial);
		}
	}
	label: RES[0] = conInd;
	RES[1] = alf;
	RES[2] = eta;
}


SEXP R2C_cenWeibullMLE(SEXP RDAT, SEXP RC, SEXP RN, SEXP RM, SEXP RConCr, SEXP RNIter)
{
	/*The interface between R and C on the Censored WeibullMLE*/
	PROTECT(RConCr = coerceVector(RConCr, REALSXP));
	double conCr = *REAL(RConCr);
	PROTECT(RNIter = coerceVector(RNIter, INTSXP));
	int nIter =*INTEGER(coerceVector(RNIter,INTSXP));
	
	/*Get the data ready*/
	PROTECT(RN = coerceVector(RN, INTSXP));
	int n =*INTEGER(coerceVector(RN,INTSXP));
	PROTECT(RM = coerceVector(RM, INTSXP));
	int m =*INTEGER(coerceVector(RM,INTSXP));
	
	PROTECT(RDAT = coerceVector(RDAT, REALSXP));
	PROTECT(RC = coerceVector(RC, REALSXP));
	
	double *ODAT;
	double Dat[m], Cx;
	ODAT = REAL(RDAT); //pointer to the orginal data;
	Cx = *REAL(RC);
	for (int i=0; i<m; i++)
		Dat[i] = ODAT[i];
	
	/* Get the output ready*/
	SEXP RRES;
	PROTECT(RRES = allocVector(REALSXP, 3));
	double *RES = REAL(RRES);
	/*Call the function to do it*/
	CWbMle(Dat, Cx, n, m, conCr, nIter, RES);
	
	UNPROTECT(7);
	return(RRES);
}
