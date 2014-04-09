/*Weibull Mixture, Uncensored mixture */
#include "allhead.h"  
double cdfWbMix(double x, double *parm, double qInt)
{
	double cdf = parm[0] * pweibull(x, parm[2], parm[4], 1, 0);
	cdf += parm[1] * pweibull(x, parm[3], parm[5], 1, 0);
	return(cdf - qInt);
}

double quanWbMix(double *parm, double qInt)
{
	/*Specified to 0.05 */
	double a=0.0, b=15.0;
	double tol = 1E-10;
	double c=(a+b)/2.0;
	double func = cdfWbMix(c, parm, qInt);
	double step = b-a;
	int icount=0, maxIte =1000;
	
	while (icount<maxIte)
	{
		icount++;
		if (func==0 || step< tol)
			break;
		if (func <0)
			a=c;
		else
			b=c;
		c = (a+b)/2.0;
		step = b-a;
		func = cdfWbMix(c, parm, qInt);
	}
	return(c);
}


double nllh_CWbMix(double *dat, double Cx, int n, int m, double *parm)
{
	//Th negativel Log-likelihood of the censored Weibull mixture
	double nllh=0;
	double pdf1, pdf2, sf1, sf2;
	for (int i=0; i<m; i++)
	{	
		pdf1 = parm[0] * dweibull(dat[i], parm[2], parm[4], 0);
		pdf2 = parm[1] * dweibull(dat[i], parm[3], parm[5], 0);
		nllh+= log(pdf1 + pdf2);
	}
	if (n>m)
	{
		sf1 = parm[0] * pweibull(Cx, parm[2], parm[4], 0, 0);
		sf2 = parm[1] * pweibull(Cx, parm[3], parm[5], 0, 0);
		nllh += log(sf1+ sf2) * (double) (n-m);
	}
	return(-nllh);
}

void eStep(double *dat, double Cx, int n, int m, double *parm, double *eVec1, double *eVec2)
{
	/*The Expectation Step of EM Algorithm*/
	double sf1 = parm[0] * pweibull(Cx, parm[2], parm[4], 0, 0);
	double sf2 = parm[1] * pweibull(Cx, parm[3], parm[5], 0, 0);
	eVec1[m] = sf1/(sf1 + sf2);
	eVec2[m] = 1.0 - eVec1[m];
	
	double pdf1, pdf2;
	int i=0;
	for(i=0; i<m; i++)
	{
		pdf1 = parm[0] * dweibull(dat[i], parm[2], parm[4], 0);
		pdf2 = parm[1] * dweibull(dat[i], parm[3], parm[5], 0);
		eVec1[i] = pdf1/(pdf1+ pdf2);
		eVec2[i] = 1- eVec1[i];
	}
}

void mStep(double *Dat, double *eVec, double Cx, int n, int m, double *Cmle, double conCr, int nIter)
{
	/*The function to calculate censored Weibull MLE with T1 or T2*/
	/*Initialize the local variable for RES*/
	double conInd =1.0;
	double alf = Cmle[1];
	double eta =1.0; 

	/*First copy the data we need*/
	double lucd[m], pucd[m], p1_d1Vec[m]; ;// Level 1 temporary variable
	double p1_n1, p1_d1; //Sum of the above variables;
	double p1_n2, p1_d2; //Level 2, variables based on Cx
	
	// Variables in alpha, but does involve alpha.
	double	p2_d, p2; 
	for (int i=0; i<m; i++)
	{
		lucd[i] = Dat[i];
	}
	cvecLog(lucd, m);
	p2_d = csum(eVec, m);
	p2 = cvecTimes_Sum(eVec, lucd, m)/p2_d;
	
	/*Temporary Variables controls on the iteration*/
	double oal, ial;
	double flag = 100.0;
	int icount = 0;
	
	/* the real calculate begins here*/
	while(flag > conCr && icount< nIter)
	{
		oal = alf;
		cvecPow(Dat, alf, m, pucd);
		cvecTimes(pucd, eVec, m, p1_d1Vec);
		
		p1_d1 = csum(p1_d1Vec, m);
		p1_n1 = cvecTimes_Sum(p1_d1Vec, lucd, m);
		
		p1_d2 = (double) (n-m) * pow(Cx, alf) * eVec[m];
		p1_n2 = p1_d2 * log(Cx);
		
		ial = (p1_n1 + p1_n2)/(p1_d1 + p1_d2) - p2;

		if (ial < 1E-10 || isnan(ial))
		{
			conInd=2.0;
			goto label;
		}
		else
		{
			alf = 1.0/ial;
			flag = dabs((ial - 1.0/oal));
		}
		icount++;
	}
	if (icount< nIter)
	{	
		conInd = 0;
		eta=pow((p1_d1 + p1_d2)/p2_d, ial);
	}
	label: Cmle[0] = conInd;
	Cmle[1] = alf;
	Cmle[2] = eta;
}

void emWbMix(double *dat, double Cx, int n, int m, double *parm, double conCr, int nIter)
{
	//The EM algorithm for (un) censored Wb mixture 
	//Implemented in a Type II setting
	//The parm is a vector of length 7, with the 6th element as the convergence flag;
	//7th as the negative log-likelihood
	double mConCr = conCr*1E-4;
	int mNIter = nIter/10;
	//Allocate the expectation vectors
	double eVec1[(m+1)], eVec2[(m+1)];
	//Allocate the temporary variables for maximization 
	double mleTemp1[3], mleTemp2[3];
	mleTemp1[0] = mleTemp2[0] =0.0;
	mleTemp1[1] = parm[2];
	mleTemp1[2] = parm[4];
	mleTemp2[1] = parm[3];
	mleTemp2[2] = parm[5];
	//The nllh and control of the overal EM algorithm
	double tnllh = dabs(nllh_CWbMix(dat, Cx, n, m, parm));
	double onllh;
	int icount = 0;
	double nllhDif= 100.0;
	double nmDif =(double) (n-m);
	
	//The EM algorithm covergence flage, 0 converge, 1, wrong nllh, 2, fail of Maximization
	double emflag =0.0;
	
	while(icount <nIter && nllhDif > conCr)
	{
		icount++;
		onllh = tnllh;
		//The E-step
		eStep(dat, Cx, n, m, parm, eVec1, eVec2);
		parm[0] = (csum(eVec1, m) + nmDif * eVec1[m])/n;
		parm[1] = 1.0- parm[0];
		//The M-step
		mStep(dat, eVec1, Cx, n, m, mleTemp1, mConCr, mNIter);
		mStep(dat, eVec2, Cx, n, m, mleTemp2, mConCr, mNIter);
		if (mleTemp1[0]==0 && mleTemp2[0]==0)
		{
			parm[2] = mleTemp1[1];
			parm[4] = mleTemp1[2];
			parm[3] = mleTemp2[1];
			parm[5] = mleTemp2[2];
			tnllh = nllh_CWbMix(dat, Cx, n, m, parm);
			parm[7] = tnllh;
			if (dabs(tnllh) > 1E100)
			{	
				emflag =1.0;
				goto label;
			}
			else
			{
				nllhDif = dabs(tnllh - onllh)/onllh;
			}
		}
		else
		{
			emflag=2.0;
			goto label;
		}
	}
	if (icount== nIter)	emflag=3.0;
	label: parm[6] =emflag;
}

SEXP R2C_CWbMix(SEXP RDAT, SEXP RC, SEXP RN, SEXP RM, SEXP RIni, SEXP RConCr, SEXP RNIter)
{
	PROTECT(RConCr = coerceVector(RConCr, REALSXP));
	double conCr = *REAL(RConCr);
	PROTECT(RNIter = coerceVector(RNIter, INTSXP));
	int nIter =*INTEGER(coerceVector(RNIter,INTSXP));
	
	PROTECT(RC = coerceVector(RC, REALSXP));
	double Cx;
	Cx = *REAL(RC);
	
	PROTECT(RN = coerceVector(RN, INTSXP));
	int n =*INTEGER(coerceVector(RN,INTSXP));
	PROTECT(RM = coerceVector(RM, INTSXP));
	int m =*INTEGER(coerceVector(RM,INTSXP));
	
	PROTECT(RDAT = coerceVector(RDAT, REALSXP));
	double *ODAT;
	double Dat[m];
	ODAT = REAL(RDAT); //pointer to the orginal data;
	for (int i=0; i<m; i++)
		Dat[i] = ODAT[i];
	
	PROTECT(RIni = coerceVector(RIni, REALSXP));
	double *OIni;
	double parm[8];
	OIni = REAL(RIni);
	for (int i=0; i<6; i++)
		parm[i] = OIni[i];
	parm[6] = 0.0;
	parm[7] = 0.0;
	
	emWbMix(Dat, Cx, n, m, parm, conCr, nIter);
	
	SEXP RRES;
	PROTECT(RRES = allocVector(REALSXP, 8));
	double *RES = REAL(RRES);
	RES[0]= parm[6];
	RES[1]= parm[7];
	RES[2]= parm[0];
	RES[3]= parm[1];
	RES[4]= parm[2];
	RES[5]= parm[3];
	RES[6]= parm[4];
	RES[7]= parm[5];
	UNPROTECT(8);
	return RRES;
}
