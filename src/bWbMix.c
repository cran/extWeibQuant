#include "allhead.h" 
void BstpCWbMix(double *dat, double qInt, int *indSet, double *iniVec,
				int n, int B, int K, double conCr, int nIter, 
				double *qnTilVec)
{	
	double rdt[n]; // the array for the bootstraped data
	int i=B;
	int resInd=0;
	double mixvec[8]={0.0};
	double tCx;
	GetRNGstate();
	do
	{
		cbootstrap(dat, n, n, rdt);
		for (int k=0; k<K; k++)
		{
			//Type II
			tCx = rdt[(indSet[k]-1)];
			for (int j=0; j<6; j++)
				mixvec[j] = iniVec[(j+k*6)];
			mixvec[6]=0.0;
			mixvec[7]=0.0;
			emWbMix(rdt, tCx, n, indSet[k], mixvec, conCr, nIter);
			if (mixvec[6]==0)
			{
				qnTilVec[resInd] = quanWbMix(mixvec, qInt);
			}
			else
			{
				qnTilVec[resInd] =-1.0;
			}
			resInd++;
		}
	} while(--i);
	PutRNGstate();
}

SEXP R2C_bstpWbMix(SEXP RDAT, SEXP RQInt, SEXP RIndSet, SEXP RIniVec, SEXP RB, SEXP RConCr, SEXP RNIter)
{
	int n = length(RDAT);
	PROTECT(RDAT = coerceVector(RDAT, REALSXP));
	double *ODAT;
	double Dat[n];
	ODAT = REAL(RDAT); //pointer to the orginal data;
	for (int i=0; i<n; i++)
		Dat[i] = ODAT[i];
		
	PROTECT(RQInt = coerceVector(RQInt, REALSXP));
	double qInt = *REAL(RQInt);
	
	int K= length(RIndSet);
	PROTECT(RIndSet = coerceVector(RIndSet, INTSXP));
	int *OInd;
	int indSet[K];
	OInd = INTEGER(RIndSet); //pointer to the orginal data;
	for (int k=0; k<K; k++)
		indSet[k] = OInd[k];	
	
	PROTECT(RIniVec = coerceVector(RIniVec, REALSXP));
	double *OIniVec;
	OIniVec = REAL(RIniVec);
	double iniVec[(K*6)];
	for (int k=0; k<(K*6); k++)
		iniVec[k]= OIniVec[k];
	
	PROTECT(RB = coerceVector(RB, INTSXP));
	int B =*INTEGER(coerceVector(RB,INTSXP));
	
	PROTECT(RConCr = coerceVector(RConCr, REALSXP));
	double conCr = *REAL(RConCr);
	PROTECT(RNIter = coerceVector(RNIter, INTSXP));
	int nIter =*INTEGER(coerceVector(RNIter,INTSXP));
		
	int M = (K*B);
	double qnTilVec[M];
	SEXP RRES;
	PROTECT(RRES = allocVector(REALSXP, M));
	double *RES = REAL(RRES);
	
	BstpCWbMix(Dat, qInt, indSet, iniVec, n, B, K, conCr, nIter, qnTilVec);
	
	for (int i=0; i<M; i++)
		RES[i] = qnTilVec[i];
	UNPROTECT(8);
	return RRES;
}
