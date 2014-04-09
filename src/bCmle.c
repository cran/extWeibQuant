/*The bootstrap CMLE*/
#include "allhead.h" 
void cbootstrap(double *odt, int n, int bn, double *ndt)
{	/*The bootstrap function with sort*/
    int i=bn;
	double *tp_ndt = ndt;
    do
	{
		*tp_ndt++ = odt[(int)(n * unif_rand())];
	} while(--i);
	R_qsort(ndt, 1, bn);
}

void BstpCWbMle(double *dat, double qInt, int *indSet,
				int n, int B, int K, double conCr, int nIter, 
				double *qnTilVec)
{	
	double rdt[n]; // the array for the bootstraped data
	int i=B;
	int resInd=0;
	double pest[3]={0.0};
	double tCx;
	GetRNGstate();
	do
	{
		cbootstrap(dat, n, n, rdt);
		for (int k=0; k<K; k++)
		{
			//Type II
			tCx = rdt[(indSet[k]-1)];
			CWbMle(rdt, tCx, n, indSet[k], conCr, nIter, pest);
			if (pest[0]==0)
			{
				qnTilVec[resInd] = qweibull(qInt, pest[1], pest[2], 1, 0);
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

SEXP R2C_boostrapCMLE(SEXP RDAT, SEXP RQInt, SEXP RIndSet, SEXP RB, SEXP RConCr, SEXP RNIter)
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
	
	BstpCWbMle(Dat, qInt, indSet, n, B, K, conCr, nIter, qnTilVec);
	
	for (int i=0; i<M; i++)
		RES[i] = qnTilVec[i];
	UNPROTECT(7);
	return RRES;
}
