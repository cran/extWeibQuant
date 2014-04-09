//The basis vector function
#include "allhead.h"  
int ccenIndex(double *org, double tc, int tlen)
{
	/*Get the index of censoring*/
	int i=0, alen;
	double ttc = tc+1e-10;
	while (org[i]<= ttc && i< tlen) i++;
	alen = (i<3) ? 0:i;
	return alen;
}

double csum(double *vec, int tlen)
{
	int i=0;
	double tsum =0;
	for (i=0; i<tlen; i++)
		tsum +=vec[i];
	return(tsum);
}
double cmean(double *vec, int tlen)
{
	double tsum, tmean;
	tsum= csum(vec, tlen);
	tmean= tsum/(double) tlen;
	return(tmean);
}

double cvar(double *vec, int tlen)
{
	double ssum =0;
	int i;
	for (i=0; i<tlen; i++)
		ssum += vec[i]*vec[i];
	double omean = cmean(vec, tlen);
	double dtlen = (double) tlen; // coerce the integer into double;
	double tsd = (ssum - dtlen*omean*omean)/(dtlen -1.0);
	return(tsd);
}

double cvar_adj(double *vec, int tlen)
{
	double ssum =0;
	int i;
	for (i=0; i<tlen; i++)
		ssum += vec[i]*vec[i];
	double omean = cmean(vec, tlen);
	double dtlen = (double) tlen; // coerce the integer into double;
	double tsd = (ssum - dtlen*omean*omean)/dtlen;
	return(tsd);
}

double csd(double *vec, int tlen)
{
	double ssum =0;
	int i;
	for (i=0; i<tlen; i++)
		ssum += vec[i]*vec[i];
	double omean = cmean(vec, tlen);
	double dtlen = (double) tlen; // coerce the integer into double;
	double tsd = (ssum - dtlen*omean*omean)/(dtlen -1.0);
	return(sqrt(tsd));
}
void cvecPow(double *vec, double y, int tlen, double *nvec)
{
	/* the function to perform power over the vec*/
	for (int i=0; i<tlen; i++)
		nvec[i] = pow(vec[i], y);
}

void cvecTimes(double *vec1, double *vec2, int tlen, double *tvec)
{
	for (int i=0; i<tlen; i++)
		tvec[i] = vec1[i]* vec2[i];
}

double cvecTimes_Sum(double *vec1, double *vec2, int tlen)
{
	double tsum =0;
	for (int i=0; i<tlen; i++)
		tsum += vec1[i]* vec2[i];
	return(tsum);
}

void cvecLog(double *vec, int tlen)
{
	for (int i=0; i<tlen; i++)
		vec[i] = log(vec[i]);
}

double dabs(double aa)
{
	double bb;
	bb= aa>0.0 ? aa:(-aa);
	return(bb);
}
