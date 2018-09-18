#include "mex.h"
/*============================================================================
 *
 *  Name        : M40_MkSum
 *
 *  Description : 
 *  calculate sum market by market
 *
 *===========================================================================*/
void M40_Sum(int nM,               /* number of markets         */
            int nCol,                 /* number of columns to sum over  */
            int nObs,                 /* total number of products  */
			double MidxL[],          /* number of products in each market */
            double DtaM[],                /* data matrix    */
			double tsum[])              /* sum for each market  */
{
	double sum_dta;
	int i, j, k, ktmp1, ktmp2, tn, id1;

	for (i=0; i<nM; i++)
	{
		id1 = MidxL[i];
		tn  = MidxL[i+1] - id1;

        for (k=0; k<nCol; k++)
        {
            ktmp1=k*nObs;
            ktmp2=k*nM;

            sum_dta=0;
            for (j=0; j<tn; j++)
            {
                sum_dta += DtaM[j+id1+ktmp1];
            }
            tsum[i+ktmp2] = sum_dta;
        }
	}
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int nM, nCol, nObs;
    double *MidxL, *DtaM, *tsum;
    
    if (nrhs!=5)
    {
        mexErrMsgTxt("Five input required.");
    }
    else if (nlhs > 1)
    {
        mexErrMsgTxt("Too many output arguments. One output only.");
    }
    
    /* assign the points to the input matrices  */
	MidxL = (mxGetPr(prhs[3]));
    DtaM = mxGetPr(prhs[4]);
    
    /* assign the int values of input   */
    nM = (int)mxGetScalar(prhs[0]);
    nCol = (int)mxGetScalar(prhs[1]);
    nObs = (int)mxGetScalar(prhs[2]);
    
    /* allocate memory for output: tsum     */
    plhs[0] = mxCreateDoubleMatrix(nM, nCol, mxREAL);
    tsum = mxGetPr(plhs[0]);
    
    /* call the C function      */
    M40_Sum(nM, nCol, nObs, MidxL, DtaM, tsum);
}
