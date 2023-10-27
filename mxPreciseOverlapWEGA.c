#include "mex.h"
#include "SavinsHelper.h"

//PreciseOverlapWEGA(fixed, weightFixed, variable, weightVariable)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int sizeFixed = mxGetN(prhs[0]);//WARNING: EXPECTING A MOLECULE BY COLUMN	
	int sizeVariable = mxGetN(prhs[2]);//WARNING: EXPECTING A MOLECULE BY COLUMN
	
	double* fixed = mxGetPr(prhs[0]);
	double* weightFixed = mxGetPr(prhs[1]);
	double* variable = mxGetPr(prhs[2]);
	double* weightVariable = mxGetPr(prhs[3]);

	double overlap = PreciseOverlapWEGA(sizeFixed, fixed, weightFixed, sizeVariable, variable, weightVariable);
	plhs[0] = mxCreateDoubleScalar(overlap);//Return the overlap
}
