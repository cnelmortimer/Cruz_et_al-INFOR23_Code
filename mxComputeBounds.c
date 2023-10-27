#include "mex.h"
#include <math.h>

//ComputeBounds(molMat1, molMat2)

int min(double a, double b){
	return (a<b)?1:0;
}

int max(double a, double b){
	return (a>b)?1:0;
}

double find(int (*func)(double, double), int numAtoms, double* atoms, int focus){//focus 0 = X, focus 1 = Y, focus 2 = Z	
	double ref = atoms[focus];
	double candidate;
	for(int i = focus + 3; i<(3*numAtoms); i+=3){
		candidate = atoms[i];//mexPrintf("Ref: %.2lf; Candidate: %.2lf\n", ref, candidate);
		if(func(candidate, ref)){
			ref = candidate;
		}
	}
	return ref;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if(nrhs != 3){
		mexErrMsgTxt("Wrong parameters for mxComputeBounds. Provide: molMat1, molMat2, dimSize.");
	}

	int sizeM1 = mxGetN(prhs[0]);//WARNING: EXPECTING A MOLECULE BY COLUMN
	double* molMat1 = mxGetPr(prhs[0]);
	int sizeM2 = mxGetN(prhs[1]);
	double* molMat2 = mxGetPr(prhs[1]);
	int sizeX = (int) mxGetScalar(prhs[2]);
	
	plhs[0] = mxCreateDoubleMatrix(2, sizeX, mxREAL);//Every column contains the lower and upper bound
	double* bounds = mxGetPr(plhs[0]);
	bounds[0] = 0.0; bounds[1] = 2.0*3.141592653589793;

	double bn2 = find(min, sizeM2, molMat2, 0); double bn3 = find(max, sizeM2, molMat2, 0);
	double bn4 = find(min, sizeM2, molMat2, 1); double bn5 = find(max, sizeM2, molMat2, 1);
	double bn6 = find(min, sizeM2, molMat2, 2); double bn7 = find(max, sizeM2, molMat2, 2);

	double tmpX [] = { fabs( bn2 - find(min, sizeM1, molMat1, 0) ), fabs( bn3 - find(max, sizeM1, molMat1, 0) )};
	double tmpY [] = { fabs( bn4 - find(min, sizeM1, molMat1, 1) ), fabs( bn5 - find(max, sizeM1, molMat1, 1) )};
	double tmpZ [] = { fabs( bn6 - find(min, sizeM1, molMat1, 2) ), fabs( bn7 - find(max, sizeM1, molMat1, 2) )};

	double diffX = (tmpX[0]>tmpX[1])?tmpX[0]:tmpX[1];
	double diffY = (tmpY[0]>tmpY[1])?tmpY[0]:tmpY[1];
	double diffZ = (tmpZ[0]>tmpZ[1])?tmpZ[0]:tmpZ[1];

	if(sizeX==10){	
		bounds[2] = bn2; bounds[3] = bn3;
		bounds[4] = bn4; bounds[5] = bn5;
		bounds[6] = bn6; bounds[7] = bn7;
		bounds[8] = bn2; bounds[9] = bn3;
		bounds[10] = bn4; bounds[11] = bn5;
		bounds[12] = bn6; bounds[13] = bn7;

		bounds[14] = -diffX; bounds[15] = diffX;
		bounds[16] = -diffY; bounds[17] = diffY;
		bounds[18] = -diffZ; bounds[19] = diffZ;
	}else{
		bounds[2] = 0.0; bounds[3] = 2.0*3.141592653589793;
		bounds[4] = 0.0; bounds[5] = 3.141592653589793*0.5;

		bounds[6] = -diffX; bounds[7] = diffX;
		bounds[8] = -diffY; bounds[9] = diffY;
		bounds[10] = -diffZ; bounds[11] = diffZ;
	}
}

//function bounds = ComputeBounds(molMat1, molMat2)
//    bounds = zeros(10, 2); % Each row is: [low_bound, up_bound]. 10 vars
//    bounds(1, :) = [0, 2*pi]; % angle
//    bounds(2, :) = [min(molMat2(:, 1)), max(molMat2(:, 1))];
//    bounds(3, :) = [min(molMat2(:, 2)), max(molMat2(:, 2))];
//    bounds(4, :) = [min(molMat2(:, 3)), max(molMat2(:, 3))];
//    bounds(5,: ) = bounds(2, :);
//    bounds(6,: ) = bounds(3, :);
//    bounds(7,: ) = bounds(4, :);
//    %------------
//    tmpX = [min(molMat1(:, 1)), max(molMat1(:, 1))];
//    tmpY = [min(molMat1(:, 2)), max(molMat1(:, 2))];
//    tmpZ = [min(molMat1(:, 3)), max(molMat1(:, 3))];
//    
//    diffX = max(abs(bounds(2,:) - tmpX));
//    diffY = max(abs(bounds(3,:) - tmpY));
//    diffZ = max(abs(bounds(4,:) - tmpZ));
//    
//    bounds(8, :) = [-diffX, diffX];
//    bounds(9, :) = [-diffY, diffY];
//    bounds(10, :) = [-diffZ, diffZ];
//end
