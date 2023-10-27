#include "mex.h"
#include "SavinsHelper.h"

//ObjFunc(x, molMat1, weights1, ov1, molMat2, weights2, ov2, bounds)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double* x = mxGetPr(prhs[0]);
	int sizeX = mxGetNumberOfElements(prhs[0]);

	int sizeM1 = mxGetN(prhs[1]);//WARNING: EXPECTING A MOLECULE BY COLUMN
	double* molMat1 = mxGetPr(prhs[1]);
	double* weights1 = mxGetPr(prhs[2]);
	double ov1 = mxGetScalar(prhs[3]);
	int sizeM2 = mxGetN(prhs[4]);
	double* molMat2 = mxGetPr(prhs[4]);
	double* weights2 = mxGetPr(prhs[5]);
	double ov2 = mxGetScalar(prhs[6]);
	double* bounds = mxGetPr(prhs[7]);

	double x_vals[10]; //Let us alloc the longest on the stack to avoid malloc-free //double* x_vals = mxMalloc(sizeof(double)*sizeX);//Either 6 or 10 variables expected
	denormalize(x, bounds, x_vals, sizeX);
	double angle=0.0, p1x=0.0, p1y=0.0, p1z=0.0, p2x=0.0, p2y=0.0, p2z=0.0, x2=0.0, y2=0.0, deltax=0.0, deltay=0.0, deltaz=0.0;
	angle = x_vals[0];
	int shift = 6;
	if (sizeX==10){
		p1x = x_vals[1];
		p1y = x_vals[2];
		p1z = x_vals[3];
		p2x = x_vals[4];
		p2y = x_vals[5];
		p2z = x_vals[6];
	}else{
		x2 = x_vals[1];
		y2 = x_vals[2];
		shift = 2;
	}
	deltax = x_vals[shift + 1];
	deltay = x_vals[shift + 2];
	deltaz = x_vals[shift + 3];
	//mxFree(x_vals); //mexPrintf("angle = %lf\np1x = %lf\np1y = %lf\np1z = %lf\np2x = %lf\np2y = %lf\np2z = %lf\ndeltax = %lf\ndeltay = %lf\ndeltaz = %lf", angle, p1x, p1y, p1z, p2x, p2y, p2z, deltax, deltay, deltaz);
	
	double fval;
	if ( sizeX==10 && ( (p1x-p2x)==0.0 && (p1y-p2y)==0.0 && (p1z-p2z)==0.0) ){
		fval = 0.0;
	}else{
		double* modMolMat2 = mxMalloc(sizeof(double)*(3*sizeM2));
		for(int i = 0; i<(3*sizeM2); i++){
			modMolMat2[i] = molMat2[i];//Let us make a deep copy of the matrix to be modified!
		}//Let's make the computation:
		Rotate(modMolMat2, sizeM2, angle, p1x, p1y, p1z, p2x, p2y, p2z, x2, y2, sizeX);// Rotate(molMat2, angle, p1x, p1y, p1z, p2x, p2y, p2z); % Rotate (Remember: This modifies modMolMat2!)
		Move(modMolMat2, sizeM2, deltax, deltay, deltaz);// Move(modMolMat2, deltax, deltay, deltaz);
		double ov3 = PreciseOverlapWEGA(sizeM1, molMat1, weights1, sizeM2, modMolMat2, weights2);//mexPrintf("ov3 = %lf\n", ov3);
		mxFree(modMolMat2);
		fval = ( ov3 / (ov1 + ov2 - ov3) );//Tanimoto		
	}
	plhs[0] = mxCreateDoubleScalar(fval);//Return the value of the objective function
}

//function value = ObjFunc(x, molMat1, weights1, ov1, molMat2, weights2, ov2, bounds)
//    % X is the vector with 10 variables handled in [0, 1] by the optimizer
//    % bounds is a matrix 10x2, lower bound in the first column, upper in the second
    
//    % Denormalize the variables for the real process:
//    angle = x(1)*(bounds(1, 2) - bounds(1, 1)) + bounds(1, 1); % [0, 1] * (max - min) + min
//    p1x = x(2)*(bounds(2, 2) - bounds(2, 1)) + bounds(2, 1);
//    p1y = x(3)*(bounds(3, 2) - bounds(3, 1)) + bounds(3, 1);
//    p1z = x(4)*(bounds(4, 2) - bounds(4, 1)) + bounds(4, 1);
//    p2x = x(5)*(bounds(5, 2) - bounds(5, 1)) + bounds(5, 1);
//    p2y = x(6)*(bounds(6, 2) - bounds(6, 1)) + bounds(6, 1);
//    p2z = x(7)*(bounds(7, 2) - bounds(7, 1)) + bounds(7, 1);
//    deltax = x(8)*(bounds(8, 2) - bounds(8, 1)) + bounds(8, 1);
//    deltay = x(9)*(bounds(9, 2) - bounds(9, 1)) + bounds(9, 1);
//    deltaz = x(10)*(bounds(10, 2) - bounds(10, 1)) + bounds(10, 1);
    
//    if isequal(p1x, p2x) && isequal(p1y, p2y) && isequal(p1z, p2z)
//        value = 0.0;
//    else
//       % Let's make the computation:
//        modMolMat2 = Rotate(molMat2, angle, p1x, p1y, p1z, p2x, p2y, p2z); % Rotate
//        modMolMat2 = Move(modMolMat2, deltax, deltay, deltaz); % Move
    
//        ov3 = PreciseOverlapWEGA(molMat1, weights1, modMolMat2, weights2);
//        value = Tanimoto(ov1, ov2, ov3); 
//    end
//end
