#include "mex.h"
#include "SavinsHelper.h"

//CalculateWeights(atoms)

//nlhs = Number of output variables
//plhs = Array of mxArray pointers to the output variables
//nrhs = Number of input variables
//prhs = Array of mxArray pointers to the input variables

//See: http://walkingrandomly.com/?p=1795

double MyCalculateWeightWEGA(double vi, int numAtoms, double* atoms, int iAtom){
	double vij = 0.0;
	double atom1X = atoms[iAtom*3 + 0];//0*numAtoms + iAtom (for column-wise)
	double atom1Y = atoms[iAtom*3 + 1];
	double atom1Z = atoms[iAtom*3 + 2];

	double atom2X, atom2Y, atom2Z;
	for(int i = 0; i<numAtoms; i++){
		if(i==iAtom){
			continue;
		}
		atom2X = atoms[i*3 + 0];
		atom2Y = atoms[i*3 + 1];
		atom2Z = atoms[i*3 + 2];
		vij = vij + PreciseCalcOverlapVolAtomsWEGA(atom1X, atom1Y, atom1Z, atom2X, atom2Y, atom2Z);
	}
	return vi / (vi + 0.8665 * vij);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int numAtoms = mxGetN(prhs[0]);//WARNING: EXPECTING A MOLECULE BY COLUMN
	double* atoms = mxGetPr(prhs[0]);//Input atoms

	plhs[0] = mxCreateDoubleMatrix(1, numAtoms, mxREAL);
	double* weights = mxGetPr(plhs[0]);

	double vi = 4.0 * 3.14159265358 * (1.8*1.8*1.8) / 3.0;

	for(int i = 0; i<numAtoms; i++){
		weights[i] = MyCalculateWeightWEGA(vi, numAtoms, atoms, i); // ML: weights(i) = CalculateWeightWEGA(atoms(i, 1), atoms(i, 2), atoms(i, 3), atoms, i);
	}	
	//mexPrintf("Atomos: %d\n", numAtoms); // En memoria esta por columnas :-/
	//for(int i = 0; i<numAtoms; i++){
	//	for(int j = 0; j<3; j++){
	//		mexPrintf("%.2lf ", atoms[j*numAtoms + i]);
	//	}
	//	mexPrintf("\n");
	//}
}
