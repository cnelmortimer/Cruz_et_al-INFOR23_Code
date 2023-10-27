#include "SavinsHelper.h"
#include <math.h>

double PreciseOverlapWEGA(int sizeFixed, double* fixed, double* weightFixed, int sizeVariable, double* variable, double* weightVariable){
	double overlap = 0.0;

	double weightFixedi, fixediX, fixediY, fixediZ;
	double variablejX, variablejY, variablejZ;
	int i3, j3;
	for(int i = 0; i<sizeFixed; i++){
		weightFixedi = weightFixed[i];
		i3 = i*3;
		fixediX = fixed[i3 + 0];
		fixediY = fixed[i3 + 1];
		fixediZ = fixed[i3 + 2];
		for(int j = 0; j<sizeVariable; j++){
			j3 = j*3;
			overlap = overlap + ( weightFixedi*weightVariable[j]*PreciseCalcOverlapVolAtomsWEGA(fixediX, fixediY, fixediZ, variable[j3 + 0], variable[j3 + 1], variable[j3 + 2]) );
		}
	}

	return overlap;
}

double PreciseCalcOverlapVolAtomsWEGA(double atom1X, double atom1Y, double atom1Z, double atom2X, double atom2Y, double atom2Z){
	return 24.428790199 * exp( ( (atom1X - atom2X)*(atom1X - atom2X) + (atom1Y - atom2Y)*(atom1Y - atom2Y) + (atom1Z - atom2Z)*(atom1Z - atom2Z) ) * -0.3731438999881213);
}

void denormalize(double* x, double* bounds, double* x_vals, int sizeX){
	double upperBound, lowerBound;//bounds is expected to store the lower and upper bounds in sequence
	for(int i = 0, iV = 0; i<sizeX; i++, iV+=2){
		lowerBound = bounds[iV];
		upperBound = bounds[iV + 1];
		x_vals[i] = x[i]*(upperBound - lowerBound) + lowerBound;//[0, 1]*(max-min) + min //mexPrintf("x[%d] = %lf; in [%lf, %lf] => %lf\n", i, x[i], lowerBound, upperBound, x_vals[i]);
	}
}

//Rotate(molMat2, angle, p1x, p1y, p1z, p2x, p2y, p2z);
void Rotate(double* matrix, int sizeMatrix, double theta, double axisAx, double axisAy, double axisAz, double axisBx, double axisBy, double axisBz, double x2, double y2, int sizeX){
	double vectorX, vectorY, vectorZ;
	if(sizeX==10){//We do not care about x2, y2	
		vectorX = axisBx - axisAx;
		vectorY = axisBy - axisAy;
		vectorZ = axisBz - axisAz;
		double modVector = sqrt(vectorX*vectorX + vectorY*vectorY + vectorZ*vectorZ);
		vectorX = vectorX / modVector;
		vectorY = vectorY / modVector;
		vectorZ = vectorZ / modVector;
	}else{
		vectorX = cos(x2)*sin(y2);
        	vectorY = sin(x2)*sin(y2);
        	vectorZ = cos(y2); // With this approach, it is already normalized!: https://www.wolframalpha.com/input?i=norm%28cos%28x%29*sin%28y%29%2C+sin%28x%29*sin%28y%29%2C+cos%28y%29%29&lang=es
	}

	double angle = theta*0.5;
	double sinTheta = sin(angle);
        double q1w = cos(angle);

	double q1x = sinTheta*vectorX;
	double q1y = sinTheta*vectorY;
	double q1z = sinTheta*vectorZ;

	double q1xConjugated = -1.0 * q1x;
	double q1yConjugated = -1.0 * q1y;
	double q1zConjugated = -1.0 * q1z;

	double atomPositionMinusAX, atomPositionMinusAY, atomPositionMinusAZ;
	double part2w, part2x, part2y, part2z;
	double part3x, part3y, part3z;

	int focus = 0;
	for(int i = 0; i<sizeMatrix; i++){
		focus = 3*i;
		atomPositionMinusAX = matrix[focus] - axisAx;
		atomPositionMinusAY = matrix[focus + 1] - axisAy;
		atomPositionMinusAZ = matrix[focus + 2] - axisAz;

		part2w = - q1x * atomPositionMinusAX - q1y * atomPositionMinusAY - q1z * atomPositionMinusAZ;
		part2x = q1w * atomPositionMinusAX + q1y * atomPositionMinusAZ - q1z * atomPositionMinusAY;
		part2y = q1w * atomPositionMinusAY - q1x * atomPositionMinusAZ + q1z * atomPositionMinusAX;
		part2z = q1w * atomPositionMinusAZ + q1x * atomPositionMinusAY - q1y * atomPositionMinusAX;
        
		part3x = part2w * q1xConjugated + part2x * q1w + part2y * q1zConjugated - part2z * q1yConjugated;
		part3y = part2w * q1yConjugated - part2x * q1zConjugated + part2y * q1w + part2z * q1xConjugated;
		part3z = part2w * q1zConjugated + part2x * q1yConjugated - part2y * q1xConjugated + part2z * q1w;

		matrix[focus] = axisAx + part3x;
		matrix[focus + 1] = axisAy + part3y;
		matrix[focus + 2] = axisAz + part3z;
	}
}

//Move(modMolMat2, deltax, deltay, deltaz);
void Move(double* matrix, int sizeMatrix, double deltaX, double deltaY, double deltaZ){
	int i3 = 0;	
	for(int i = 0; i<sizeMatrix; i++){
		i3 = i*3;
		matrix[i3] += deltaX;
		matrix[i3 + 1] += deltaY;
		matrix[i3 + 2] += deltaZ;
	}
}
