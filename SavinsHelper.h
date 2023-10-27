#ifndef SAVINS_HELPER
#define SAVINS_HELPER

double PreciseOverlapWEGA(int sizeFixed, double* fixed, double* weightFixed, int sizeVariable, double* variable, double* weightVariable);

double PreciseCalcOverlapVolAtomsWEGA(double atom1X, double atom1Y, double atom1Z, double atom2X, double atom2Y, double atom2Z);

void denormalize(double* x, double* bounds, double* x_vals, int sizeX);

void Rotate(double* matrix, int sizeMatrix, double theta, double axisAx, double axisAy, double axisAz, double axisBx, double axisBy, double axisBz, double x2, double y2, int sizeX);

void Move(double* matrix, int sizeMatrix, double deltaX, double deltaY, double deltaZ);

#endif
