#pragma once
#include "auxfunc.h"
#include "model.h"
#define DARK_BEAM_PL 1000000
struct QuickCalcMem{
	double CosBetaPi180;
	double SinBetaPi180;
	double CosAlphaPi180;
	double SinAlphaPi180;
	double Cos2AlphaPi180;
	double CosZ_ZsPi180;
	double SinZ_ZsPi180;
	double rc2;
	double rv2;
	double rc2rv2;
	double CosAlphaPi1802;
	double dc2;
};

double PathLength(TreeGeometry Tree, double TreeInterval, VecPoint P, double Alpha, double Azimuth, double Slope, double Aspect);
QuickCalcMem QuickCalc(double Alpha,double Slope,double Azimuth, double Aspect,TreeGeometry Tree);
QuickCalcMem QuickCalcTrunk(double Altitude,double Slope,double Azimuth, double Aspect,TreeGeometry Tree);
double SkyViewFactor(VecPoint P, bool center, double Slope, TreeGeometry Tree, double TreeInterval, int iDensity, InputData Input, double *TrunkPortion, double ***SkyMap, int* SkyMapPos);

