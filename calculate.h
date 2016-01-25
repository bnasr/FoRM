#pragma once
#include <fstream>
#include <iostream>
#include "model.h"

using namespace std;
//void CalculateTemporal(InputData Input, double **DOY, double **Tair, double **Tair4, double **Tcan4, double **Ttrunk4, double **Tsnow4, double **DeclinationAngle, double **HourAngle,double **AltitudeAngle, double **AzimuthAngle, double **ExtraterrestrialSolarRadiation, double **tb, double **td, int **month, ofstream *fsLOG);
void CalculateOpenRadiation(InputData Input, double **AltitudeAngle, double **HourAngle, double **tb, double **td, double **DeclinationAngle ,double**ExtraterrestrialSolarRadiation, double Slope, double Aspect, double **SolarOpenRadiationDir, double**SolarOpenRadiationDif, ofstream *fsLOG);
void CalculateTemporal(InputData Input, double **DOY, double **Tair, double **Tair4, double **Tsnow4, double **DeclinationAngle, double **HourAngle,double **AltitudeAngle, double **AzimuthAngle, double **ExtraterrestrialSolarRadiation, double **tb, double **td, int **month, ofstream *fsLOG);
void CalculateSpatial(InputData Input, double Slope, double**SVF_AVG, double **TrunkVF_AVG, double **MultiReflect_dir, double **MultiReflectSVF_dif, double ***SkyMap, int* SkyMapPos, ofstream *fsLOG);
void CalculateEnergyComponents(InputData Input, ofstream *fsLOG, ofstream *fsRadTime,ofstream *fsRadAvg, ofstream *fsPL, ofstream *fsSDIR, double **DOY,double **AltitudeAngle, double **AzimuthAngle,double Slope, double Aspect, double **SVF_AVG, double **TrunkVF_AVG, double **MultiReflect_dir, double **MultiReflectSVF_dif, double **Tair, double **Tair4, double **Tsnow4, double **SolarOpenRadiationDir, double **SolarOpenRadiationDif, int **month, int *PLMatrixPos, double ***PLMatrix);
//void DeleteTemporal( double **DOY, double **Tair, double **Tair4, double **Tcan4, double **Ttrunk4, double **Tsnow4, double **DeclinationAngle, double **HourAngle,double **AltitudeAngle, double **AzimuthAngle, double **ExtraterrestrialSolarRadiation, double **tb, double **td, int **month);
void DeleteTemporal( double **DOY, double **Tair, double **Tair4, double **Tsnow4, double **DeclinationAngle, double **HourAngle,double **AltitudeAngle, double **AzimuthAngle, double **ExtraterrestrialSolarRadiation, double **tb, double **td, int **month, ofstream *fsLOG);
void DeleteSpatial(double**SVF_AVG, double **TrunkVF_AVG, double **MultiReflect_dir, double **MultiReflectSVF_dif, ofstream *fsLOG);
//void DeleteSpatioTemporal(InputData Input, double ***Tcan4, double ***Ttrunk4, ofstream *fsLOG);
double DewPointTemperature(double T, double RH);