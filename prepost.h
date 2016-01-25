#pragma once
#include "model.h"
#include <time.h>

InputData ReadInputData(string ProjectName, ofstream *fsLOG);
void Error(ofstream *fsLOG, string log);
void ReleaseLog(ofstream * fsLOG, string log);
void CpuTime(ofstream *fsLOG, clock_t start);
void SpatialOut(InputData Input, double Slope, double **SVF_AVG, double **TrunkVF_AVG, double **MultiReflect_dir, double **MultiReflectSVF_dif, ofstream *fsSVF, ofstream *fsLOG);
void SkyMapOut(InputData Input, double ***SkyMap, int SkyMapPosLength, ofstream *fsSKYMAP, ofstream *fsLOG);
void PLOut(InputData Input, double ***PLMatrix, int PLMatrixLength, ofstream *fsPL, ofstream *fsLOG);
