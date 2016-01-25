#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "prepost.h"
#include "auxfunc.h"
#include "calculate.h"
#include "svfpl.h"

int main(int argc, char*argv[]){


	string ProjectFolder, strRadAvg="/RadiationAverage",strRadTime="/RadiationTimeSeries",strPL="/PL",strSDIR="/SDIR", strSKYMAP="/SKYMAP", strLOG="/log.txt", strSVF="/SVF";
	double *DOY, *Tair, *Tair4, *Tsnow4, *DeclinationAngle, *HourAngle, *AltitudeAngle, *AzimuthAngle, *ExtraterrestrialSolarRadiation, *tb, *td;
	int *month;
	double *SVF_AVG, *TrunkVF_AVG, *MultiReflect_dir, *MultiReflectSVF_dif;
	double **SkyMap;
	double **PLMatrix;
	double *SolarOpenRadiationDir, *SolarOpenRadiationDif;
	stringstream log;
	time_t starttime, endtime;
	clock_t start;
	
	int SkyMapPos=0;
	int PLMatrixPos=0;
	
	time(&starttime);
	start = clock();
	
	if(argc==1)
		ProjectFolder="";
	else
		ProjectFolder=argv[1];
		
	if (ProjectFolder[ProjectFolder.length()-1] == '\r')
		ProjectFolder.erase(ProjectFolder.size() - 1);

	ofstream fsLOG ((ProjectFolder+strLOG).c_str());
	if(!fsLOG.is_open())
		Error(&fsLOG,"fsLOG is not open!");
		
	log.str("");		
	log<<"\n----------------------\nProject: "<<ProjectFolder;		
	ReleaseLog(&fsLOG, log.str());

	InputData Input=ReadInputData(ProjectFolder, &fsLOG);
	//VecPoint p;
	//p.x=0;
	//p.y=5;
	//p.z=0;
	//double pl=PathLength(Input.Tree,5,p, 0, 0,0,0);

	ofstream fsSVF;	
	if(Input.SVF_OUT==1){
		fsSVF.open((ProjectFolder+strSVF).c_str());
		if(!fsSVF.is_open())
			Error(&fsLOG,"fsSVF is not open!");
		else
			fsSVF<<Input.SLOPEN<<'\t'<<Input.NUMBER_OF_DENSITIES<<endl;		
	}//if
	ofstream fsPL ((ProjectFolder+strPL).c_str());
	ofstream fsSDIR ((ProjectFolder+strSDIR).c_str());
	ofstream fsSKYMAP ((ProjectFolder+strSKYMAP).c_str());
	ofstream fsRadTime ((ProjectFolder+strRadTime).c_str());
	ofstream fsRadAvg ((ProjectFolder+strRadAvg).c_str());
	
	if (!fsPL.is_open())
		Error(&fsLOG,"fsPL is not open!");
	if (!fsSDIR.is_open())
		Error(&fsLOG,"fsSDIR is not open!");
	if (!fsSKYMAP.is_open())
		Error(&fsLOG,"fsSKYMAP is not open!");
	if (!fsRadTime.is_open())
		Error(&fsLOG,"fsRadTime is not open!");
	else if(!fsRadAvg.is_open())
		Error(&fsLOG,"fsRadAvg is not open!");
		
	CalculateTemporal( Input, &DOY, &Tair, &Tair4, &Tsnow4, &DeclinationAngle, &HourAngle, &AltitudeAngle, &AzimuthAngle, &ExtraterrestrialSolarRadiation, &tb, &td, &month, &fsLOG);
	CpuTime(&fsLOG, start);
//	CalculateSpatioTemporal( Input, &AltitudeAngle, &DOY, &Tair, &Tcan4, &Ttrunk4, &fsLOG);
	//CpuTime(&fsLOG, start);

	if(Input.SKYMAP_OUT==1)
		Create2DArray<double>( Input.SLOPEN*Input.NUMBER_OF_DENSITIES*Input.SKYMAP_GRID_AZIMUTH*Input.SKYMAP_GRID_ALTITUDE, 7, &SkyMap, 0);
	if(Input.PL_OUT==1)
		Create2DArray<double>( Input.SLOPEN*Input.ASPECTN*Input.NUMBER_OF_DENSITIES*(Input.DOY_END-Input.DOY_START)*Input.NUMBER_OF_TIMESTEP_IN_A_DAY*Input.NUMBER_OF_SPACE_GRID*Input.NUMBER_OF_SPACE_GRID, 11, &PLMatrix, 0);

	for(int iSlope=0;iSlope<Input.SLOPEN;iSlope++){
		double Slope=Input.SLOPE0+iSlope*Input.SLOPE_DEL;
		ReleaseLog(&fsLOG, "-------------");		
		log.str("");		
		log<<"Slope="<<Slope;		
		ReleaseLog(&fsLOG, log.str());
		CalculateSpatial( Input, Slope, &SVF_AVG, &TrunkVF_AVG, &MultiReflect_dir, &MultiReflectSVF_dif, &SkyMap, &SkyMapPos, &fsLOG);
		CpuTime(&fsLOG, start);

		if(Input.SVF_OUT==1)
			SpatialOut(Input, Slope, &SVF_AVG, &TrunkVF_AVG, &MultiReflect_dir, &MultiReflectSVF_dif, &fsSVF, &fsLOG);
		for(int iAspect=0;iAspect<Input.ASPECTN;iAspect++){
			double Aspect=Input.ASPECT0+iAspect*Input.ASPECT_DEL;
			ReleaseLog(&fsLOG, "\t---");			
			log.str("");			
			log<<"\tAspect="<<Aspect;			
			ReleaseLog(&fsLOG, log.str());
			
			CalculateOpenRadiation(Input, &AltitudeAngle, &HourAngle, &tb, &td, &DeclinationAngle ,&ExtraterrestrialSolarRadiation,  Slope,  Aspect, &SolarOpenRadiationDir, &SolarOpenRadiationDif, &fsLOG);
			CpuTime(&fsLOG, start);
			//cout<<"1============="<<PLMatrixPos<<endl;
			CalculateEnergyComponents(Input, &fsLOG, &fsRadTime, &fsRadAvg, &fsPL, &fsSDIR, &DOY,&AltitudeAngle, &AzimuthAngle, Slope, Aspect, &SVF_AVG, &TrunkVF_AVG, &MultiReflect_dir, &MultiReflectSVF_dif, &Tair, &Tair4, &Tsnow4, &SolarOpenRadiationDir, &SolarOpenRadiationDif, &month, &PLMatrixPos, &PLMatrix);
			//cout<<"2============="<<PLMatrixPos<<endl;
			CpuTime(&fsLOG, start);

			delete [] SolarOpenRadiationDir;
			delete [] SolarOpenRadiationDif;
		}//Aspect	
		DeleteSpatial(&SVF_AVG, &TrunkVF_AVG, &MultiReflect_dir, &MultiReflectSVF_dif, &fsLOG);
	}//Slope
	if(Input.SKYMAP_OUT==1)
		SkyMapOut( Input, &SkyMap,  SkyMapPos, &fsSKYMAP, &fsLOG);

	if(Input.PL_OUT==1)
		PLOut( Input, &PLMatrix,  PLMatrixPos, &fsPL, &fsLOG);

	//DeleteSpatioTemporal(Input, &Tcan4, &Ttrunk4, &fsLOG);
	DeleteTemporal(&DOY, &Tair, &Tair4, &Tsnow4, &DeclinationAngle, &HourAngle, &AltitudeAngle, &AzimuthAngle, &ExtraterrestrialSolarRadiation, &tb, &td, &month, &fsLOG);
	ReleaseLog(&fsLOG, "---");
	time(&endtime);
	log.str("");		
	log<<"Done in "<<SecondsToDayHourMinuteSecond(difftime (endtime,starttime))<<"\n----------------------";	
	ReleaseLog(&fsLOG, log.str());
	ReleaseLog(&fsLOG, "---");

	fsPL.close();
	fsRadTime.close();
	fsRadAvg.close();
	fsLOG.close();
	fsSVF.close();
	//getchar();
	return 0;
}//main()
