#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <cstring>

using namespace std;

#include "calculate.h"
#include "prepost.h"
#include "auxfunc.h"
#include <ctime>

inline void mySleep(double sec){
	clock_t start_time = clock();
	clock_t end_time = clock_t(sec * 1000) + start_time;
	while(clock() <end_time);
} 

void ReleaseLog(ofstream *fsLOG, string log){
		time_t rawtime;
		string strtime;

		time(&rawtime);
		strtime=ctime(&rawtime);
		strtime.replace(strtime.size()-1,1,"\t");
		*fsLOG<<strtime<<log<<endl;
		
		time(&rawtime);
		strtime=ctime(&rawtime);
		strtime.replace(strtime.size()-1,1,"\t");
		cout<<strtime<<log<<endl;
		mySleep(0.001);
}
void Error(ofstream *fsLOG, string message){
	ReleaseLog(fsLOG, message);
	exit(1);
}
void CpuTime(ofstream *fsLOG, clock_t start){
	stringstream log;		
	log<<"CPU Time:\t"<<SecondsToDayHourMinuteSecond(double(clock()-start)/CLOCKS_PER_SEC);	
	ReleaseLog(fsLOG, log.str());
}
template <class T> void ReadInput(ifstream *fsInputData, string Def, T *Param, ofstream *fsLOG){
	string DEFINITION;
	*fsInputData>>DEFINITION;
	if(strcmp(DEFINITION.c_str(),Def.c_str())){
		stringstream log;
		log<<"Error in the input file!\n"<<DEFINITION<<" is not same as "<<Def<<"." ;
		Error(fsLOG, log.str());
	}
	*fsInputData>>*Param;
	stringstream log;
	log<<"\t"<<DEFINITION<<" -->> "<<*Param;
	ReleaseLog(fsLOG, log.str());	
}
InputData ReadInputData(string FolderName, ofstream *fsLOG){
	InputData Data;
	string strInput="/Input", streCLR="/eCLR", streCLD="/eCLD", strtDIR="/tDIR", strtDIF="/tDIF", strkb_kD="/kb_kD", strRH="/RH", DEFINITION;
	ReleaseLog(fsLOG, "---");
	ReleaseLog(fsLOG, "Reading input file:");
	string InputFileName=FolderName+strInput;
	string eCLRFileName=FolderName+streCLR;
	string eCLDFileName=FolderName+streCLD;
	string tDIRFileName=FolderName+strtDIR;
	string tDIFFileName=FolderName+strtDIF;
	string kb_kDFileName=FolderName+strkb_kD;
	string RHFileName=FolderName+strRH;
	
	ifstream fskb_kD(kb_kDFileName.c_str());

	if (fskb_kD.is_open())
	{		
		fskb_kD>>Data.kb_kD_Slope;
		fskb_kD>>Data.kb_kD_Intercept;
		fskb_kD.close();
	}else
		Error(fsLOG, "kb_kD file cannot open.");
	
	ifstream fsInputData(InputFileName.c_str());

	if (fsInputData.is_open())
	{		
		ReadInput<double>(&fsInputData, "TREE_HEIGHT", &(Data.TREE_HEIGHT), fsLOG);//
		ReadInput<double>(&fsInputData, "CROWN_RADIUS", &(Data.CROWN_RADIUS), fsLOG);//
		ReadInput<double>(&fsInputData, "CROWN_DEPTH", &(Data.CROWN_DEPTH), fsLOG);//
		ReadInput<double>(&fsInputData, "TRUNK_RADIUS", &(Data.TRUNK_RADIUS), fsLOG);//
		ReadInput<double>(&fsInputData, "LAMBDA", &(Data.LAMBDA), fsLOG);//
		ReadInput<double>(&fsInputData, "G_CONSTANT", &(Data.G_CONSTANT), fsLOG);//
		ReadInput<double>(&fsInputData, "CANOPY_ALBEDO", &(Data.CANOPY_ALBEDO), fsLOG);//
		ReadInput<double>(&fsInputData, "CANOPY_EMISSIVITY", &(Data.CANOPY_EMISSIVITY), fsLOG);//
		ReadInput<double>(&fsInputData, "SNOW_ALBEDO_DIR", &(Data.SNOW_ALBEDO_DIR), fsLOG);//
		ReadInput<double>(&fsInputData, "SNOW_ALBEDO_DIF", &(Data.SNOW_ALBEDO_DIF), fsLOG);//
		ReadInput<double>(&fsInputData, "SNOW_EMISSIVITY", &(Data.SNOW_EMISSIVITY), fsLOG);//
		ReadInput<double>(&fsInputData, "ELEVATION", &(Data.ELEVATION), fsLOG);//
		ReadInput<double>(&fsInputData, "LATITUDE", &(Data.LATITUDE), fsLOG);//
		ReadInput<double>(&fsInputData, "LONGITUDE", &(Data.LONGITUDE), fsLOG);//
		ReadInput<double>(&fsInputData, "LONGITUDE_STD", &(Data.LONGITUDE_STD), fsLOG);//
		ReadInput<double>(&fsInputData, "DAYLIGHT_SAVING", &(Data.DAYLIGHT_SAVING), fsLOG);//
		ReadInput<double>(&fsInputData, "MEAN_TEMPERATURE", &(Data.MEAN_TEMPERATURE), fsLOG);//
		ReadInput<double>(&fsInputData, "TEMPERATURE_DAILY_RANGE", &(Data.TEMPERATURE_DAILY_RANGE), fsLOG);//
		ReadInput<double>(&fsInputData, "TEMPERATURE_YEARLY_RANGE", &(Data.TEMPERATURE_YEARLY_RANGE), fsLOG);//
		ReadInput<double>(&fsInputData, "HOTTEST_DOY", &(Data.HOTTEST_DOY), fsLOG);//
		ReadInput<int>(&fsInputData, "DOY_START", &(Data.DOY_START), fsLOG);//
		ReadInput<int>(&fsInputData, "DOY_END", &(Data.DOY_END), fsLOG);//
		ReadInput<int>(&fsInputData, "NUMBER_OF_SPACE_GRID", &(Data.NUMBER_OF_SPACE_GRID), fsLOG);//
		ReadInput<int>(&fsInputData, "NUMBER_OF_TIMESTEP_IN_A_DAY", &(Data.NUMBER_OF_TIMESTEP_IN_A_DAY), fsLOG);//
		ReadInput<int>(&fsInputData, "NUMBER_OF_DENSITIES", &(Data.NUMBER_OF_DENSITIES), fsLOG);//
		ReadInput<int>(&fsInputData, "SKYMAP_GRID_AZIMUTH", &(Data.SKYMAP_GRID_AZIMUTH), fsLOG);//
		ReadInput<int>(&fsInputData, "SKYMAP_GRID_ALTITUDE", &(Data.SKYMAP_GRID_ALTITUDE), fsLOG);//
		ReadInput<double>(&fsInputData, "VERY_SPARSE_SPACING", &(Data.VERY_SPARSE_SPACING), fsLOG);//
		ReadInput<double>(&fsInputData, "VERY_DENSE_SPACING", &(Data.VERY_DENSE_SPACING), fsLOG);//
		ReadInput<double>(&fsInputData, "SLOPE0", &(Data.SLOPE0), fsLOG);//
		ReadInput<double>(&fsInputData, "SLOPE_DEL", &(Data.SLOPE_DEL), fsLOG);//
		ReadInput<int>(&fsInputData, "SLOPEN", &(Data.SLOPEN), fsLOG);//
		ReadInput<double>(&fsInputData, "ASPECT0", &(Data.ASPECT0), fsLOG);//
		ReadInput<double>(&fsInputData, "ASPECT_DEL", &(Data.ASPECT_DEL), fsLOG);//
		ReadInput<int>(&fsInputData, "ASPECTN", &(Data.ASPECTN), fsLOG);//
		ReadInput<int>(&fsInputData, "SVF_OUT", &(Data.SVF_OUT), fsLOG);//
		ReadInput<int>(&fsInputData, "PL_OUT", &(Data.PL_OUT), fsLOG);//
		ReadInput<int>(&fsInputData, "SDIR_OUT", &(Data.SDIR_OUT), fsLOG);//
		ReadInput<int>(&fsInputData, "SKYMAP_OUT", &(Data.SKYMAP_OUT), fsLOG);//
		ReadInput<int>(&fsInputData, "SKY_CLEAR_CLOUDY", &(Data.SKY_CLEAR_CLOUDY), fsLOG);//

		ReadInput<double>(&fsInputData, "DEL_T_CROWN_DENSE", &(Data.DEL_T_CROWN_DENSE), fsLOG);//
		ReadInput<double>(&fsInputData, "DEL_T_CROWN_SPARSE", &(Data.DEL_T_CROWN_SPARSE), fsLOG);//
		ReadInput<double>(&fsInputData, "DEL_T_TRUNK_DENSE", &(Data.DEL_T_TRUNK_DENSE), fsLOG);//
		ReadInput<double>(&fsInputData, "DEL_T_TRUNK_SPARSE", &(Data.DEL_T_TRUNK_SPARSE), fsLOG);//
		ReadInput<double>(&fsInputData, "DEL_T_RADIATION_DENSE", &(Data.DEL_T_RADIATION_DENSE), fsLOG);//
		ReadInput<double>(&fsInputData, "DEL_T_RADIATION_SPARSE", &(Data.DEL_T_RADIATION_SPARSE), fsLOG);//
		ReadInput<double>(&fsInputData, "SNOW_DEL_TEMPERATURE", &(Data.SNOW_DEL_TEMPERATURE), fsLOG);//
		
		fsInputData.close();
	}else
		Error(fsLOG, "Input file cannot open.");


	Data.Tree.H=Data.TREE_HEIGHT;
	Data.Tree.Rc=Data.CROWN_RADIUS;
	Data.Tree.Rt=Data.TRUNK_RADIUS;
	Data.Tree.Dc=Data.CROWN_DEPTH;
	Data.Tree.Hb=Data.TREE_HEIGHT-Data.CROWN_DEPTH;
	Data.Tree.Lambda=Data.LAMBDA;
	Data.Tree.G=Data.G_CONSTANT;
	
	if(Data.SKY_CLEAR_CLOUDY==1){
		ifstream fseCLD(eCLDFileName.c_str());
		if (fseCLD.is_open()){
			for(int d=0;d<365;d++)
				fseCLD>>Data.SKY_EMISSIVITY[d];	
			fseCLD.close();
		}else
			Error(fsLOG, "eCLD file cannot open.");

		ifstream fstDIR(tDIRFileName.c_str());
		if (fstDIR.is_open()){
			for(int d=0;d<365;d++)
				fstDIR>>Data.tDIR[d];
			fstDIR.close();
		}else
		Error(fsLOG, "tDIR file cannot open.");

		ifstream fstDIF(tDIFFileName.c_str());
		if (fstDIF.is_open()){
			for(int d=0;d<365;d++)
				fstDIF>>Data.tDIF[d];
			fstDIF.close();
		}else
		Error(fsLOG, "tDIF file cannot open.");

	}else{
		ifstream fseCLR(eCLRFileName.c_str());
		if (fseCLR.is_open()){
			for(int d=0;d<365;d++)
				fseCLR>>Data.SKY_EMISSIVITY[d];
			fseCLR.close();
		}else
			Error(fsLOG, "eCLR file cannot open.");
	}
	
	ifstream fsRH(RHFileName.c_str());
	if (fsRH.is_open()){
		for(int d=0;d<365;d++)
			fsRH>>Data.RH[d];	
		fsRH.close();
	}else
		Error(fsLOG, "RH file cannot open.");

	ReleaseLog(fsLOG, "---");
	
	return Data;
}
void SpatialOut(InputData Input, double Slope, double **SVF_AVG, double **TrunkVF_AVG, double **MultiReflect_dir, double **MultiReflectSVF_dif, ofstream *fsSVF, ofstream *fsLOG){
	*fsSVF<<Slope<<'\n';
	for(int iDensity=0;iDensity<Input.NUMBER_OF_DENSITIES;iDensity++){
		*fsSVF<<(*SVF_AVG)[iDensity]<<'\t'<<(*TrunkVF_AVG)[iDensity]<<'\t'<<(*MultiReflect_dir)[iDensity]<<'\t'<<(*MultiReflectSVF_dif)[iDensity]<<'\n';
	}
	*fsSVF<<endl;
	ReleaseLog(fsLOG, "Spatial parameres were printed!");
}
void SkyMapOut(InputData Input, double ***SkyMap, int SkyMapPosLength, ofstream *fsSKYMAP, ofstream *fsLOG){
	
	ReleaseLog(fsLOG, "SkyMap is being printed...");
	for(int SkyMapPos=0;SkyMapPos<SkyMapPosLength;SkyMapPos++){
		///if((*SkyMap)[SkyMapPos][1]==3)
		*fsSKYMAP<<(*SkyMap)[SkyMapPos][0]<<'\t'<<(*SkyMap)[SkyMapPos][1]<<'\t'<<(*SkyMap)[SkyMapPos][4]<<'\t'<<(*SkyMap)[SkyMapPos][5]<<'\t'<<(*SkyMap)[SkyMapPos][6]<<endl;
	}
	Delete2DArray<double>( Input.SLOPEN*Input.NUMBER_OF_DENSITIES*Input.SKYMAP_GRID_AZIMUTH*Input.SKYMAP_GRID_ALTITUDE, 7, SkyMap);
	ReleaseLog(fsLOG, "SkyMap was printed!");
}
void PLOut(InputData Input, double ***PLMatrix, int PLMatrixLength, ofstream *fsPL, ofstream *fsLOG){
	ReleaseLog(fsLOG, "PL Matrix is being printed...");
	int PLLength=Input.SLOPEN*Input.ASPECTN*Input.NUMBER_OF_DENSITIES*(Input.DOY_END-Input.DOY_START)*Input.NUMBER_OF_TIMESTEP_IN_A_DAY*Input.NUMBER_OF_SPACE_GRID*Input.NUMBER_OF_SPACE_GRID;
	
	if(PLMatrixLength!=PLLength){
		cout<<"\n*******\nERROR in PL Length\n\n";
		cout<<PLLength<<endl;
		cout<<PLMatrixLength<<endl;
	}
	
	for(int i=0; i<PLMatrixLength;i++){
		*fsPL<<(*PLMatrix)[i][0]<<'\t';
		*fsPL<<(*PLMatrix)[i][1]<<'\t';
		*fsPL<<(*PLMatrix)[i][2]<<'\t';
		*fsPL<<(*PLMatrix)[i][3]<<'\t';
		*fsPL<<(*PLMatrix)[i][4]<<'\t';
		*fsPL<<(*PLMatrix)[i][5]<<'\t';
		*fsPL<<(*PLMatrix)[i][6]<<'\t';
		*fsPL<<(*PLMatrix)[i][7]<<'\t';
		*fsPL<<(*PLMatrix)[i][8]<<'\t';
		*fsPL<<(*PLMatrix)[i][9]<<'\t';
		*fsPL<<(*PLMatrix)[i][10]<<endl;
	}
	Delete2DArray<double>( PLLength,11, PLMatrix);
	ReleaseLog(fsLOG, "PL Matrix was printed!");
}