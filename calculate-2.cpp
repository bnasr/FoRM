#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

#include <time.h>
#include <math.h>
#include <omp.h>

#include "calculate.h"
#include "auxfunc.h"
#include "svfpl.h"
#include "sunearth.h"
#include "prepost.h"

void CalculateOpenRadiation(InputData Input, double **AltitudeAngle, double **HourAngle, double **tb, double **td, double **DeclinationAngle ,double**ExtraterrestrialSolarRadiation, double Slope, double Aspect, double **SolarOpenRadiationDir, double**SolarOpenRadiationDif, ofstream *fsLOG){

	int NumberOfTimeSteps=(Input.DOY_END-Input.DOY_START)*Input.NUMBER_OF_TIMESTEP_IN_A_DAY;
	(*SolarOpenRadiationDir)=new double[NumberOfTimeSteps];
	(*SolarOpenRadiationDif)=new double[NumberOfTimeSteps];
	ReleaseLog(fsLOG, "\tOpen radiation components...");
	for(int TimeStep=0;TimeStep<NumberOfTimeSteps;TimeStep++){
		double IncidenceAngle = SolarIncidenceAngle(Input.LATITUDE,(*DeclinationAngle)[TimeStep],Slope,Aspect ,(*HourAngle)[TimeStep]);
		if((*AltitudeAngle)[TimeStep]>0){
			(*SolarOpenRadiationDir)[TimeStep]= (*tb)[TimeStep]*(*ExtraterrestrialSolarRadiation)[TimeStep]*cos(IncidenceAngle*Pi180);
			(*SolarOpenRadiationDif)[TimeStep]= (*td)[TimeStep]*(*ExtraterrestrialSolarRadiation)[TimeStep]*(1+cos(Slope*Pi180))/2*sin((*AltitudeAngle)[TimeStep]*Pi180);
		}else{
			(*SolarOpenRadiationDir)[TimeStep]=0;
			(*SolarOpenRadiationDif)[TimeStep]=0;
		}//if
	}//TimeStep
	ReleaseLog(fsLOG, "\tOpen radiation components were calculated.");
}
//void CalculateTemporal(InputData Input, double **DOY, double **Tair, double **Tair4, double **Tcan4, double **Ttrunk4, double **Tsnow4, double **DeclinationAngle, double **HourAngle,double **AltitudeAngle, double **AzimuthAngle, double **ExtraterrestrialSolarRadiation, double **tb, double **td, int **month, ofstream *fsLOG){
void CalculateTemporal(InputData Input, double **DOY, double **Tair, double **Tair4, double **Tsnow4, double **DeclinationAngle, double **HourAngle,double **AltitudeAngle, double **AzimuthAngle, double **ExtraterrestrialSolarRadiation, double **tb, double **td, int **month, ofstream *fsLOG){
	int NumberOfTimeSteps=(Input.DOY_END-Input.DOY_START)*Input.NUMBER_OF_TIMESTEP_IN_A_DAY;
	*DOY=new double[NumberOfTimeSteps];
	*Tair=new double[NumberOfTimeSteps];
	*DeclinationAngle=new double[NumberOfTimeSteps];
	*HourAngle=new double[NumberOfTimeSteps];
	*AltitudeAngle=new double[NumberOfTimeSteps];
	*AzimuthAngle=new double[NumberOfTimeSteps];
	*ExtraterrestrialSolarRadiation=new double[NumberOfTimeSteps];
	*tb=new double[NumberOfTimeSteps];
	*td=new double[NumberOfTimeSteps];
	*Tair4=new double[NumberOfTimeSteps];
	//*Tcan4=new double[NumberOfTimeSteps];
	//*Ttrunk4=new double[NumberOfTimeSteps];
	*Tsnow4=new double[NumberOfTimeSteps];
	*month=new int[NumberOfTimeSteps];
	
	ReleaseLog(fsLOG, "Temporal parameteres...");
	
	double A=Input.ELEVATION/1000.0;//
	double a0=.4237-.00821*(6-A)*(6-A);
	double a1=.5055+.00595*(6.5-A)*(6.5-A);
	double k=.2711+.01858*(2.5-A)*(2.5-A);

	// ---------------temporal variables
	for(int TimeStep=0;TimeStep<NumberOfTimeSteps;TimeStep++){
		(*DOY)[TimeStep]=mod((356+Input.DOY_START+(double)TimeStep/Input.NUMBER_OF_TIMESTEP_IN_A_DAY),365);
		(*Tair)[TimeStep]=AirDiurnalTemperature((*DOY)[TimeStep], Input.MEAN_TEMPERATURE,  Input.TEMPERATURE_YEARLY_RANGE, Input.TEMPERATURE_DAILY_RANGE, Input.HOTTEST_DOY);
		(*Tair4)[TimeStep]=pow((*Tair)[TimeStep]+273.15,4);
		//(*Tcan4)[TimeStep]=pow((*Tair)[TimeStep]+273.15,4);
		//(*Ttrunk4)[TimeStep]=pow((*Tair)[TimeStep]+273.15,4);
		(*Tsnow4)[TimeStep]=pow(min((*Tair)[TimeStep],0)+273.15+Input.SNOW_DEL_TEMPERATURE,4);
		(*DeclinationAngle)[TimeStep] = SolarDeclinationAngle((*DOY)[TimeStep]);
		double AST=AparentSolarTime((*DOY)[TimeStep],Input.LONGITUDE_STD,Input.LONGITUDE,Input.DAYLIGHT_SAVING);
		(*HourAngle)[TimeStep] = SolarHourAngle(AST);
		(*AltitudeAngle)[TimeStep] = SolarAltitudeAngle( Input.LATITUDE, (*DeclinationAngle)[TimeStep], (*HourAngle)[TimeStep]);
		(*AzimuthAngle)[TimeStep] = SolarAzimuthAngle((*DeclinationAngle)[TimeStep], (*HourAngle)[TimeStep], (*AltitudeAngle)[TimeStep], Input.LATITUDE,AST);

		(*ExtraterrestrialSolarRadiation)[TimeStep] =   SOLAR_CONSTANT *(1 + 0.033 *cos((360* ((*DOY)[TimeStep])/365)*Pi180));

		if ((*AltitudeAngle)[TimeStep]>0){
			if(Input.SKY_CLEAR_CLOUDY==0){
				(*tb)[TimeStep]=a0+a1*exp(-k/sin((*AltitudeAngle)[TimeStep]*Pi180));
				(*td)[TimeStep]=Input.kb_kD_Intercept+Input.kb_kD_Slope*(*tb)[TimeStep];
			}else{
				(*tb)[TimeStep]=Input.tDIR[(int)(*DOY)[TimeStep]];
				(*td)[TimeStep]=Input.tDIF[(int)(*DOY)[TimeStep]];
			}
			//cout<<(*tb)[TimeStep]<<'\t'<<(*td)[TimeStep]<<'\n';
		}else{
			(*tb)[TimeStep]=0;
			(*td)[TimeStep]=0;
		}
		(*month)[TimeStep]=DOY2Month((*DOY)[TimeStep]);
		//cout<<(*DOY)[TimeStep]<<'\t'<<(*tb)[TimeStep]<<'\t'<<(*td)[TimeStep]<<'\n';
	}//TimeStep
	ReleaseLog(fsLOG, "Temporal parameteres were calculated.");
}
void CalculateSpatial(InputData Input, double Slope, double**SVF_AVG, double **TrunkVF_AVG, double **MultiReflect_dir, double **MultiReflectSVF_dif, double ***SkyMap, int* SkyMapPos, ofstream *fsLOG){

	*SVF_AVG=new double[Input.NUMBER_OF_DENSITIES];
	*TrunkVF_AVG=new double[Input.NUMBER_OF_DENSITIES];
	*MultiReflect_dir=new double[Input.NUMBER_OF_DENSITIES];
	*MultiReflectSVF_dif=new double[Input.NUMBER_OF_DENSITIES];

	double ***SVF;	
	Create3DArray<double>(Input.NUMBER_OF_SPACE_GRID,Input.NUMBER_OF_SPACE_GRID,Input.NUMBER_OF_DENSITIES,&SVF,0);
	
	double ***TrunkVF;	
	Create3DArray<double>(Input.NUMBER_OF_SPACE_GRID,Input.NUMBER_OF_SPACE_GRID,Input.NUMBER_OF_DENSITIES,&TrunkVF,0);
	
	ReleaseLog(fsLOG, "Spatial parameteres...");
	for(int iDensity=0;iDensity<Input.NUMBER_OF_DENSITIES;iDensity++){
		stringstream log;
		log<<"\tiDensity=\t"<<iDensity+1<<" of "<<Input.NUMBER_OF_DENSITIES;;
		ReleaseLog(fsLOG, log.str());
		(*SVF_AVG)[iDensity]=0;
		(*TrunkVF_AVG)[iDensity]=0;
		(*MultiReflect_dir)[iDensity]=0;
		(*MultiReflectSVF_dif)[iDensity]=0;
		double TreeInterval=1./(1/Input.VERY_SPARSE_SPACING + (1.0/Input.VERY_DENSE_SPACING - 1.0/Input.VERY_SPARSE_SPACING) * iDensity / (Input.NUMBER_OF_DENSITIES-1));
		double dx=TreeInterval/Input.NUMBER_OF_SPACE_GRID;
		double dy=TreeInterval/Input.NUMBER_OF_SPACE_GRID;
		int iCenter=(int)(Input.NUMBER_OF_SPACE_GRID/2);
		int jCenter=(int)(Input.NUMBER_OF_SPACE_GRID/2);
//			***************		SPATIAL VARIABLES	at each Slope
		
		VecPoint P;
		int j;
		#pragma omp parallel for shared( SVF,TrunkVF,iDensity,SVF_AVG,TrunkVF_AVG, MultiReflect_dir,MultiReflectSVF_dif) private(j,P)
		for(int i=0;i<Input.NUMBER_OF_SPACE_GRID;i++){
			for(j=0;j<Input.NUMBER_OF_SPACE_GRID;j++){
				P.x=dx*(i-iCenter);
				P.y=dy*(j-jCenter);//(dy*(j-jCenter)>=0?dy*(j-jCenter):TreeInterval+dy*(j-jCenter));
				
				SVF[i][j][iDensity]=SkyViewFactor( P, (i==0)&(j==0), Slope,  Input.Tree,  TreeInterval, iDensity, Input, &(TrunkVF[i][j][iDensity]), SkyMap,  SkyMapPos);	
				//cout<<i<<"\t"<<j<<"\t"<<iDensity<<"\t"<<SVF[i][j][iDensity]<<"\t"<<(TrunkVF[i][j][iDensity])<<endl;

				(*SVF_AVG)[iDensity]+=(SVF)[i][j][iDensity]/(Input.NUMBER_OF_SPACE_GRID*Input.NUMBER_OF_SPACE_GRID);
				(*TrunkVF_AVG)[iDensity]+=(TrunkVF)[i][j][iDensity]/(Input.NUMBER_OF_SPACE_GRID*Input.NUMBER_OF_SPACE_GRID);
				(*MultiReflect_dir)[iDensity]+=1/(1-Input.CANOPY_ALBEDO*Input.SNOW_ALBEDO_DIR*(1-SVF[i][j][iDensity]))/(Input.NUMBER_OF_SPACE_GRID*Input.NUMBER_OF_SPACE_GRID);
				(*MultiReflectSVF_dif)[iDensity]+=SVF[i][j][iDensity]/(1-Input.CANOPY_ALBEDO*Input.SNOW_ALBEDO_DIF*(1-SVF[i][j][iDensity]))/(Input.NUMBER_OF_SPACE_GRID*Input.NUMBER_OF_SPACE_GRID);
			}//j
		}//i
	}//iDensity
	Delete3DArray<double>(Input.NUMBER_OF_SPACE_GRID,Input.NUMBER_OF_SPACE_GRID,Input.NUMBER_OF_DENSITIES,&SVF);
	Delete3DArray<double>(Input.NUMBER_OF_SPACE_GRID,Input.NUMBER_OF_SPACE_GRID,Input.NUMBER_OF_DENSITIES,&TrunkVF);
	ReleaseLog(fsLOG, "Spatial parameteres were calculated.");
}

void CalculateEnergyComponents(InputData Input, ofstream *fsLOG, ofstream *fsRadTime, ofstream *fsRadAvg, ofstream *fsPL, ofstream *fsSDIR, double **DOY,double **AltitudeAngle, double **AzimuthAngle,double Slope, double Aspect, double **SVF_AVG, double **TrunkVF_AVG, double **MultiReflect_dir, double **MultiReflectSVF_dif, double **Tair, double **Tair4, double **Tsnow4, double **SolarOpenRadiationDir, double **SolarOpenRadiationDif, int **month, int *PLMatrixPos, double ***PLMatrix){
	int NumberOfTimeSteps=(Input.DOY_END-Input.DOY_START)*Input.NUMBER_OF_TIMESTEP_IN_A_DAY;

	ReleaseLog(fsLOG, "\tEnergy components ...");
	for(int iDensity=0;iDensity<Input.NUMBER_OF_DENSITIES;iDensity++){
		stringstream log;
		log<<"\t\tiDensity=\t"<<iDensity+1<<" of "<<Input.NUMBER_OF_DENSITIES;
		ReleaseLog(fsLOG, log.str());
		double TreeInterval=1./(1/Input.VERY_SPARSE_SPACING + (1.0/Input.VERY_DENSE_SPACING - 1.0/Input.VERY_SPARSE_SPACING) * iDensity / (Input.NUMBER_OF_DENSITIES-1));
		double dx=TreeInterval/Input.NUMBER_OF_SPACE_GRID;
		double dy=TreeInterval/Input.NUMBER_OF_SPACE_GRID;
		int iCenter=(int)(Input.NUMBER_OF_SPACE_GRID/2);
		int jCenter=(int)(Input.NUMBER_OF_SPACE_GRID/2);

		double *expPathLengthAVG=new double[NumberOfTimeSteps];
		double *SFexpPathLengthAVG=new double[NumberOfTimeSteps];
		double *lcan=new double[NumberOfTimeSteps];
		double *ltrunk=new double[NumberOfTimeSteps];
		double *lsky=new double[NumberOfTimeSteps];
		double *lsnow=new double[NumberOfTimeSteps];
		double *sdir=new double[NumberOfTimeSteps];
		double *sdif=new double[NumberOfTimeSteps];
		double *ssnow=new double[NumberOfTimeSteps];
				
		int i,j;
		VecPoint P;
		double PL;
		//int PLMatrixPosOMP=0;
//------------SHORTWAVE
		#pragma omp parallel for private(i,j,P,PL) shared(expPathLengthAVG,SFexpPathLengthAVG,sdir,sdif,ssnow) 
		for(int TimeStep=0;TimeStep<NumberOfTimeSteps;TimeStep++){
			expPathLengthAVG[TimeStep]=0;
			SFexpPathLengthAVG[TimeStep]=0;
			for( i=0;i<Input.NUMBER_OF_SPACE_GRID;i++){
				for( j=0;j<Input.NUMBER_OF_SPACE_GRID;j++){
					P.x=dx*(i-iCenter);
					P.y=dy*(j-jCenter);
					PL=(*AltitudeAngle)[TimeStep]<0?0:PathLength( Input.Tree, TreeInterval, P, (*AltitudeAngle)[TimeStep], (*AzimuthAngle)[TimeStep], Slope, Aspect);
					if(Input.PL_OUT==1){
						//#pragma omp critical
						//*fsPL<<Slope<<'\t'<<Aspect<<'\t'<<iDensity+1<<'\t'<<TreeInterval<<'\t'<<TimeStep<<'\t'<<(*DOY)[TimeStep]<<'\t'<<i<<'\t'<<j<<'\t'<<P.x<<'\t'<<P.y<<'\t'<<PL<<endl;
						//cout<<"==============="<<Slope<<'\t'<<Aspect<<'\t'<<iDensity+1<<'\t'<<TreeInterval<<'\t'<<TimeStep<<'\t'<<(*DOY)[TimeStep]<<'\t'<<i<<'\t'<<j<<'\t'<<P.x<<'\t'<<P.y<<'\t'<<PL<<endl;
						int NN=(*PLMatrixPos)+TimeStep*Input.NUMBER_OF_SPACE_GRID*Input.NUMBER_OF_SPACE_GRID+i*Input.NUMBER_OF_SPACE_GRID+j;
						//cout<<"======NN\t"<<NN<<"\n";

						(*PLMatrix)[NN][0]=Slope;
						(*PLMatrix)[NN][1]=Aspect;
						(*PLMatrix)[NN][2]=iDensity+1;
						(*PLMatrix)[NN][3]=TreeInterval;
						(*PLMatrix)[NN][4]=TimeStep;
						(*PLMatrix)[NN][5]=(*DOY)[TimeStep];
						(*PLMatrix)[NN][6]=i;
						(*PLMatrix)[NN][7]=j;
						(*PLMatrix)[NN][8]=P.x;
						(*PLMatrix)[NN][9]=P.y;
						(*PLMatrix)[NN][10]=PL;
						//PLMatrixPosOMP++;
					}
					expPathLengthAVG[TimeStep]+=(double)exp(-Input.G_CONSTANT*Input.LAMBDA*PL)/(Input.NUMBER_OF_SPACE_GRID*Input.NUMBER_OF_SPACE_GRID);
					SFexpPathLengthAVG[TimeStep]=(*AltitudeAngle)[TimeStep]<0?SFexpPathLengthAVG[TimeStep]:SFexpPathLengthAVG[TimeStep]+(double)exp(-Input.G_CONSTANT*Input.LAMBDA*PL)/(Input.NUMBER_OF_SPACE_GRID*Input.NUMBER_OF_SPACE_GRID);

				}//j
			}//i
			
			if(Input.SDIR_OUT==1){
				VecPoint P1, P2;	
				P1.x=0;
				P1.y=0.50*TreeInterval;
				P2.x=0.75*TreeInterval;
				P2.y=0;
				double PL1=(*AltitudeAngle)[TimeStep]<0?0:PathLength( Input.Tree, TreeInterval, P1, (*AltitudeAngle)[TimeStep], (*AzimuthAngle)[TimeStep], Slope, Aspect);
				double PL2=(*AltitudeAngle)[TimeStep]<0?0:PathLength( Input.Tree, TreeInterval, P2, (*AltitudeAngle)[TimeStep], (*AzimuthAngle)[TimeStep], Slope, Aspect);
				double sdir1=(*MultiReflect_dir)[iDensity]*exp(-Input.G_CONSTANT*Input.LAMBDA*PL1)*(*SolarOpenRadiationDir)[TimeStep];
				double sdir2=(*MultiReflect_dir)[iDensity]*exp(-Input.G_CONSTANT*Input.LAMBDA*PL2)*(*SolarOpenRadiationDir)[TimeStep];
				double SDIR_AVG=(*MultiReflect_dir)[iDensity]*expPathLengthAVG[TimeStep]*(*SolarOpenRadiationDir)[TimeStep];

				*fsSDIR<<Slope<<"\t"<<Aspect<<"\t"<<iDensity<<"\t"<<TreeInterval<<"\t"<<TimeStep<<"\t"<<(*DOY)[TimeStep]<<"\t"<<(*SolarOpenRadiationDir)[TimeStep]<<"\t"<<SDIR_AVG<<"\t"<<PL1<<"\t"<<PL2<<"\t"<<sdir1<<"\t"<<sdir2<<"\n";
			}
			
			sdir[TimeStep]=(*MultiReflect_dir)[iDensity]*expPathLengthAVG[TimeStep]*(*SolarOpenRadiationDir)[TimeStep];
			sdif[TimeStep]=(*MultiReflectSVF_dif)[iDensity]*(*SolarOpenRadiationDif)[TimeStep];
			ssnow[TimeStep]=Input.SNOW_ALBEDO_DIR*sdir[TimeStep]+Input.SNOW_ALBEDO_DIF*sdif[TimeStep];
		}//TimeStep
		
		//	cout<<"**********"<<(*PLMatrixPos)<<"\t";
		(*PLMatrixPos)+=NumberOfTimeSteps*Input.NUMBER_OF_SPACE_GRID*Input.NUMBER_OF_SPACE_GRID;
		//scout<<"**********"<<(*PLMatrixPos)<<"\n";
		double SDIR=0, SDIF=0, SSNOW=0, SF=0, NumberOfDayTimeSteps=0;
		#pragma omp parallel for shared(AltitudeAngle, NumberOfDayTimeSteps, SFexpPathLengthAVG) reduction(+:SDIR) reduction(+:SDIF) reduction(+:SSNOW) reduction(+:SF) 
		for(int TimeStep=0;TimeStep<NumberOfTimeSteps;TimeStep++){
			SDIR+=sdir[TimeStep]/NumberOfTimeSteps;
			SDIF+=sdif[TimeStep]/NumberOfTimeSteps;
			SSNOW+=ssnow[TimeStep]/NumberOfTimeSteps;
			SF+=SFexpPathLengthAVG[TimeStep];
			if((*AltitudeAngle)[TimeStep]>0)
				NumberOfDayTimeSteps=NumberOfDayTimeSteps+1;
			//*fsRadTime<<Slope<<'\t'<<Aspect<<'\t'<<iDensity+1<<'\t'<<1./TreeInterval<<'\t'<<(*DOY)[TimeStep]<<'\t'<<TimeStep+1<<'\t'<<lcan[TimeStep]<<'\t'<<ltrunk[TimeStep]<<'\t'<<lsky[TimeStep]<<'\t'<<lsnow[TimeStep]<<'\t'<<sdir[TimeStep]<<'\t'<<sdif[TimeStep]<<'\t'<<ssnow[TimeStep]<<'\n';
		}//TimeStep
		SF=1-SF/NumberOfDayTimeSteps;
//-------------------
//------------LONGWAVE
		double DEL_T_CROWN=(Input.DEL_T_CROWN_SPARSE-Input.DEL_T_CROWN_DENSE)/(Input.DEL_T_RADIATION_SPARSE-Input.DEL_T_RADIATION_DENSE)*(SDIR-Input.DEL_T_RADIATION_DENSE)+Input.DEL_T_CROWN_DENSE;
		double DEL_T_TRUNK=(Input.DEL_T_TRUNK_SPARSE-Input.DEL_T_TRUNK_DENSE)/(Input.DEL_T_RADIATION_SPARSE-Input.DEL_T_RADIATION_DENSE)*(SDIR-Input.DEL_T_RADIATION_DENSE)+Input.DEL_T_TRUNK_DENSE;
		double DEL_T_CROWN_M, DEL_T_TRUNK_M;

		#pragma omp parallel for shared( AltitudeAngle, Input, lcan,ltrunk,lsky,lsnow) private(DEL_T_CROWN_M, DEL_T_TRUNK_M)
		for(int TimeStep=0;TimeStep<NumberOfTimeSteps;TimeStep++){

			DEL_T_CROWN_M=(*AltitudeAngle)[TimeStep]<=0?0.25*DEL_T_CROWN:DEL_T_CROWN;
			DEL_T_TRUNK_M=(*AltitudeAngle)[TimeStep]<=0?0.25*DEL_T_TRUNK:DEL_T_TRUNK;

			lcan[TimeStep]=(1-(*SVF_AVG)[iDensity]-(*TrunkVF_AVG)[iDensity])*Sigma*Input.CANOPY_EMISSIVITY*pow((*Tair)[TimeStep]+DEL_T_CROWN_M+273.15,4);
			ltrunk[TimeStep]=(*TrunkVF_AVG)[iDensity]*Sigma*Input.CANOPY_EMISSIVITY*pow((*Tair)[TimeStep]+DEL_T_TRUNK_M+273.15,4);
			lsky[TimeStep]=(*SVF_AVG)[iDensity]*Sigma*Input.SKY_EMISSIVITY[(int)(*DOY)[TimeStep]]*(*Tair4)[TimeStep];
			lsnow[TimeStep]=Sigma*Input.SNOW_EMISSIVITY*(*Tsnow4)[TimeStep];					
		}//TimeStep
			
		double LCAN=0, LTRUNK=0, LSKY=0, LSNOW=0;
		#pragma omp parallel for reduction(+:LCAN) reduction(+:LTRUNK) reduction(+:LSKY) reduction(+:LSNOW) 
		for(int TimeStep=0;TimeStep<NumberOfTimeSteps;TimeStep++){
			LCAN+=lcan[TimeStep]/NumberOfTimeSteps;
			LTRUNK+=ltrunk[TimeStep]/NumberOfTimeSteps;
			LSKY+=lsky[TimeStep]/NumberOfTimeSteps;
			LSNOW+=lsnow[TimeStep]/NumberOfTimeSteps;
			//*fsRadTime<<Slope<<'\t'<<Aspect<<'\t'<<iDensity+1<<'\t'<<1./TreeInterval<<'\t'<<(*DOY)[TimeStep]<<'\t'<<TimeStep+1<<'\t'<<lcan[TimeStep]<<'\t'<<ltrunk[TimeStep]<<'\t'<<lsky[TimeStep]<<'\t'<<lsnow[TimeStep]<<'\t'<<sdir[TimeStep]<<'\t'<<sdif[TimeStep]<<'\t'<<ssnow[TimeStep]<<'\n';
		}//TimeStep
//-------------------

		delete [] expPathLengthAVG;
		delete [] lcan;
		delete [] ltrunk;
		delete [] lsky;
		delete [] lsnow;
		delete [] sdir;
		delete [] sdif;
		delete [] ssnow;
		 
		*fsRadAvg<<Slope<<'\t'<<Aspect<<'\t'<<iDensity+1<<'\t'<<1./TreeInterval<<'\t'<<LCAN<<'\t'<<LTRUNK<<'\t'<<LSKY<<'\t'<<LSNOW<<'\t'<<SDIR<<'\t'<<SDIF<<'\t'<<SSNOW<<'\t'<<SF<<'\n';
	//cout<<"------------------------------"<<PLMatrixPosOMP<<endl;
	}//iDensity
	ReleaseLog(fsLOG, "\tEnergy components were calculated.");
}
//void DeleteTemporal( double **DOY, double **Tair, double **Tair4, double **Tcan4, double **Ttrunk4, double **Tsnow4, double **DeclinationAngle, double **HourAngle,double **AltitudeAngle, double **AzimuthAngle, double **ExtraterrestrialSolarRadiation, double **tb, double **td, int **month){
void DeleteTemporal( double **DOY, double **Tair, double **Tair4, double **Tsnow4, double **DeclinationAngle, double **HourAngle,double **AltitudeAngle, double **AzimuthAngle, double **ExtraterrestrialSolarRadiation, double **tb, double **td, int **month, ofstream *fsLOG){
	delete [] *DOY;
	delete [] *Tair;
	delete [] *DeclinationAngle;
	delete [] *HourAngle;
	delete [] *AltitudeAngle;
	delete [] *AzimuthAngle;
	delete [] *ExtraterrestrialSolarRadiation;
	delete [] *tb;
	delete [] *td;
	delete [] *Tair4;
	//delete [] *Tcan4;
	//delete [] *Ttrunk4;
	delete [] *Tsnow4;
	delete [] *month;
	ReleaseLog(fsLOG, "Temporal parameteres were cleared!");

}
void DeleteSpatial(double**SVF_AVG, double **TrunkVF_AVG, double **MultiReflect_dir, double **MultiReflectSVF_dif, ofstream *fsLOG){
	delete [] *SVF_AVG;
	delete [] *TrunkVF_AVG;
	delete [] *MultiReflect_dir;
	delete [] *MultiReflectSVF_dif; 
	ReleaseLog(fsLOG, "Spatial parameteres were cleared!");

}
/*
void DeleteSpatioTemporal(InputData Input, double ***Tcan4, double ***Ttrunk4, ofstream *fsLOG){
	int NumberOfTimeSteps=(Input.DOY_END-Input.DOY_START)*Input.NUMBER_OF_TIMESTEP_IN_A_DAY;
	Delete2DArray<double>(NumberOfTimeSteps, Input.NUMBER_OF_DENSITIES, Tcan4);
	Delete2DArray<double>(NumberOfTimeSteps, Input.NUMBER_OF_DENSITIES, Ttrunk4);
	ReleaseLog(fsLOG, "Spatiotemporal parameteres were cleared!");
}

*/