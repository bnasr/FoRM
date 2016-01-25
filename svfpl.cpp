#include <iostream>
#include <fstream>
using namespace std;
#include <math.h>
#include "auxfunc.h"
#include "svfpl.h"
#include "auxfunc.h"
#include "model.h"

QuickCalcMem QuickCalc(double Altitude,double Slope,double Azimuth, double Aspect,TreeGeometry Tree){
	QuickCalcMem QCM;
	
	QCM.CosBetaPi180=cos(Slope*Pi180);
	QCM.SinBetaPi180=sin(Slope*Pi180);
	QCM.CosAlphaPi180=cos(Altitude*Pi180);
	QCM.SinAlphaPi180=sin(Altitude*Pi180);
	QCM.Cos2AlphaPi180=cos(2*Altitude*Pi180);
	QCM.CosZ_ZsPi180=cos((Azimuth-Aspect)*Pi180);
	QCM.SinZ_ZsPi180=sin((Azimuth-Aspect)*Pi180);
	QCM.rc2=Tree.Rc*Tree.Rc;
	QCM.rv2=(Tree.H-Tree.Hb)*(Tree.H-Tree.Hb)/4.0;
	QCM.rc2rv2=QCM.rc2*QCM.rv2;
	QCM.CosAlphaPi1802=QCM.CosAlphaPi180*QCM.CosAlphaPi180;
	QCM.dc2=(Tree.H-Tree.Hb)*(Tree.H-Tree.Hb);

	return QCM;
}

QuickCalcMem QuickCalcTrunk(double Altitude,double Slope,double Azimuth, double Aspect,TreeGeometry Tree){
	QuickCalcMem QCM;
	
	QCM.CosBetaPi180=cos(Slope*Pi180);
	QCM.SinBetaPi180=sin(Slope*Pi180);
	QCM.CosAlphaPi180=cos(Altitude*Pi180);
	QCM.SinAlphaPi180=sin(Altitude*Pi180);
	QCM.Cos2AlphaPi180=cos(2*Altitude*Pi180);
	QCM.CosZ_ZsPi180=cos((Azimuth-Aspect)*Pi180);
	QCM.SinZ_ZsPi180=sin((Azimuth-Aspect)*Pi180);
	QCM.rc2=Tree.Rt*Tree.Rt;
	QCM.rv2=(Tree.Hb)*(Tree.Hb)/4.0;
	QCM.rc2rv2=QCM.rc2*QCM.rv2;
	QCM.CosAlphaPi1802=QCM.CosAlphaPi180*QCM.CosAlphaPi180;
	QCM.dc2=(Tree.Hb)*(Tree.Hb);

	return QCM;
}

double PathLength(TreeGeometry Tree,double TreeInterval,VecPoint P,double Altitude,double Azimuth,double Slope,double Aspect){
	double PL=0;
	int NumberOfComputingTreesInOneDirection;
	int i,j;
	double xs=P.x;
	double ys=P.y;

	QuickCalcMem QCM=QuickCalc(Altitude, Slope, Azimuth,  Aspect, Tree);
	QuickCalcMem QCMTrunk=QuickCalcTrunk(Altitude, Slope, Azimuth,  Aspect, Tree);
	NumberOfComputingTreesInOneDirection=(int)min(max((int)(Tree.H/(TreeInterval*(QCM.CosBetaPi180*tan(Altitude*Pi180)+QCM.SinBetaPi180))),7),30);
	/*
	if(NumberOfComputingTreesInOneDirection==20)
		return 2*Tree.Rc*Tree.H/(TreeInterval*(QCM.CosBetaPi180*tan(Altitude*Pi180)+QCM.SinBetaPi180));
	*/
	for(i=-NumberOfComputingTreesInOneDirection;i<=NumberOfComputingTreesInOneDirection;i++){
		for(j=-NumberOfComputingTreesInOneDirection;j<=NumberOfComputingTreesInOneDirection;j++){
			Tree.x=i*TreeInterval;
			Tree.y=j*TreeInterval*QCM.CosBetaPi180;
			Tree.z=j*TreeInterval*QCM.SinBetaPi180;
			double a,b,c,delta;
			double PL1=0, PLtrunk=0;
			
// Trunk calculations
			//if(Tree.Rt>0){
				a=(QCMTrunk.CosAlphaPi1802)/QCMTrunk.rc2;	
				b=(2*QCMTrunk.CosAlphaPi180*(QCMTrunk.CosZ_ZsPi180*(Tree.y-ys*QCMTrunk.CosBetaPi180)+(xs-Tree.x)*QCMTrunk.SinZ_ZsPi180))/QCMTrunk.rc2;
				c=(-QCMTrunk.rc2+(xs-Tree.x)*(xs-Tree.x)+Tree.y*Tree.y+ys*QCMTrunk.CosBetaPi180*(-2*Tree.y+ys*QCMTrunk.CosBetaPi180))/QCMTrunk.rc2;
				delta=b*b-4*a*c;
				if(delta>=0){
					double l1=(-b-sqrt(delta))/2/a;
					double l2=(-b+sqrt(delta))/2/a;
					double zl1=l1*QCMTrunk.SinAlphaPi180 + ys*QCMTrunk.SinBetaPi180-Tree.z;
					double zl2=l2*QCMTrunk.SinAlphaPi180 + ys*QCMTrunk.SinBetaPi180-Tree.z;
					double zl1b=max(min(zl1,Tree.Hb),0);
					double zl2b=max(min(zl2,Tree.Hb),0);
					if(zl1!=zl2)
						PL1=(zl2b-zl1b)/(zl2-zl1)*(l2-l1);
					else
						PL1=((zl1<0)||(zl1>Tree.Hb)?0:l2-l1);
				}//if
//				else
//					PL1=0;
			//}//if
			if(PL1!=0)
				PLtrunk=DARK_BEAM_PL;
				//return DARK_BEAM_PL;
			//else{
		
#ifdef CONICAL_CROWN
				a=((QCM.dc2-QCM.rc2+(QCM.dc2+QCM.rc2)*QCM.Cos2AlphaPi180))/(2*QCM.rc2);
				b=(2*(QCM.dc2*QCM.CosZ_ZsPi180*QCM.CosAlphaPi180*(Tree.y-ys*QCM.CosBetaPi180)+QCM.dc2*(xs-Tree.x)*QCM.CosAlphaPi180*QCM.SinZ_ZsPi180+QCM.rc2*QCM.SinAlphaPi180*(Tree.H-ys*QCM.SinBetaPi180)))/QCM.rc2;
				c=(-Tree.H*Tree.H*QCM.rc2+QCM.dc2*((xs-Tree.x)*(xs-Tree.x)+Tree.y*Tree.y)+ys*(-2*QCM.dc2*Tree.y*QCM.CosBetaPi180+QCM.dc2*ys*QCM.CosBetaPi180*QCM.CosBetaPi180+QCM.rc2*QCM.SinBetaPi180*(2*Tree.H-ys*QCM.SinBetaPi180)))/QCM.rc2;
				delta=b*b-4*a*c;
				if(delta>=0){
					double l1=(-b-sqrt(delta))/2/a;
					double l2=(-b+sqrt(delta))/2/a;
					double zl1=l1*QCM.SinAlphaPi180 + ys*QCM.SinBetaPi180-Tree.z;
					double zl2=l2*QCM.SinAlphaPi180 + ys*QCM.SinBetaPi180-Tree.z;
					if((zl1>=Tree.Hb&&zl1<=Tree.H)&&(zl2>=Tree.Hb&&zl2<=Tree.H)){
						PL1=fabs(l2-l1);
					}
					else if((zl1>=Tree.Hb&&zl1<=Tree.H)&&(!(zl2>=Tree.Hb&&zl2<=Tree.H))){
						PL1=(zl1-Tree.z-Tree.Hb)/QCM.SinAlphaPi180;
					}else if((!(zl1>=Tree.Hb&&zl1<=Tree.H))&&(zl2>=Tree.Hb&&zl2<=Tree.H)){
						PL1=(zl2-Tree.z-Tree.Hb)/QCM.SinAlphaPi180;
					}else{
						PL1=0;
					}//if
				}else
					PL1=0;	
				if(PL1<0) PL1=0;
#endif
#ifdef ELLIPSOIDAL_CROWN
				c=(-1+((xs-Tree.x)*(xs-Tree.x)+Tree.y*Tree.y+ys*QCM.CosBetaPi180*(-2*Tree.y+ys*QCM.CosBetaPi180))/(QCM.rc2)+(-Tree.H+0.5*Tree.Dc+ys*QCM.SinBetaPi180)*(-Tree.H+0.5*Tree.Dc+ys*QCM.SinBetaPi180)/(QCM.rv2));
				b=(2*(QCM.rv2*QCM.CosZ_ZsPi180*QCM.CosAlphaPi180*(Tree.y-ys*QCM.CosBetaPi180)+QCM.rv2*(xs-Tree.x)*QCM.CosAlphaPi180*QCM.SinZ_ZsPi180+QCM.rc2*QCM.SinAlphaPi180*(-Tree.H+0.5*Tree.Dc+ys*QCM.SinBetaPi180)))/(QCM.rc2rv2);
				a=((QCM.rc2+QCM.rv2+(-QCM.rc2+QCM.rv2)*QCM.Cos2AlphaPi180))/(2*QCM.rc2rv2);	
				delta=b*b-4*a*c;
				PL1 = (delta > 0? pow(delta,.5)/a: 0);
#endif
#ifdef CYLINDRICAL_CROWN
				a=(QCM.CosAlphaPi1802)/QCM.rc2;	
				b=(2*QCM.CosAlphaPi180*(QCM.CosZ_ZsPi180*(Tree.y-ys*QCM.CosBetaPi180)+(xs-Tree.x)*QCM.SinZ_ZsPi180))/QCM.rc2;
				c=(-QCM.rc2+(xs-Tree.x)*(xs-Tree.x)+Tree.y*Tree.y+ys*QCM.CosBetaPi180*(-2*Tree.y+ys*QCM.CosBetaPi180))/QCM.rc2;
				delta=b*b-4*a*c;
				if(delta>=0){
					double l1=(-b-sqrt(delta))/2/a;
					double l2=(-b+sqrt(delta))/2/a;
					double zl1=l1*QCM.SinAlphaPi180 + ys*QCM.SinBetaPi180-Tree.z;
					double zl2=l2*QCM.SinAlphaPi180 + ys*QCM.SinBetaPi180-Tree.z;
					double zl1b=max(min(zl1,Tree.H),Tree.Hb);
					double zl2b=max(min(zl2,Tree.H),Tree.Hb);
					if(zl1!=zl2)
						PL1=(zl2b-zl1b)/(zl2-zl1)*(l2-l1);
					else
						PL1=((zl1<Tree.Hb)||(zl1>Tree.H)?0:l2-l1);
				}else
					PL1=0;		
#endif
			//}
			PL+=(PL1+PLtrunk);
		}
	}
	return PL;
}

double SkyViewFactor(VecPoint P, bool center, double Slope, TreeGeometry Tree, double TreeInterval, int iDensity, InputData Input, double *TrunkPortion, double ***SkyMap, int* SkyMapPos){
	double OpenSkyPortion=0;
	int TrunkBeams=0;
//	cout<<iDensity<<"\t"<<P.x<<"\t"<<P.y<<endl;
	for(int iAzimuth=0;iAzimuth<Input.SKYMAP_GRID_AZIMUTH;iAzimuth++){
		double Azimuth=iAzimuth*360./Input.SKYMAP_GRID_AZIMUTH;
		//double HorizoneAltitude=((Azimuth>=90&&Azimuth<=270)?+1:-1)*acos(1/sqrt(1+cos(Azimuth*Pi180)*cos(Azimuth*Pi180)*tan(Slope*Pi180)*tan(Slope*Pi180)));
		for(int iAltitude=0;iAltitude<Input.SKYMAP_GRID_ALTITUDE;iAltitude++){
			double Altitude=-90+iAltitude*180.0/Input.SKYMAP_GRID_ALTITUDE;
			//double PL=(Altitude<HorizoneAltitude?1000:PathLength( Tree, TreeInterval, P, Altitude, Azimuth, Slope, 0));
			double PL=PathLength( Tree, TreeInterval, P, Altitude, Azimuth, Slope, 0);
			//cout<<Altitude<<"\t"<<PL<<endl;
			OpenSkyPortion+=exp(-Input.G_CONSTANT*Input.LAMBDA*PL);

			if(PL>DARK_BEAM_PL-1)
				TrunkBeams++;
			if((Input.SKYMAP_OUT==1)&(center)){
				(*SkyMap)[*SkyMapPos][0]=Slope;
				(*SkyMap)[*SkyMapPos][1]=iDensity;
				(*SkyMap)[*SkyMapPos][2]=iAzimuth;
				(*SkyMap)[*SkyMapPos][3]=iAltitude;
				(*SkyMap)[*SkyMapPos][4]=Azimuth;
				(*SkyMap)[*SkyMapPos][5]=Altitude;
				(*SkyMap)[*SkyMapPos][6]=PL;
				(*SkyMapPos)++;
			}
		}
	}
	//cout<<OpenSkyPortion<<"\t"<<TrunkBeams<<endl;

	*TrunkPortion=(double)(TrunkBeams/(double)(Input.SKYMAP_GRID_AZIMUTH*Input.SKYMAP_GRID_ALTITUDE));
	double svf=OpenSkyPortion/(Input.SKYMAP_GRID_AZIMUTH*Input.SKYMAP_GRID_ALTITUDE);
	//cout<<" -----------------------"<<svf<<"\t"<<*TrunkPortion<<"\t"<<endl;
	return svf;
}


