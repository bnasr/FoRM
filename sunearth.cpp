#include <math.h>
#include "auxfunc.h"
#include "sunearth.h"

double SolarAltitudeAngle(double L, double Delta,double HourAngle){
	double Alpha =   asin((sin(L*Pi180)*sin(Delta *Pi180)+ cos(L *Pi180)*cos(Delta *Pi180)*cos(HourAngle*Pi180)))*180/Pi;
	return Alpha;
}
double SolarDeclinationAngle( double DOY){
	double Delta =   23.45 *sin(360.0/365 * (284 + DOY) *Pi180);

	return Delta;
}
double EquationOfTime ( double DOY){
	double B = (DOY - 81)*360/365;
	double ET = 9.87*sin(2 *B *Pi180) - 7.53 *cos(B *Pi180) -  1.5 *sin(B*Pi180);
	return ET;
}

double Day_Length(double L, double Delta){
	//(*Sunset Time*)
	double Hss = acos(-tan(L  *Pi180)* tan(Delta *Pi180))/15*180/Pi;
	//(*Sunrise Time*)
	double Hsr = -acos(-tan(L  *Pi180)* tan(Delta *Pi180))/15*180/Pi;
	//(*Day Length*)
	double DayLength = Hss - Hsr;
	return DayLength;
}
double SolarAzimuthAngle(double Delta,double HourAngle,double Alpha,double L,double AST){
	//(*Solar Azimuth Angle*)
	//(*The angle of sun's ray measured in the horizental plane from due south*)
	double Az1 = asin(cos(Delta *Pi180)* sin(HourAngle *Pi180)/cos(Alpha *Pi180 ))*180/Pi;
	double Az2 = (AST < 12*60? -180 + fabs(Az1): 180 - Az1);
	double Az = ((cos(HourAngle  *Pi180) > tan(Delta *Pi180)/ tan(L *Pi180))? Az1: Az2);
	return Az;
}
double SolarIncidenceAngle(double L,double Delta,double Beta,double AzSurface,double HourAngle){
//(*Solar Incidence Angle*)
	//(*The angle between sun's ray and the normal on a surface*)
	double Theta= acos(
			sin(L*Pi180)*sin(Delta*Pi180)*cos(Beta*Pi180)
			-cos(L*Pi180)*sin(Delta*Pi180)*sin(Beta*Pi180)*cos(AzSurface*Pi180)
			+cos(L*Pi180)*cos(Delta*Pi180)*cos(HourAngle*Pi180)*cos(Beta*Pi180)
			+sin(L*Pi180)*cos(Delta*Pi180)*cos(HourAngle*Pi180)*sin(Beta*Pi180)*cos(AzSurface*Pi180)
			+cos(Delta*Pi180)*sin(HourAngle*Pi180)*sin(Beta*Pi180)*sin(AzSurface*Pi180))*180/Pi;

	//(*Beta: Surface tilt angle from the horizental*)
	//(*AzSurface: Surface Azimuth angle:the angle betweenthe normal to the surface from ture south, westward is positive*)
	return Theta;
}
double LocalStandardTime(double DOY){
	double LST = (DOY*24*60)-((int)((DOY*24*60)/( 24*60)))*(24*60);
	return LST;
}
double AparentSolarTime(double DOY,double SL,double LL,double DS){
	double AST = LocalStandardTime(DOY) + EquationOfTime(DOY) + 4 *(SL - LL) - DS ;
	return AST;
}
double SolarHourAngle(double AST){
	double HourAngle = (AST - 12*60)/4;
	return HourAngle;
}

double AirDiurnalTemperature(double DOY,double Mean, double YearlyRange,double DailyRange,double HottestDOY){
	double Tair=Mean+0.5*YearlyRange*cos((DOY-HottestDOY)/365*2*Pi)+0.5*DailyRange*cos(DOY*2*Pi-Pi/2);

	return Tair;
}
