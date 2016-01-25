#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <sstream>

using namespace std;

#include "auxfunc.h"

double mod(double a, double b){
	return a-((int)(a/b))*b;
}
double max(double a, double b){
	return a>b?a:b;
}
double min(double a, double b){
	return a<b?a:b;
}
struct VecPoint Vector(struct VecPoint P1,struct VecPoint P2){
	struct VecPoint V;
	V.x=P2.x - P1.x;
	V.y=P2.y - P1.y;
	V.z=P2.z - P1.z;

	return V;
}
double Norm(struct VecPoint V){
	return sqrt(V.x*V.x+V.y*V.y+V.z*V.z);
}
double DotProduct(struct VecPoint V1,struct VecPoint V2){
	double prd;
	prd=V1.x*V2.x+V1.y*V2.y+V1.z*V2.z;
	return prd;
}

double VectorAngle(struct VecPoint V1, struct VecPoint V2){
	double a=DotProduct(V1,V2)/(Norm(V1)*Norm(V2));
	 double angle= acos(min(max(a,-1),1));
	 return angle;
 }
double Distance(struct VecPoint P1, struct VecPoint P2){
	double distance=sqrt((P1.x-P2.x)*(P1.x-P2.x)+(P1.y-P2.y)*(P1.y-P2.y)+(P1.z-P2.z)*(P1.z-P2.z));
	return distance;
}
int DOY2Month(double doy){
	int monthends[12]={31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 366};
	int month=1;
	for(int m=12;m>0;m--){
		if(doy>monthends[m-1]){
			month=m+1;
			break;
		}
	}
	if (month==13)
		month=12;
	return month;
}
string SecondsToDayHourMinuteSecond(double seconds){
	int Day, Hour, Minute, Second;
	stringstream str;
	
	Day=(int)(seconds/24/60/60);
	Hour=(int)(mod(seconds,24*60*60.)/60/60);
	Minute=(int)(mod(seconds,60*60.)/60);
	Second=(int)mod(seconds,60.);
	if(Day!=0)
		str<<Day<<" days and ";
	if(Hour!=0)
		str<<Hour<<" hours and ";
	if(Minute!=0)
		str<<Minute<<" minutes and ";
	str<<Second<<" seconds.";
	return str.str();	
}
