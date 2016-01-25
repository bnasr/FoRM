#pragma once
#include <string>
using namespace std;

#define Pi (3.1415926)
#define Pi180 (0.017453293)

struct VecPoint{
	double x;
	double y;
	double z;
};

double VectorAngle(struct VecPoint V1, struct VecPoint V2);
double DotProduct(struct VecPoint V1,struct VecPoint V2);
double Distance(struct VecPoint P1, struct VecPoint P2);
double Norm(struct VecPoint V);
struct VecPoint Vector(struct VecPoint P1,struct VecPoint P2);
double min(double a, double b);
double max(double a, double b);
double mod(double a, double b);
int DOY2Month(double doy);
string SecondsToDayHourMinuteSecond(double seconds);
template <class T> void Create3DArray(int N1,int N2,int N3,T ****Array, T Value){

	int i,j,k;
	(*Array)= new T**[N1];

	for(i=0;i<N1;i++) 
		(*Array)[i]= new T*[N2];

	for(i=0;i<N1;i++)
		for(j=0;j<N2;j++)
			(*Array)[i][j]= new T[N3];

	for(i=0;i<N1;i++)
		for (j=0;j<N2;j++)
			for (k=0;k<N3;k++)
				(*Array)[i][j][k]=Value;
}
template <class T> void Create2DArray(int N1,int N2,T ***Array, T Value){
	int i,j;
	(*Array)= new T*[N1];

	for(i=0;i<N1;i++) 
		(*Array)[i]= new T[N2];

	for(i=0;i<N1;i++)
		for (j=0;j<N2;j++)
			(*Array)[i][j]=Value;
}
template <class T> void Delete3DArray(int N1,int N2,int N3,T ****Array){

	int i,j;
	for(i=0;i<N1;i++)
		for(j=0;j<N2;j++)
			delete [] (*Array)[i][j];

	for(i=0;i<N1;i++) 
		delete [] (*Array)[i];

	delete [] (*Array);
}
template <class T> void Delete2DArray(int N1,int N2,T ***Array){

	int i;
	for(i=0;i<N1;i++) 
		delete [] (*Array)[i];

	delete [] (*Array);
}
