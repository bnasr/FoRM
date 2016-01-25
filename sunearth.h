#pragma once

double SolarIncidenceAngle(double L,double Delta,double Beta,double AzSurface,double HourAngle);
double SolarAzimuthAngle(double Delta,double HourAngle,double Alpha,double L,double AST);
double Day_Length(double L,double Delta);
double SolarAltitudeAngle(double L,double Delta,double HourAngle);
double SolarDeclinationAngle(double DOY);
double EquationOfTime(double DOY);
double LocalStandardTime(double DOY);
double AparentSolarTime(double DOY,double SL,double LL,double DS);
double SolarHourAngle(double AST);
double AirDiurnalTemperature(double DOY,double Mean, double YearlyRange,double DailyRange,double HottestDOY);
