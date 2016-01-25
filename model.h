#pragma once

#define CYLINDRICAL_CROWN
//#define ELLIPSOIDAL_CROWN
//#define CONICAL_CROWN

#define Sigma (5.67E-8)
#define SOLAR_CONSTANT (1366.1) //W/m2

enum TreeShape{
	Cube,
	Conical,
	Cylinder,	
	Ellipsoid
};

struct TreeGeometry{
	double H;
	double Rc;
	double Rt;
	double Hb;
	double Dc;
	double Lambda;
	double G;
	double x,y,z;
	TreeShape Shape;
};


struct InputData{
	double CROWN_RADIUS, CROWN_DEPTH, TREE_HEIGHT, LAMBDA, G_CONSTANT, TRUNK_RADIUS;
	double CANOPY_ALBEDO, CANOPY_EMISSIVITY;
	double SNOW_ALBEDO_DIR, SNOW_ALBEDO_DIF, SNOW_EMISSIVITY;
	int DOY_START, DOY_END;
	
	double LATITUDE, LONGITUDE, ELEVATION, LONGITUDE_STD, DAYLIGHT_SAVING;

	double MEAN_TEMPERATURE, TEMPERATURE_DAILY_RANGE, TEMPERATURE_YEARLY_RANGE, HOTTEST_DOY;
	
	double SLOPE0, SLOPE_DEL;
	int SLOPEN;

	double ASPECT0, ASPECT_DEL;
	int ASPECTN;

	double VERY_SPARSE_SPACING, VERY_DENSE_SPACING;

	int NUMBER_OF_SPACE_GRID, NUMBER_OF_TIMESTEP_IN_A_DAY, NUMBER_OF_DENSITIES, SKYMAP_GRID_AZIMUTH, SKYMAP_GRID_ALTITUDE;
	TreeGeometry Tree;
	int SVF_OUT, PL_OUT, SKYMAP_OUT, SDIR_OUT;
	int SKY_CLEAR_CLOUDY;
	double SKY_EMISSIVITY[365];
	double tDIR[365];
	double tDIF[365];
	double RH[365];
	double kb_kD_Slope, kb_kD_Intercept;

	double DEL_T_CROWN_DENSE, DEL_T_CROWN_SPARSE;
	double DEL_T_TRUNK_DENSE, DEL_T_TRUNK_SPARSE;
	double DEL_T_RADIATION_DENSE, DEL_T_RADIATION_SPARSE;
	
	double SNOW_DEL_TEMPERATURE;
};
