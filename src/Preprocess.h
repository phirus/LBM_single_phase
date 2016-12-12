/// preprocesses the physical properties to LB properties that are then stored in a ParamSet

#ifndef PREPROCESS_H
#define PREPROCESS_H

#include"ParamSet.h"
#include"Constants.h"

using namespace std;

class Preprocess
{
public: 
	/// Lifecycle
	Preprocess(double Re = 10, double res = 30, double rho_ini = 1, double mu_rate = 2, double s3 = 1, double s5 = 1, bool shear = false, double shear_rate = 0, int xCells_ini = 50, int yCells_ini = 50 ,int zCells_ini = 50);

	/// operations

	// get the parameter set
	const ParamSet getParamSet()const;

	/// accessors
	inline const double getReynoldsMax()const{return ReynoldsMax;};
	inline const double getResolution()const{return resolution;};
	inline const double getRho()const{return rho;};
	inline const double getMuRatio()const{return muRatio;};
	inline const double getS_3()const{return s_3;};
	inline const double getS_5()const{return s_5;};
	inline const bool getIsShearFlow()const{return isShearFlow;};
	inline const double getShearRate()const{return shearRate;};
	
	inline const double getSpacestep()const{return spacestep;};
	inline const double getTimestep()const{return timestep;};
	inline const double getTau()const {return tau;};
	inline const double getSoundspeed()const{return c_s;};
	inline const double getNu()const{return nu;};
	inline const double getS2()const {return s_2;};

	inline const int getXCells()const{return xCells;};
	inline const int getYCells()const{return yCells;};
	inline const int getZCells()const{return zCells;};


	void setReynoldsMax(double val){ReynoldsMax = val;};

	/// operators
	const bool operator==(const Preprocess& other)const;

private:
	// given
	double ReynoldsMax ; 	/// < maximum Reynolds-Number
	double resolution; 		/// < width of bubble in cells
	double rho ;			/// < liquid density
  	double muRatio;		/// ratio of second to first viscosity mu'/mu
	double s_3, s_5;
	bool isShearFlow; 		/// < switch that determines active shear flow
	double shearRate;		/// < shear rate in (1/s)

	// stored
	int xCells, yCells, zCells;
	

    // deduced
    double spacestep;  /// < spacestep /m	
	double timestep;   /// < timestep /s

    double tau;
   	double c_s;
	double nu;
    double s_2;

	/// operations
	inline void calcTau(){tau = (resolution * MACH_MAX   * sqrt(3) / ReynoldsMax ) + 0.5;};
	inline void calcSpacestep(){spacestep = 1;}	//diameter / resolution;};
	inline void calcTimestep(){timestep = 1;} 	//spacestep / (sqrt(3) * c_s);};
	inline void calcSoundspeed(){c_s = 1.0/sqrt(3);};
	inline void calcNu(){nu = c_s * c_s * timestep * (tau - 0.5);};
	inline void calcS2(){s_2 = 1.0/( (nu * muRatio) / (c_s*c_s * timestep) + 0.5);};

	void deduceAll();
};

#endif
