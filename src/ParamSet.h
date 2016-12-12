#ifndef PARAMSET_H
#define PARAMSET_H

#include<cmath>

#include"Definitions.h"
#include<string.h>

using namespace std;

class ParamSet
{
public:
    /// Lifecycle
    ParamSet(double om = 1, double rho = 1, double t_step = 1e-3, double s_step = 1e-3, RelaxationPar2D rel = RelaxationPar2D(1,1,1)):omega(om),rho(rho),timestep(t_step),spacestep(s_step),relax(rel){}; /// < consructor

    /// get-methods, including calculations if necessary
    const double getOmega()const{return omega;};          /// < return omega, based on inter and the color field
    const RelaxationPar2D getRelaxation2D()const{return relax;};
    const double getRho()const{return rho;}; 
    const double getDeltaT()const{return timestep;};
    const double getDeltaX()const{return spacestep;};

    /// set-methods, including calculations if necessary
    void setOmega(double om){omega = om;};
    void setRelaxation(double s_2, double s_3, double s_5);

    /// overloaded == Operator
    const bool operator==(const ParamSet& other)const;

private:
    double omega;
    double rho;
    double timestep;            /// < LB timestep
    double spacestep;            /// < LB timestep
    RelaxationPar2D relax;

    const boost::array<double,4> getEverything()const;
};

#endif
