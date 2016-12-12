#include"Preprocess.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Preprocess::Preprocess(double Re, double res, double rho_ini, double mu_rate, double s3, double s5, bool shear, double shear_rate, int xCells_ini, int yCells_ini ,int zCells_ini):
ReynoldsMax(Re), 
resolution(res), rho(rho_ini),  
muRatio(mu_rate), s_3(s3), s_5(s5), 
isShearFlow(shear), shearRate(shear_rate),
xCells(xCells_ini), yCells(yCells_ini), zCells(zCells_ini)
{
    deduceAll();
}

//=========================== OPERATIONS ===========================

const ParamSet Preprocess::getParamSet()const{
    const double omega = 1/tau;
    const double rho_r = 1;  // normalized
    const RelaxationPar2D relax(s_2,s_3,s_5);
    ParamSet param(omega, rho_r, timestep, spacestep, relax);
    return param;
}

//=========================== OPERATOR ===========================

const bool Preprocess::operator==(const Preprocess& other)const
{
    bool exit = true;
    if(ReynoldsMax != other.getReynoldsMax()) exit = false;
    if(resolution != other.getResolution()) exit = false;
    if(rho != other.getRho()) exit = false;
    if(muRatio != other.getMuRatio()) exit = false;
    if(s_3 != other.getS_3()) exit = false;
    if(s_5 != other.getS_5()) exit = false;
    if(isShearFlow != other.getIsShearFlow()) exit = false;
    if(shearRate != other.getShearRate()) exit = false;
    if(spacestep != other.getSpacestep()) exit = false;
    if(timestep != other.getTimestep()) exit = false;   
    if(tau != other.getTau()) exit = false;
    if(c_s != other.getSoundspeed()) exit = false;
    if(nu != other.getNu()) exit = false;
    if(s_2 != other.getS2()) exit = false;

    if(xCells != other.getXCells()) exit = false;
    if(yCells != other.getYCells()) exit = false;
    if(zCells != other.getZCells()) exit = false;

    
    return exit;
}
    
///////////////////////////// PRIVATE /////////////////////////////

//=========================== OPERATIONS ===========================

void Preprocess::deduceAll(){
	calcTau();
    calcSpacestep();
    calcTimestep();
    calcSoundspeed();
	calcNu();
    calcS2();
}