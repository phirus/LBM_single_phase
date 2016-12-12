#include"ParamSet.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== OPERATIONS ===========================


const boost::array<double,4> ParamSet::getEverything()const{
    boost::array<double,4> pinkie;
    pinkie[0] = omega;
    pinkie[1] = rho;
    pinkie[2] = timestep;
    pinkie[3] = spacestep;
    
    return pinkie;
}

void ParamSet::setRelaxation(double s_2, double s_3, double s_5)
{
    relax.s_2 = s_2;
    relax.s_3 = s_3;
    relax.s_5 = s_5;
}

const bool ParamSet::operator==(const ParamSet& other)const{
    bool control = true;
    {
        boost::array<double,4> foo, bar;
        foo = getEverything();
        bar = other.getEverything();

        if(foo != bar) control = false;
    }
    
    {
        RelaxationPar2D brelax = other.getRelaxation2D();

        if(relax.s_2 != brelax.s_2 || relax.s_3 != brelax.s_3 || relax.s_5 != brelax.s_5) control = false;
    }
    return control;
}
