#include "Cell.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Cell::Cell(double fzero, bool solid):isSolid(solid)
{
    f[0] = fzero;
    for (int i=1; i<=8; i++)
    {
        f[i]=0;
    }
    rho = 0;

    u.x = 0;
    u.y = 0;
}

// like a copy constructor for the bulk phase
Cell::Cell(const array2D& fini):isSolid(false)
{
    f = fini;
    rho = 0;
    u.x = 0;
    u.y = 0;
}

//=========================== OPERATORS ===========================

const bool Cell::operator==(const Cell& other)const {
    bool exit = true;
    if (isSolid != other.getIsSolid()) exit = false;

    if(rho != other.getRho() ) exit = false;

    array2D fOther = other.getF();
    for(int q=0;q<9;q++)
    {
        if(f[q] != fOther[q]) exit = false;
    }
    return exit;
}

//=========================== OPERATIONS ===========================

void Cell::calcRho()
{
    // initialize density
    rho = 0;

    if (isSolid == false)
    {
        // initialize velocities
        u.x = 0;
        u.y = 0;

        // iterate
        for (int i=0; i<9; i++)
        {
            rho += f[i];
            u.x += ( f[i]) * DIRECTION_2D[i].x;
            u.y += ( f[i]) * DIRECTION_2D[i].y;
        }

        if(rho >0) 
        {
            u.x /= rho;
            u.y /= rho;
        }
    }
}