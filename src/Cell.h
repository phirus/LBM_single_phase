/** all distributions are provided as a Cell
* with an additional element representing the density to reduce redundance 
* when calculating rho or values like the color field or deltaRho
*/

#ifndef CELL_H
#define CELL_H

#include"Constants.h"
using namespace std;

class Cell
{
public:
    /// Lifecylce
    Cell(double fzero=1, bool solid = false); /// < construcor
    Cell(const array2D& fini); // constr

    /// operators
    const bool operator==(const Cell& other)const;

    /// operations
    void calcRho();                                         /// < calculates both densities, the velocity and delta rho
 
    /// acsessors
    inline void setF(const array2D& newF){f = newF;};
    inline void setIsSolid(bool tmp){isSolid = tmp;};
    inline void setSolidVelocity(const Vector2D& newU){if(isSolid == true) u = newU;};

    inline const array2D getF()const{return f;};
    inline const double getRho()const{return rho;};
    inline const bool getIsSolid()const{return isSolid;};
    inline const Vector2D getU()const{return u;};
   
private:
    array2D f;         /// < set of two distributions
    double rho;                    /// < rho_r = rho[0], rho_b = rho[1]
    Vector2D u;                   /// < velocity vector resultińg from the distribution
    bool isSolid;                  /// < used to mark solid cells
};

#endif // CELL_H
