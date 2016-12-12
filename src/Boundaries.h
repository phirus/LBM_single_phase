/// organize the available information on boundary conditions

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include"Definitions.h"

using namespace std;

enum boundary_type {periodic, bounceback, pressure, velocity};

class BoundaryInformation
{
public: 
    /// Lifecycle
    BoundaryInformation(boundary_type bt = periodic, double density = 0, Vector2D velocity = Vector2D()):type(bt),rho(density),u(velocity){};
    BoundaryInformation(int bt, double density = 0, Vector2D velocity = Vector2D()):type(static_cast<boundary_type>(bt)),rho(density),u(velocity){};
    
    /// accessors
    inline const boundary_type getType()const{return type;};
    inline const double getRho()const{return rho;};
    inline const Vector2D getVelocity()const{return u;};

    /// operators
    const bool operator==(const BoundaryInformation& other)const{return (type == other.getType())&& (rho == other.getRho())&&(u == other.getVelocity());};

private:
    boundary_type type;
    double rho;                    
    Vector2D u;
};

struct Boundaries
{
    BoundaryInformation north;
    BoundaryInformation south;
    BoundaryInformation east;
    BoundaryInformation west;
    const bool operator==(const Boundaries& o)const{return (north == o.north)&&(south == o.south)&&(east == o.east)&&(west == o.west);};
};

#endif 