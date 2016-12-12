#ifndef LATTICE_H
#define LATTICE_H

#include<iostream>
#include<vector>
#include<omp.h>

#include"Cell.h"
#include"ParamSet.h"
#include"Boundaries.h"

/// custom typedef for the whole field of cells
typedef boost::multi_array<Cell,2> field2D;

/// contains the domain of the simulation with LB operators
class Lattice
{
public:
    /// Lifecycle
    Lattice(int x_size=10, int y_size=10,double fzero=1);
    Lattice(const Lattice& other);
    ~Lattice();

    /// operations
    /// calculations
    void equilibriumIni(); /// < replace all distribution functions with the equilibrium distribution
    void balance(double& mass, double& momentum)const; /// < monitor overall mass and momentum
    void overallRho();
    direction2D directions(int x, int y)const; /// < calculates positions of neighboring sites (periodical)
            
    /// LB steps
    void streamAll(int threads = 1); /// < streaming step
    bool collideAll(int threads = 1); /// < collision step
    void evaluateBoundaries(int threads = 1);

    /// walls
    void closedBox(); /// < initialize the Lattice (set up walls and calculate rho)
    void bottomWall(); /// < initialize the Lattice (set up walls and calculate rho)
    void genericWall(std::vector<double> x, std::vector<double> y,  const Vector2D& u_w);
    void lidDrivenCavity(const Vector2D& u_w); /// < initialize the Lattice with moving top wall
    void shearWall(const Vector2D& u_w);    /// < initialize the Lattice with moving left wall
    void setShearProfile(double gradient, double offset); /// < initialize linear shear profile according to V_y = m x + n

    /// accessors
    const DimSet2D getSize()const; /// < get the extend of the Lattice
    inline const field2D getData()const{return *data;}; /// < get the data field2D
    inline const Cell getCell(int x, int y)const{return (*data)[x][y];};  /// < get a Cell
    inline const ParamSet getParams()const{return param;}; /// < get the paramter set
    inline const Boundaries getBoundaries()const{return bound;};
    inline const array2D getF(int x, int y)const{return (*data)[x][y].getF();};          /// < get F

    void setData(const field2D& ndata, int x, int y); /// < set the data field2D (and size)
    void setCell(int y, int x, const Cell& ncell);    /// < set a Cell
    void setF(int x, int y, const array2D& nf);
    void setF(int x, int y, int index, double value);
    void setBoundaries(const Boundaries& newBound);
    inline void setParams(const ParamSet& newParam){param = newParam;}; /// < set a new parameter set

    /// operators
    Lattice& operator=(const Lattice& other);
    const bool operator==(const Lattice& other)const;

private:
    int xsize, ysize;   /// < extent of the Lattice
    field2D * data;    
    ParamSet param;     /// < set of parameters used during the simulation
    Boundaries bound;

    inline void linearIndex(int index, int& x, int& y)const;
    void streamAndBouncePull(Cell& tCell, const direction2D& dir)const; /// < internal streaming mechanism with bounce back
    const bool isBoundary(int x, int y)const;
    void buildWalls(); /// < and the Mexicans pay for it
};

/// calculates the equilibrium distribution based of a cell
const array2D eqDistro(const double rho, const Vector2D& u);

#endif // LATTICE_H
