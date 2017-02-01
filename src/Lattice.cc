#include"Lattice.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Lattice::Lattice(int x_size, int y_size,double fzero):
xsize(x_size)
,ysize(y_size)
,data(new field2D(boost::extents[xsize][ysize]))
,param()
,bound()
{
    for (int x = 0; x<xsize; x++)
    {
        for (int y=0; y<ysize; y++)
        {
            (*data)[x][y] = Cell(fzero);
        }
    }
}

Lattice::Lattice(const Lattice& other):
xsize(other.getSize()[0])
,ysize(other.getSize()[1])
,data(new field2D(boost::extents[xsize][ysize]))
,param(other.getParams())
,bound(other.getBoundaries())
{
    (*data) = other.getData();
}


Lattice::~Lattice(){
    delete data;
    data = NULL;
}

//=========================== OPERATIONS ===========================

void Lattice::equilibriumIni()
{
    Cell tmp;
    array2D eqDis;

    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            tmp = (*data)[i][j];
            tmp.calcRho();
            double rho = tmp.getRho();
            Vector2D u = tmp.getU();
            eqDis = eqDistro(rho,u);
            tmp.setF(eqDis);
            (*data)[i][j] = tmp;
        }
    }
}

void Lattice::balance(double& mass, double& momentum)const
{
    Vector2D u;
    double rho;

    mass = 0;
    momentum = 0;

    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            (*data)[i][j].calcRho();
            rho = (*data)[i][j].getRho();
            u = (*data)[i][j].getU();

            mass += rho;
            momentum += rho * sqrt(u*u);
        }
    }
}

void Lattice::overallRho()
{
    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            (*data)[i][j].calcRho();
        }
    }
}

direction2D Lattice::directions(int x, int y)const
{
    direction2D dir;
    int tmp;
    for (int q=0; q<13; q++)
    {
        tmp = x + DIRECTION_2D[q].x;
        if (tmp<0) tmp += xsize;
        if (tmp>= xsize) tmp -= xsize;
        dir[q].x = tmp;

        tmp = y + DIRECTION_2D[q].y;
        if (tmp<0) tmp += ysize;
        if (tmp>= ysize) tmp -= ysize;
        dir[q].y = tmp;
    }
    return dir;
}

void Lattice::streamAll(int threads)
{
    field2D *newData = new field2D(boost::extents[xsize][ysize]);

    omp_set_num_threads (threads);
    //const int range = xsize * ysize;

    #pragma omp parallel
    {
        #pragma omp for collapse(2) schedule(dynamic, 100)
        for (int x=0; x<xsize; x++)
        {
            for (int y = 0; y < ysize; y++)
            {
                const direction2D dir = directions(x,y);
                Cell tmpCell = (*data)[x][y];

                if (tmpCell.getIsSolid() == false) streamAndBouncePull(tmpCell,dir);

                tmpCell.calcRho();
                (*newData)[x][y] = tmpCell;
            }
        }
    }
    delete data;
    data = newData;
}

bool Lattice::collideAll(int threads)
{
    bool success(true);
    field2D *newData = new field2D(boost::extents[xsize][ysize]);

    omp_set_num_threads (threads);

    const RelaxationPar2D relax = param.getRelaxation2D();

    const double omega = param.getOmega();

    const Matrix relaxation_matrix(relax,false,omega);
    const Matrix single_relax = INV_TRAFO_MATRIX * relaxation_matrix * TRAFO_MATRIX;

    #pragma omp parallel firstprivate(single_relax)
    {
        #pragma omp for collapse(2) schedule(dynamic, 100)
        for (int x=0; x<xsize; x++)
        {
            for (int y = 0; y < ysize; y++)
            {
                Cell tmpCell = (*data)[x][y];
                // Cell tmpCell = realData[x][y];

                if (tmpCell.getIsSolid() == false && isBoundary(x, y) == false)
                {

                    array2D  fTmp;
                    const array2D fCell = tmpCell.getF();

                    const double rho = tmpCell.getRho();
                    Vector2D u = tmpCell.getU();

                    const array2D diff = array_diff_2D(fCell, eqDistro(rho, u));
                    const array2D single_phase_col =  single_relax * diff;
                    
                    for (int q=0; q<9; q++)
                    {
                        fTmp[q] =  fCell[q] - single_phase_col[q];
                    } // end for
 
                    tmpCell.setF(fTmp);
                }
                (*newData)[x][y] = tmpCell;
            }
        }
    }
        
    if(success == true)
    {
        delete data;
        data = newData;
    }
    else
    {
        delete newData;
    }     

    return success;
}

void Lattice::evaluateBoundaries(int threads)
{
    omp_set_num_threads (threads);

    //#pragma omp parallel
    //{
        // north boundary 
        if(bound.north.getType() == pressure)
        {
            double rho = bound.north.getRho();

            int lowerX = 0;
            int upperX = xsize;

            if (bound.west.getType() != periodic) lowerX++;
            if (bound.east.getType() != periodic) upperX--;

            if(bound.west.getType() == pressure || bound.west.getType() == velocity)
            {
                (*data)[0][ysize-1] = cornerNorthWest((*data)[0][ysize-1], rho);
            }

            if(bound.east.getType() == pressure || bound.east.getType() == velocity)
            {
                (*data)[xsize-1][ysize-1] = cornerNorthEast((*data)[xsize-1][ysize-1], rho);
            }

            //#pragma omp for schedule(static,10) nowait
            for (int x=lowerX; x<upperX; x++)
            {
                (*data)[x][ysize-1] = boundaryNorthPres((*data)[x][ysize-1], rho);
            }
        } // end if north pressure

        if(bound.north.getType() == velocity)
        {
            double u_y = bound.north.getVelocity().y;

            int lowerX = 0;
            int upperX = xsize;

            if (bound.west.getType() != periodic) lowerX++;
            if (bound.east.getType() != periodic) upperX--;

            if(bound.west.getType() == pressure || bound.west.getType() == velocity)
            {
                (*data)[0][ysize-1] = cornerNorthWest((*data)[0][ysize-1], bound.north.getRho());
            }

            if(bound.east.getType() == pressure || bound.east.getType() == velocity)
            {
                (*data)[xsize-1][ysize-1] = cornerNorthEast((*data)[xsize-1][ysize-1], bound.north.getRho());
            }

            //#pragma omp for schedule(static,10) nowait
            for (int x=lowerX; x<upperX; x++)
            {
                (*data)[x][ysize-1] = boundaryNorthVelo((*data)[x][ysize-1], u_y);
            }
        } // end if north velocity

        // south boundary
        if(bound.south.getType() == pressure)
        {
            double rho = bound.south.getRho();
            int lowerX = 0;
            int upperX = xsize;
    
            if (bound.west.getType() != periodic) lowerX++;
            if (bound.east.getType() != periodic) upperX--;
    
            if(bound.west.getType() == pressure || bound.west.getType() == velocity)
            {
                (*data)[0][0] = cornerSouthWest((*data)[0][0], bound.south.getRho());
            }
    
            if(bound.east.getType() == pressure || bound.east.getType() == velocity)
            {
                (*data)[xsize-1][0] = cornerSouthEast((*data)[xsize-1][0], bound.south.getRho());
            }

            //#pragma omp for schedule(static,10) nowait
            for (int x=lowerX; x<upperX; x++)
            {
                (*data)[x][0] = boundarySouthPres((*data)[x][0],rho);
            }
        } // end if south pressure

        // south boundary
        if(bound.south.getType() == velocity)
        {
            double u_y = bound.south.getVelocity().y;

            int lowerX = 0;
            int upperX = xsize;
    
            if (bound.west.getType() != periodic) lowerX++;
            if (bound.east.getType() != periodic) upperX--;
    
            if(bound.west.getType() == pressure || bound.west.getType() == velocity)
            {
                (*data)[0][0] = cornerSouthWest((*data)[0][0], bound.south.getRho());
            }
    
            if(bound.east.getType() == pressure || bound.east.getType() == velocity)
            {
                (*data)[xsize-1][0] = cornerSouthEast((*data)[xsize-1][0], bound.south.getRho());
            }

            //#pragma omp for schedule(static,10) nowait
            for (int x=lowerX; x<upperX; x++)
            {
                (*data)[x][0] = boundarySouthVelo((*data)[x][0],u_y);
            }
        } // end if south velocity   

        // west boundary
        if(bound.west.getType() == pressure)
        {
            double rho = bound.west.getRho();
    
            int lowerY = 0;
            int upperY = ysize;
            if (bound.south.getType() != periodic) lowerY++;
            if (bound.north.getType() != periodic) upperY--;
        
            //#pragma omp for schedule(static,10) nowait
            for (int y=lowerY; y< upperY; y++)
            {
                (*data)[0][y] = boundarySouthVelo((*data)[0][y],rho);
            }
        }
        
        if(bound.west.getType() == velocity)
        {
            double u_x = bound.west.getVelocity().x;
    
            int lowerY = 0;
            int upperY = ysize;
            if (bound.south.getType() != periodic) lowerY++;
            if (bound.north.getType() != periodic) upperY--;
        
            //#pragma omp for schedule(static,10) nowait
            for (int y=lowerY; y< upperY; y++)
            {
                (*data)[0][y] = boundaryWestVelo((*data)[0][y], u_x);
            }
        }

        // east boundary
        if(bound.east.getType() == pressure)
        {
            double rho = bound.east.getRho();
    
            int lowerY = 0;
            int upperY = ysize;
    
            if (bound.south.getType() != periodic) lowerY++;
            if (bound.north.getType() != periodic) upperY--;

            //#pragma omp for schedule(static,10) nowait
            for (int y=lowerY; y< upperY; y++)
            {
                (*data)[xsize-1][y] = boundaryEastPres((*data)[xsize-1][y], rho);
            }
        }
        
        if(bound.east.getType() == velocity)
        {
            double u_x = bound.east.getVelocity().x;
    
            int lowerY = 0;
            int upperY = ysize;
    
            if (bound.south.getType() != periodic) lowerY++;
            if (bound.north.getType() != periodic) upperY--;

            //#pragma omp for schedule(static,10) nowait
            for (int y=lowerY; y< upperY; y++)
            {
                (*data)[xsize-1][y] = boundaryEastVelo((*data)[xsize-1][y], u_x) ;
            }
        }

    //}
}

void Lattice::closedBox()
{
    const Cell wall(0,true);

    for (int x=0; x<xsize; x++)
    {
        (*data)[x][0] = wall;
        (*data)[x][ysize-1] = wall;
    }
    for (int y=0; y<ysize; y++)
    {
        (*data)[0][y] = wall;
        (*data)[xsize-1][y] = wall;
    }

    for (int y=0; y<ysize; y++)
    {
        for (int x=0; x<xsize; x++)
        {
            (*data)[x][y].calcRho();
        }
    }
}

void Lattice::bottomWall()
{
    const Cell wall(0,true);

    for (int x=0; x<xsize; x++)
    {
        (*data)[x][0] = wall;
    }

    for (int y=0; y<ysize; y++)
    {
        for (int x=0; x<xsize; x++)
        {
            (*data)[x][y].calcRho();
        }
    }
}

void Lattice::genericWall(std::vector<double> x, std::vector<double> y, const Vector2D& u_w)
{
    if(x.size() == y.size())
    {
        Cell wall(0,true);
        wall.setSolidVelocity(u_w);
        
        for(unsigned int i = 0; i< x.size(); i++)
        {
            (*data)[x[i]][y[i]] = wall;
        }
    }
    else throw "vector size mismatch";    
}

void Lattice::lidDrivenCavity(const Vector2D& u_w)
{
    Cell wall(0,true);

    for (int y=0; y<ysize; y++)
    {
        (*data)[0][y] = wall;   // left wall
        (*data)[xsize-1][y] = wall; // right wall
    }

    for (int x=0; x<xsize; x++) 
    {
        (*data)[x][0] = wall; // bottom wall
    }

    wall.setSolidVelocity(u_w);
    for (int x=0; x<xsize; x++) 
    {
         (*data)[x][ysize-1] = wall;    // top wall
    }

    for (int y=0; y<ysize; y++)
    {
        for (int x=0; x<xsize; x++)
        {
            (*data)[x][y].calcRho();
        }
    }
}

void Lattice::shearWall(const Vector2D& u_w)
{
    Cell wall(0,true);

    for (int y=0; y<ysize; y++)
    {
        (*data)[xsize-1][y] = wall; // right wall
    }

    wall.setSolidVelocity(u_w);
    for (int y=0; y<ysize; y++)
    {
        (*data)[0][y] = wall;   // left wall
    }

    for (int y=0; y<ysize; y++)
    {
        for (int x=0; x<xsize; x++)
        {
            (*data)[x][y].calcRho();
        }
    }
}

void Lattice::setShearProfile(double gradient, double offset)
{
    const int range = xsize * ysize;
    const double m = gradient;
    const double n = offset - (m * ( xsize *  param.getDeltaX() ) );
    for (int index = 0;  index < range; index++)
        {
            int x,y;
            linearIndex(index,x,y);
            Cell tmpCell = (*data)[x][y];
            tmpCell.calcRho();
            Vector2D u;
            u = Vector2D(0,m*x + n);
            //const Vector2D u(0 , m*x + n);   
            const double rho = tmpCell.getRho();
            
            tmpCell.setF(eqDistro(rho,u));
            (*data)[x][y] = tmpCell;
        }

}

//=========================== ACCESSORS ===========================

const DimSet2D Lattice::getSize()const
{
    DimSet2D pony = {{xsize, ysize}}; 
    return pony;
}

void Lattice::setData(const field2D& ndata, int x, int y){
    data->resize(boost::extents[x][y]);
    *data = ndata;
    xsize = x;
    ysize = y;
}

void Lattice::setCell(int x, int y, const Cell& ncell)
{
    if (y >= 0 && y < ysize && x >= 0 && x < xsize) (*data)[x][y] = ncell;
}

void Lattice::setF(int x, int y, const array2D& nf)
{
    array2D f = (*data)[x][y].getF();
    f = nf;
    (*data)[x][y].setF(f);
}

void Lattice::setF(int x, int y, int pos, double value)
{
    array2D f = (*data)[x][y].getF();
    f[pos] = value;
    (*data)[x][y].setF(f);
}

void Lattice::setBoundaries(const Boundaries& newBound)
{
    bound = newBound;
    buildWalls();
}

//=========================== OPERATOR ===========================

Lattice& Lattice::operator=(const Lattice& other)
{
    this->setData(other.getData(), other.getSize()[0], other.getSize()[1]);
    this->setParams(other.getParams());
    this->setBoundaries(other.getBoundaries());

    return *this;
}

const bool Lattice::operator==(const Lattice& other)const
{
    bool exit = true;
    DimSet2D extent = other.getSize();
    if (xsize == extent[0] && ysize == extent[1])
    {
        ParamSet pOther = other.getParams();
        if (!(param == pOther)) exit = false;

        Boundaries bOther = other.getBoundaries();
        if (!(bound == bOther)) exit = false;

        field2D otherData = other.getData();
        for (int x = 0; x< xsize;x++)
        {
            for (int y = 0; y<ysize;y++)
            {
                if (!((*data)[x][y]==otherData[x][y]))
                {
                    exit = false;
                    break;
                }
            }
        }
    }
    else exit = false;
   
    return exit;
}

///////////////////////////// PRIVATE /////////////////////////////

//=========================== BOUNDARY TREATMENT ===========================

const bool Lattice::isBoundary(int x, int y)const
{
    return  (x == 0 && bound.west.getType()!= periodic) || (x == xsize-1 && bound.east.getType()!= periodic) || (y == ysize-1 && bound.north.getType()!= periodic) || (y == 0 && bound.south.getType()!= periodic);
}

void Lattice::buildWalls()
{
    Cell wall(0,true);

    if(bound.north.getType() == bounceback)
    {
        Vector2D u_wall = bound.north.getVelocity();
        wall.setSolidVelocity(u_wall);

        for (int x=0; x<xsize; x++)
        {
            (*data)[x][ysize-1] = wall;
        }
    }

    if(bound.south.getType() == bounceback)
    {
        Vector2D u_wall = bound.south.getVelocity();
        wall.setSolidVelocity(u_wall);

        for (int x=0; x<xsize; x++)
        {
            (*data)[x][0] = wall;
        }
    }

    if(bound.west.getType() == bounceback)
    {
        Vector2D u_wall = bound.west.getVelocity();
        wall.setSolidVelocity(u_wall);

        for (int y=0; y<ysize; y++)
        {
            (*data)[0][y] = wall;
        }
    }

    if(bound.east.getType() == bounceback)
    {
        Vector2D u_wall = bound.east.getVelocity();
        wall.setSolidVelocity(u_wall);

        for (int y=0; y<ysize; y++)
        {
            (*data)[xsize-1][y] = wall;
        }
    }    

}

const Cell Lattice::boundaryNorthPres(Cell tmp, double rho)const
{
    array2D f = tmp.getF();
    double uy;

    if(rho>0)
    {
        uy = -1.0 + (f[0] + f[1] + f[5] + 2* (f[2] + f[3] + f[4])) / rho;
        f[7] = f[3] - 2.0/3.0 * rho*uy;
        f[6] = - rho*uy/6.0 + f[2] + (f[1] - f[5])/2.0;
        f[8] = - rho*uy/6.0 + f[4] + (f[5] - f[1])/2.0;
    }
    else
    {
        uy = 0;
        f[7] = 0;
        f[6] = 0;
        f[8] = 0;
    }
    tmp.setF(f);
    return tmp;
}

const Cell Lattice::boundaryNorthVelo(Cell tmp, double uy)const
{
    //array2D f = tmp.getF();
    //double rho = (f[0] + f[1] + f[5] + 2* (f[2] + f[3] + f[4])) / (uy + 1);

    //f[7] = f[3] - 2.0/3.0 * rho*uy;
    //f[6] = -rho*uy/6.0 + f[2] + (f[1] - f[5])/2.0;
    //f[8] = -rho*uy/6.0 + f[4] + (f[5] - f[1])/2.0;

//    tmp.setF(f);
    tmp.setF(eqDistro(1.0, Vector2D(0,uy)));
    return tmp;
}

const Cell Lattice::boundarySouthPres(Cell tmp, double rho)const
{
    array2D f = tmp.getF();
    double uy;
    
    if(rho>0 )
    {
        uy = 1.0 - (f[0] + f[1] + f[5] + 2* (f[6] + f[7] + f[8])) / rho;
        f[3] = f[7] + 2.0/3.0 * rho*uy;
        f[2] = rho*uy/6.0 + f[6] + (f[5] - f[1])/2.0;
        f[4] = rho*uy/6.0 + f[8] + (f[1] - f[5])/2.0;
    }
    else
    {
        uy = 0;
        f[3] = 0;
        f[2] = 0;
        f[4] = 0;
    }
    tmp.setF(f);
    return tmp;
}

const Cell Lattice::boundarySouthVelo(Cell tmp, double uy)const
{
    // array2D f = tmp.getF(); 
    // double rho = (f[0] + f[1] + f[5] + 2* (f[6] + f[7] + f[8])) / (1 - uy);

    // f[3] = f[7] + 2.0/3.0 * rho*uy;
    // f[2] = rho*uy/6.0 + f[6] + (f[5] - f[1])/2.0;
    // f[4] = rho*uy/6.0 + f[8] + (f[1] - f[5])/2.0;
    
    // tmp.setF(f);
    tmp.setF(eqDistro(1.0, Vector2D(0,uy)));
    return tmp;
}

const Cell Lattice::boundaryWestPres(Cell tmp, double rho)const
{
    array2D f = tmp.getF();
    double ux;
    
    if(rho>0)
    {
        ux = 1.0 - (f[0] + f[3] + f[7] + 2* (f[4] + f[5] + f[6])) / rho;
        f[1] = f[5] + 2.0/3.0 * rho*ux;
        f[2] = rho*ux/6.0 + f[6] + (f[7] - f[3])/2.0;
        f[8] = rho*ux/6.0 + f[4] + (f[3] - f[7])/2.0;
    }
    else
    {
        ux = 0;
        f[1] = 0;
        f[2] = 0;
        f[8] = 0;
    }

    tmp.setF(f);
    return tmp;
}

const Cell Lattice::boundaryWestVelo(Cell tmp, double ux)const
{
    // array2D f = tmp.getF();
    // double rho = (f[0] + f[3] + f[7] + 2* (f[4] + f[5] + f[6])) / (1 - ux);
    
    // f[1] = f[5] + 1.5 * rho*ux;
    // f[2] = rho*ux/6.0 + f[6] + (f[7] - f[3]) * 0.5;
    // f[8] = rho*ux/6.0 + f[4] + (f[3] - f[7]) * 0.5;
    
    // tmp.setF(f);
    tmp.setF(eqDistro(1.0, Vector2D(ux,0)));
    return tmp;
}

const Cell Lattice::boundaryEastPres(Cell tmp, double rho)const
{
    array2D f = tmp.getF();
    double ux;
    
    if(rho>0)
    {
        ux = -1.0 + (f[0] + f[3] + f[7] + 2* (f[1] + f[2] + f[8])) / rho;
        f[5] = f[1] - 2.0/3.0 * rho*ux;
        f[6] = - rho*ux/6.0 + f[2] + (f[3] - f[7])/2.0;
        f[4] = - rho*ux/6.0 + f[8] + (f[7] - f[3])/2.0;
    }

    else
    {
        ux = 0;
        f[5] = 0;
        f[6] = 0;
        f[4] = 0;
    }

    tmp.setF(f);
    return tmp;
}

const Cell Lattice::boundaryEastVelo(Cell tmp, double ux)const
{
    // array2D f = tmp.getF();

    // double rho = (f[0] + f[3] + f[7] + 2* (f[1] + f[2] + f[8])) / (1 + ux);
    // f[5] = f[1] - 2.0/3.0 * rho*ux;
    // f[6] = - rho*ux/6.0 + f[2] + (f[3] - f[7])/2.0;
    // f[4] = - rho*ux/6.0 + f[8] + (f[7] - f[3])/2.0;
    // tmp.setF(f);

    tmp.setF(eqDistro(1.0, Vector2D(ux,0)));
    return tmp;
}

const Cell Lattice::cornerNorthWest(Cell tmp, double rho)const
{
    array2D f = tmp.getF();
    if(rho>0 )
    {
        f[1] = f[5];
        f[7] = f[3];
        f[8] = f[4];
        f[2] = 0.5 * (rho - f[0] )  - f[3] - f[4] - f[5];
        f[6] = f[2];
    }
    else
    {
        f[1] = 0;
        f[2] = 0;
        f[6] = 0;
        f[7] = 0;
        f[8] = 0;
    }

    tmp.setF(f);
    return tmp;    
}

const Cell Lattice::cornerNorthEast(Cell tmp, double rho)const
{
    array2D f = tmp.getF();
    if(rho>0 )
    {
        f[5] = f[1];
        f[6] = f[2];
        f[7] = f[3];
        f[4] = 0.5 * (rho - f[0] )  - f[1] - f[2] - f[3];
        f[8] = f[4];
    }
    else
    {
        f[4] = 0;
        f[5] = 0;
        f[6] = 0;
        f[7] = 0;
        f[8] = 0;
    }

    tmp.setF(f);
    return tmp;
}

const Cell Lattice::cornerSouthWest(Cell tmp, double rho)const
{
    array2D f = tmp.getF();
    if(rho>0)
    {
        f[1] = f[5];
        f[2] = f[6];
        f[3] = f[7];
        f[4] = 0.5 * (rho - f[0] )  - f[5] - f[6] - f[7];
        f[8] = f[4];
    }
    else
    {
        f[1] = 0;
        f[2] = 0;
        f[3] = 0;
        f[4] = 0;
        f[8] = 0;
    }

    tmp.setF(f);
    return tmp;   
}

const Cell Lattice::cornerSouthEast(Cell tmp, double rho)const
{
    array2D f = tmp.getF();
    if(rho>0 )
    {
        f[3] = f[7];
        f[4] = f[8];
        f[5] = f[1];
        f[2] = 0.5 * (rho - f[0] )  - f[1] - f[7] - f[8];
        f[6] = f[2];
    }
    else
    {
        f[3] = 0;
        f[4] = 0;
        f[5] = 0;
        f[2] = 0;
        f[6] = 0;
    }

    tmp.setF(f);
    return tmp;
}


//=========================== OPERATIONS ===========================

inline void Lattice::linearIndex(int index, int& x, int& y)const
{
    x = (index)%xsize;
    y = (index)/xsize;
}

void Lattice::streamAndBouncePull(Cell& tCell, const direction2D& dir)const
{
    const array2D f = tCell.getF();
    array2D ftmp;
    for (int color = 0; color<=1;color++)
    {
        ftmp[0] = (*data)[ dir[0].x ][ dir[0].y ].getF()[0];

        for (int i=1;i<9;i++)
        {
            const Cell neighbor = (*data)[ dir[PULL_INDEX_2D[i]].x ][ dir[PULL_INDEX_2D[i]].y ];
            if (neighbor.getIsSolid() == false) // if(neighbor not solid?) -> stream
            {
                ftmp[i] = neighbor.getF()[i];    
            }
            else  // else -> bounce back
            {
                const double rho = tCell.getRho();
                ftmp[i] = f[PULL_INDEX_2D[i]] - (2.0 * 3.0 * WEIGHTS_2D[i] * rho * (DIRECTION_2D[i] * neighbor.getU()) ) ;
            } 
        } // end for i
    } // end for color
    tCell.setF(ftmp);
}





///////////////////////////// C-STYLE /////////////////////////////

//=========================== OPERATIONS ===========================

const array2D eqDistro(const double rho, const Vector2D& u)
{
    array2D feq;
    double usqr = u*u;
   
    for (int i=0; i<9; i++)
    {
        double scal = u*DIRECTION_2D[i];
        feq[i] = rho * ( WEIGHTS_2D[i] * (1 + 3 * scal + 4.5 * (scal*scal) - 1.5 * usqr));
    }
    return feq;
}
