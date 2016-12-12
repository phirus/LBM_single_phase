#ifndef TESTS_H
#define TESTS_H

#include"gtest/gtest.h"
#include"../src/BinaryIO2D.h"


using namespace std;

TEST(BasicIO,queryTest)
{
    double value;
    EXPECT_FALSE( input_query("existiertnicht","test",value) );
    EXPECT_FALSE(input_query("queryTest","noflag",value));
    EXPECT_TRUE(input_query("queryTest","test",value));
    EXPECT_DOUBLE_EQ(13.4, value); 
}

TEST(BinaryIO2D,output)
{
    const Lattice lattice(150,150);
    write_binary2D(lattice);
    Lattice vergleich;
    EXPECT_FALSE(read_binary2D(vergleich,"existiertnicht.txt"));
    EXPECT_TRUE(read_binary2D(vergleich));
    EXPECT_EQ(lattice,vergleich);
}

TEST(BinaryIO2D,restart)
{
    const Lattice lattice(150,150);

    Timetrack time(2e5,100,1000);
    time.timestep();

    const Preprocess newProcess = read_preprocess_file("preprocessFile");
    write_restart_file2D(lattice, newProcess,time);
    
    Lattice vergleichL;
    Preprocess vergleichP;
    Timetrack vergleichT;
    EXPECT_TRUE(read_restart_file2D(vergleichL,vergleichP,vergleichT));
    EXPECT_EQ(lattice, vergleichL);
    EXPECT_EQ(newProcess, vergleichP);
    EXPECT_EQ(time, vergleichT);
}

TEST(Boundaries,constructor)
{

    const boundary_type bt0 = periodic;
    const double density0 = 0;
    const Vector2D velocity0 = Vector2D();

    const boundary_type bt1 = pressure;
    const double density1 = 4;
    const Vector2D velocity1 = Vector2D(1,2);

    const BoundaryInformation bi0 = BoundaryInformation();
    const BoundaryInformation bi1 = BoundaryInformation(bt1, density1, velocity1);
    const BoundaryInformation bi2 = BoundaryInformation(1);
    const BoundaryInformation bi3 = BoundaryInformation(2);


    EXPECT_EQ(bi0.getType(),bt0);
    EXPECT_EQ(bi0.getRho(),density0);
    EXPECT_EQ(bi0.getVelocity(),velocity0);

    EXPECT_EQ(bi1.getType(),bt1);
    EXPECT_EQ(bi1.getRho(),density1);
    EXPECT_EQ(bi1.getVelocity(),velocity1);

    EXPECT_EQ(bi2.getType(),bounceback);
    EXPECT_EQ(bi3.getType(),pressure);
}

TEST(Boundaries,FileInput)
{
    const Boundaries b = read_boundaries_file("BoundaryInput");

    double density;
    Vector2D velocity;

    density = 1;
    velocity = Vector2D(3,4);
    EXPECT_EQ(b.north.getType(), bounceback);
    EXPECT_EQ(b.north.getRho(),density);
    EXPECT_EQ(b.north.getVelocity(),velocity);

    density = 6;
    velocity = Vector2D(8,9);
    EXPECT_EQ(b.south.getType(), periodic);
    EXPECT_EQ(b.south.getRho(),density);
    EXPECT_EQ(b.south.getVelocity(),velocity);

    density = 11;
    velocity = Vector2D(13,14);
    EXPECT_EQ(b.east.getType(), pressure);
    EXPECT_EQ(b.east.getRho(),density);
    EXPECT_EQ(b.east.getVelocity(),velocity);

    density = 16;
    velocity = Vector2D(18,19);
    EXPECT_EQ(b.west.getType(), periodic);
    EXPECT_EQ(b.west.getRho(),density);
    EXPECT_EQ(b.west.getVelocity(),velocity);
}

TEST(Cell,constructor0)
{
    const array2D f = {{1,0,0,0,0,0,0,0,0}};
    const Cell cell;
    EXPECT_EQ(f,cell.getF());
    EXPECT_FALSE(cell.getIsSolid());
}

TEST(Cell,constructorSolid)
{
    const Cell cell(0,true);
    EXPECT_EQ(true,cell.getIsSolid());
}

TEST(Cell,constructorIni)
{
    const array2D f1 = {{4,2,5.4,0,1,12,6,7,8}};
    const Cell cell(f1);
    EXPECT_EQ(f1,cell.getF());
     EXPECT_FALSE(cell.getIsSolid());
}

TEST(Cell,rho)
{
    const array2D f1 = {{10,20,30,40,50,60,70,80,90}};
    Cell cell(f1);
    EXPECT_EQ(0,cell.getRho());
    cell.calcRho();
    EXPECT_EQ(450,cell.getRho());
    cell.setIsSolid(true);
    cell.calcRho();
    EXPECT_EQ(0,cell.getRho());
}

TEST(Constants,W)
{
    EXPECT_DOUBLE_EQ(4,WEIGHTS_2D.at(0)*9);
    EXPECT_DOUBLE_EQ(1,WEIGHTS_2D.at(1)*9);
    EXPECT_DOUBLE_EQ(1,WEIGHTS_2D.at(2)*36);
}

TEST(Definitions,array2D_diff_add)
{
    const array2D one = {{100, 5.5, 1.4, -2, 0, 0, 1, 91.6, 45}};
    const array2D two = {{5,   0.6, 1.9, -2, 1.8, 100, 0.5, 91.6, 25}};
    const array2D vergleich = {{95, 4.9, -0.5, 0, -1.8, -100, 0.5, 0, 20}};

    EXPECT_EQ (vergleich, array_diff_2D(one, two));
    EXPECT_EQ (one, array_add_2D(vergleich, two));
    EXPECT_EQ (one, array_add_2D(two, vergleich));
}

TEST(Definitions2D,array2D_times)
{
    const array2D one = {{100, 5, 14, -2, 0, 0, 1, 916, 45}};
    const array2D vergleich = {{110, 5.5, 15.4, -2.2, 0, 0, 1.1, 1007.6, 49.5}};
    const array2D result = array_times_2D(one, 1.1);

    for(int i=0;i<9;i++)
    {
        EXPECT_DOUBLE_EQ (vergleich[i], result[i]);
    }    
}

TEST(Lattice,constructor)
{
    const Lattice lattice;
    EXPECT_EQ(10, lattice.getSize()[1]);
    EXPECT_EQ(10, lattice.getSize()[0]);
    EXPECT_EQ(1,lattice.getF(1,1)[0]);
}

TEST(Lattice,stream)
{
    // Phase 0: Distributions stream towards the center
    // Phase 1: Distribution stream from the center to the neighbouring sites
    Lattice lattice(3,3,0);

    const array2D f1 = {{0,1,0,0,0,0,0,0,0}};
    const array2D f2 = {{0,0,1,0,0,0,0,0,0}};
    const array2D f3 = {{0,0,0,1,0,0,0,0,0}};
    const array2D f4 = {{0,0,0,0,1,0,0,0,0}};
    const array2D f5 = {{0,0,0,0,0,1,0,0,0}};
    const array2D f6 = {{0,0,0,0,0,0,1,0,0}};
    const array2D f7 = {{0,0,0,0,0,0,0,1,0}};
    const array2D f8 = {{0,0,0,0,0,0,0,0,1}};

    lattice.setF(0,1,f1);
    lattice.setF(0,0,f2);
    lattice.setF(1,0,f3);
    lattice.setF(2,0,f4);
    lattice.setF(2,1,f5);
    lattice.setF(2,2,f6);
    lattice.setF(1,2,f7);
    lattice.setF(0,2,f8);
        
    const array2D fcenter = {{0,1,1,1,1,1,1,1,1}};

    lattice.streamAll();

    EXPECT_EQ(fcenter,lattice.getF(1,1));

    lattice.streamAll();

    EXPECT_EQ(f6,lattice.getF(0,0));
    EXPECT_EQ(f7,lattice.getF(1,0));
    EXPECT_EQ(f8,lattice.getF(2,0));
    EXPECT_EQ(f5,lattice.getF(0,1));
    EXPECT_EQ(f1,lattice.getF(2,1));
    EXPECT_EQ(f4,lattice.getF(0,2));
    EXPECT_EQ(f3,lattice.getF(1,2));
    EXPECT_EQ(f2,lattice.getF(2,2));
}


TEST(Lattice,bounceSmall)
{
    Lattice lattice(3,3,0);
    lattice.closedBox();

    const array2D fcenter = {{0,1,1,1,1,1,1,1,1}};
    const array2D fzero = {{0,0,0,0,0,0,0,0,0}};

    lattice.setF(1,1,fcenter);

    lattice.streamAll();

    EXPECT_EQ(fzero,lattice.getF(0,0));
    EXPECT_EQ(fzero,lattice.getF(0,1));
    EXPECT_EQ(fzero,lattice.getF(0,1));
    EXPECT_EQ(fzero,lattice.getF(1,0));
    EXPECT_EQ(fzero,lattice.getF(1,2));
    EXPECT_EQ(fzero,lattice.getF(2,0));
    EXPECT_EQ(fzero,lattice.getF(2,1));
    EXPECT_EQ(fzero,lattice.getF(2,2));

    EXPECT_EQ(fcenter,lattice.getF(1,1));

}

TEST(Lattice,bounceClosed)
{
    Lattice lattice(5,5,0);
    lattice.closedBox();
    lattice.setCell(2,2,Cell(0,true));
    const Cell tmp = lattice.getCell(2,2);
    EXPECT_EQ(true,tmp.getIsSolid());

    const array2D f1 = {{0,1,0,0,0,0,0,0,0}};
    const array2D f2 = {{0,0,1,0,0,0,0,0,0}};
    const array2D f3 = {{0,0,0,1,0,0,0,0,0}};
    const array2D f4 = {{0,0,0,0,1,0,0,0,0}};
    const array2D f5 = {{0,0,0,0,0,1,0,0,0}};
    const array2D f6 = {{0,0,0,0,0,0,1,0,0}};
    const array2D f7 = {{0,0,0,0,0,0,0,1,0}};
    const array2D f8 = {{0,0,0,0,0,0,0,0,1}};

    lattice.setF(1,1,f2);
    lattice.setF(2,1,f3);
    lattice.setF(3,1,f4);
    lattice.setF(1,2,f1);
    lattice.setF(3,2,f5);
    lattice.setF(1,3,f8);
    lattice.setF(2,3,f7);
    lattice.setF(3,3,f6);

    lattice.streamAll();

    EXPECT_EQ(f6,lattice.getF(1,1));
    EXPECT_EQ(f7,lattice.getF(2,1));
    EXPECT_EQ(f8,lattice.getF(3,1));
    EXPECT_EQ(f5,lattice.getF(1,2));
    EXPECT_EQ(f1,lattice.getF(3,2));
    EXPECT_EQ(f4,lattice.getF(1,3));
    EXPECT_EQ(f3,lattice.getF(2,3));
    EXPECT_EQ(f2,lattice.getF(3,3));

}

TEST(Lattice,bounceClosed2)
{
    Lattice lattice(5,5,0);
    lattice.closedBox();
    const array2D fcenter = {{0,1,1,1,1,1,1,1,1}};
    lattice.setF(2,2,fcenter);
    lattice.streamAll();
    lattice.streamAll();
    lattice.streamAll();

    EXPECT_EQ(fcenter,lattice.getF(2,2));
}

TEST(Lattice,periodic)
{
    Lattice lattice(5,5,0);
    const array2D fcenter = {{0,1,1,1,1,1,1,1,1}};
    lattice.setF(2,2,fcenter);

    // 5 mal
    lattice.streamAll();
    lattice.streamAll();
    lattice.streamAll();
    lattice.streamAll();
    lattice.streamAll();

    EXPECT_EQ(fcenter,lattice.getF(2,2));
}

TEST(Lattice,streamRho)
{
    Cell tmp;

    Lattice lattice(5,5,0);
    lattice.closedBox();

    const array2D fcenter = {{0,1,1,1,1,1,1,1,1}};
    lattice.setF(2,2,fcenter);
    lattice.streamAll();

    tmp = lattice.getCell(2,2);
    EXPECT_EQ(0,tmp.getRho());
    tmp = lattice.getCell(1,2);
    EXPECT_EQ(1,tmp.getRho());
    tmp = lattice.getCell(3,2);
    EXPECT_EQ(1,tmp.getRho());

    lattice.streamAll();
    tmp = lattice.getCell(1,2);
    EXPECT_EQ(1,tmp.getRho());
    tmp = lattice.getCell(3,2);
    EXPECT_EQ(1,tmp.getRho());
}

TEST(Lattice,collideSingle)
{
    Lattice lattice(5,5,1);
    lattice.closedBox();
    
    lattice.collideAll();
    Cell cell = lattice.getCell(2,2);
    cell.calcRho();
    const double rho = cell.getRho();
    const Vector2D u = cell.getU();
    const double usqr = u*u;

    EXPECT_DOUBLE_EQ(1.0,rho);
    EXPECT_NEAR(0,usqr,1e-10);
}

TEST(Lattice, directions)
{
    const Lattice lattice(5,5);
    const direction2D dir = lattice.directions(0,0);

    const boost::array<int,13> x = {{0,1,1,0,4,4,4,0,1,2,0,3,0}};
    const boost::array<int,13> y = {{0,0,1,1,1,0,4,4,4,0,2,0,3}};
    for(int q=0; q<13;q++){
    EXPECT_EQ(y[q],dir[q].y);
    EXPECT_EQ(x[q],dir[q].x);
    }
}

TEST(Lattice,collisionBalanceAll)
{
    Lattice lattice(5,5,1);
    lattice.closedBox();
    double mass, momentum;

    lattice.balance(mass, momentum);
    EXPECT_DOUBLE_EQ(9,mass);
    EXPECT_DOUBLE_EQ(0,momentum);

    lattice.collideAll();
    lattice.streamAll();
    lattice.balance(mass, momentum);
    EXPECT_DOUBLE_EQ(9,mass);
    EXPECT_NEAR(0,momentum,1e-10);

    lattice.collideAll();
    lattice.streamAll();

    lattice.balance(mass, momentum);
    EXPECT_DOUBLE_EQ(9,mass);
    EXPECT_NEAR(0,momentum,1e-10);
}

TEST(Lattice, copy_constr)
{
    const Lattice lBig(100,20);
    const Lattice tmp(lBig);
    EXPECT_EQ(lBig,tmp);
}

TEST(Lattice, assign)
{ 
    const Lattice lBig(100,20);
    const Lattice tmp= lBig;
    EXPECT_EQ(lBig,tmp);
}


TEST(Matrix,trafo)
{
    const array2D verteilung = {{1,2,3,4,5,6,7,8,9}};
    const array2D vergleich = {{45,24,-12,-4,8,-12,0,-4,-4}};

    const array2D trafo = TRAFO_MATRIX * verteilung;

    EXPECT_EQ(vergleich, trafo);
}

TEST(Matrix,backtrafo)
{
    const array2D vergleich = {{1,2,3,4,5,6,7,8,9}};
    const array2D verteilung = {{45,24,-12,-4,8,-12,0,-4,-4}};

    const array2D backtrafo = INV_TRAFO_MATRIX * verteilung;

    for(int i = 0; i<9;i++)
    {
        EXPECT_DOUBLE_EQ(vergleich[i],backtrafo[i]);
    }
}

TEST(Matrix,multiply)
{
    const Matrix S(RelaxationPar2D(1,10,100),false);
    const array2D f = {{1,2,3,4,5,6,7,8,9}};
    const array2D vergleich = {{ 1, 2, 30, 4, 500, 6, 700, 8, 9}};

    array2D test1 = S*f;
    array2D test2 = S.diagMult(f);

    for(int i = 0; i<9;i++)
    {
        EXPECT_DOUBLE_EQ(vergleich[i],test1[i])<<"i = "<<i ;
        EXPECT_DOUBLE_EQ(vergleich[i],test2[i])<<"i = "<<i ;
    }
}

TEST(Matrix,multiply_matrix)
{
    const Matrix Identity(true);
    
    const Matrix test1 = TRAFO_MATRIX * INV_TRAFO_MATRIX;
    const Matrix test2 = INV_TRAFO_MATRIX * TRAFO_MATRIX;

    const boost::multi_array<double,2> id = Identity.getData();
    const boost::multi_array<double,2> t1 = test1.getData();
    const boost::multi_array<double,2> t2 = test2.getData();

    for (int i = 0;i <9;i++)
    {
        for (int j=0;j<9;j++)
        {
            EXPECT_NEAR(id[i][j],t1[i][j],1e-15);
            EXPECT_NEAR(id[i][j],t2[i][j],1e-15);
        }
    }
}

TEST(Matrix,multiply_matrix_left_right_alright)
{
    const Matrix zero(false);

    boost::multi_array<double,2> a = zero.getData();
    boost::multi_array<double,2> b = zero.getData();
    boost::multi_array<double,2> v1 = zero.getData();
    boost::multi_array<double,2> v2 = zero.getData();

    a[0][0] = 1;
    a[0][1] = 2;
    a[1][0] = 3;
    a[1][1] = 4;
    const Matrix A = Matrix(a);

    b[0][0] = 5;
    b[0][1] = 6;
    b[1][0] = 7;
    b[1][1] = 8;
    const Matrix B = Matrix(b);

    v1[0][0] = 19;
    v1[0][1] = 22;
    v1[1][0] = 43;
    v1[1][1] = 50;

    v2[0][0] = 23;
    v2[0][1] = 34;
    v2[1][0] = 31;
    v2[1][1] = 46;

    const Matrix test1 = A * B;
    const Matrix test2 = B * A;
    const boost::multi_array<double,2> t1 = test1.getData();
    const boost::multi_array<double,2> t2 = test2.getData();

    for (int i = 0;i <9;i++)
    {
        for (int j=0;j<9;j++)
        {
            EXPECT_NEAR(v1[i][j],t1[i][j],1e-15);
            EXPECT_NEAR(v2[i][j],t2[i][j],1e-15);
        }
    }
}

TEST(Matrix,plus_times){
    const Matrix Identity = Matrix(1.0);
    
    EXPECT_EQ(Identity+Identity, Identity*2);

    EXPECT_EQ(TRAFO_MATRIX+TRAFO_MATRIX+TRAFO_MATRIX, TRAFO_MATRIX*3);
    EXPECT_EQ(TRAFO_MATRIX+TRAFO_MATRIX, (TRAFO_MATRIX*3)-TRAFO_MATRIX);
}

TEST(ParamSet,IO)
{
    const RelaxationPar2D rel = RelaxationPar2D(0.8,1.2,1.2);
    const ParamSet params(1.1, 1.1, 1e-4, 1e-3,rel);   
    EXPECT_NO_THROW(write_param_log(params));
    EXPECT_NO_THROW(write_param_log_csv(params));

    const ParamSet read_in = read_paramset_file();
    EXPECT_EQ(params, read_in);
}

TEST(Preprocess,constr)
{
    const Preprocess newProcess;

    // test the given parameters (default values)
    EXPECT_DOUBLE_EQ(10,newProcess.getReynoldsMax());
    EXPECT_EQ(30,newProcess.getResolution());
    EXPECT_DOUBLE_EQ(1,newProcess.getRho());
    EXPECT_DOUBLE_EQ(2,newProcess.getMuRatio());

    // test the deduced parameters
    EXPECT_DOUBLE_EQ(1,newProcess.getSpacestep());
    EXPECT_DOUBLE_EQ(1,newProcess.getTimestep());

    EXPECT_DOUBLE_EQ(1.0196152422706632,newProcess.getTau());
    EXPECT_DOUBLE_EQ(0.5773502691896258,newProcess.getSoundspeed());
    EXPECT_DOUBLE_EQ(0.17320508075688779,newProcess.getNu());
 
    EXPECT_EQ(50,newProcess.getXCells());
    EXPECT_EQ(50,newProcess.getYCells());
    EXPECT_EQ(50,newProcess.getZCells());
}

TEST(Preprocess,FileInput)
{
    const Preprocess newProcess = read_preprocess_file("preprocessFile");

    // test the given parameters (default values)
    EXPECT_DOUBLE_EQ(15,newProcess.getReynoldsMax());
    EXPECT_EQ(35,newProcess.getResolution());
    EXPECT_DOUBLE_EQ(1.1,newProcess.getRho());
    EXPECT_DOUBLE_EQ(1.8,newProcess.getMuRatio());
    EXPECT_TRUE(newProcess.getIsShearFlow());
    EXPECT_DOUBLE_EQ(0.4, newProcess.getShearRate());

    // test the stored parameters
    EXPECT_DOUBLE_EQ(1,newProcess.getSpacestep());
    EXPECT_DOUBLE_EQ(1,newProcess.getTimestep());
    EXPECT_EQ(70,newProcess.getXCells());
    EXPECT_EQ(80,newProcess.getYCells());
    EXPECT_EQ(90,newProcess.getZCells());

    // test the deduced parameters
    EXPECT_DOUBLE_EQ(0.9041451884327381,newProcess.getTau());
    EXPECT_DOUBLE_EQ(0.5773502691896258,newProcess.getSoundspeed());
    EXPECT_DOUBLE_EQ(0.13471506281091272,newProcess.getNu());
}

TEST(Preprocess,ParameterCheck)
{
    const Preprocess newProcess = read_preprocess_file("realParameters");
    const ParamSet params = newProcess.getParamSet();

    // test the given parameters (default values)
    EXPECT_DOUBLE_EQ(75,newProcess.getReynoldsMax());
    
    EXPECT_DOUBLE_EQ(1.6558400984359005, params.getOmega());
    EXPECT_DOUBLE_EQ(1,params.getRho());
}

TEST(timetrack,FileInput)
{
    const Timetrack newTimetrack = read_timetrack_file("preprocessFile");
     // test the given parameters 
    EXPECT_DOUBLE_EQ(1e6,newTimetrack.getMaxCount());
    EXPECT_DOUBLE_EQ(1000,newTimetrack.getOutputInt());
    EXPECT_DOUBLE_EQ(10000,newTimetrack.getRestartInt());
}

TEST(Vector2D,scalar)
{
    const Vector2D v0, v1(1,2), v2(3,4);
    EXPECT_DOUBLE_EQ(0, v0*v1);
    EXPECT_DOUBLE_EQ(0, v2*v0);
    EXPECT_DOUBLE_EQ(11, v1*v2);
    EXPECT_DOUBLE_EQ(11, v2*v1);
}

#endif