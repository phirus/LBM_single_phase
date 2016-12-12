#include "Constants.h"

const boost::multi_array<double,2> define_trafo_matrix()
{
    boost::multi_array<double,2> matrix(boost::extents[9][9]);

    matrix[0][0] = 1;
    matrix[0][1] = 1;
    matrix[0][2] = 1;
    matrix[0][3] = 1;
    matrix[0][4] = 1;
    matrix[0][5] = 1;
    matrix[0][6] = 1;
    matrix[0][7] = 1;
    matrix[0][8] = 1;
    
    matrix[1][0] = -4;
    matrix[1][1] = -1;
    matrix[1][2] =  2;
    matrix[1][3] = -1;
    matrix[1][4] =  2;
    matrix[1][5] = -1;
    matrix[1][6] =  2;
    matrix[1][7] = -1;
    matrix[1][8] =  2;

    matrix[2][0] =  4;
    matrix[2][1] = -2;
    matrix[2][2] =  1;
    matrix[2][3] = -2;
    matrix[2][4] =  1;
    matrix[2][5] = -2;
    matrix[2][6] =  1;
    matrix[2][7] = -2;
    matrix[2][8] =  1;

    matrix[3][0] =  0;
    matrix[3][1] =  1;
    matrix[3][2] =  1;
    matrix[3][3] =  0;
    matrix[3][4] = -1;
    matrix[3][5] = -1;
    matrix[3][6] = -1;
    matrix[3][7] =  0;
    matrix[3][8] =  1;

    matrix[4][0] =  0;
    matrix[4][1] = -2;
    matrix[4][2] =  1;
    matrix[4][3] =  0;
    matrix[4][4] = -1;
    matrix[4][5] =  2;
    matrix[4][6] = -1;
    matrix[4][7] =  0;
    matrix[4][8] =  1;

    matrix[5][0] =  0;
    matrix[5][1] =  0;
    matrix[5][2] =  1;
    matrix[5][3] =  1;
    matrix[5][4] =  1;
    matrix[5][5] =  0;
    matrix[5][6] = -1;
    matrix[5][7] = -1;
    matrix[5][8] = -1;

    matrix[6][0] =  0;
    matrix[6][1] =  0;
    matrix[6][2] =  1;
    matrix[6][3] = -2;
    matrix[6][4] =  1;
    matrix[6][5] =  0;
    matrix[6][6] = -1;
    matrix[6][7] =  2;
    matrix[6][8] = -1;

    matrix[7][0] =  0;
    matrix[7][1] =  1;
    matrix[7][2] =  0;
    matrix[7][3] = -1;
    matrix[7][4] =  0;
    matrix[7][5] =  1;
    matrix[7][6] =  0;
    matrix[7][7] = -1;
    matrix[7][8] =  0;

    matrix[8][0] =  0;
    matrix[8][1] =  0;
    matrix[8][2] =  1;
    matrix[8][3] =  0;
    matrix[8][4] = -1;
    matrix[8][5] =  0;
    matrix[8][6] =  1;
    matrix[8][7] =  0;
    matrix[8][8] = -1;

    return matrix;
}


const boost::multi_array<double,2> define_inverse_trafo_matrix()
{
    boost::multi_array<double,2> matrix(boost::extents[9][9]);

    double t = 1;
    t/=36;

    matrix[0][0] =  4*t;
    matrix[0][1] = -4*t;
    matrix[0][2] =  4*t;
    matrix[0][3] =  0;
    matrix[0][4] =  0;
    matrix[0][5] =  0;
    matrix[0][6] =  0;
    matrix[0][7] =  0;
    matrix[0][8] =  0;

    matrix[1][0] =  4*t;
    matrix[1][1] = -1*t;
    matrix[1][2] = -2*t;
    matrix[1][3] =  6*t;
    matrix[1][4] = -6*t;
    matrix[1][5] =  0;
    matrix[1][6] =  0;
    matrix[1][7] =  9*t;
    matrix[1][8] =  0;

    matrix[2][0] =  4*t;
    matrix[2][1] =  2*t;
    matrix[2][2] =  1*t;
    matrix[2][3] =  6*t;
    matrix[2][4] =  3*t;
    matrix[2][5] =  6*t;
    matrix[2][6] =  3*t;
    matrix[2][7] =  0;
    matrix[2][8] =  9*t;

    matrix[3][0] =  4*t;
    matrix[3][1] = -1*t;
    matrix[3][2] = -2*t;
    matrix[3][3] =  0;
    matrix[3][4] =  0;
    matrix[3][5] =  6*t;
    matrix[3][6] = -6*t;
    matrix[3][7] = -9*t;
    matrix[3][8] =  0;

    matrix[4][0] =  4*t;
    matrix[4][1] =  2*t;
    matrix[4][2] =  1*t;
    matrix[4][3] = -6*t;
    matrix[4][4] = -3*t;
    matrix[4][5] =  6*t;
    matrix[4][6] =  3*t;
    matrix[4][7] =  0;
    matrix[4][8] = -9*t;

    matrix[5][0] =  4*t;
    matrix[5][1] = -1*t;
    matrix[5][2] = -2*t;
    matrix[5][3] = -6*t;
    matrix[5][4] =  6*t;
    matrix[5][5] =  0;
    matrix[5][6] =  0;
    matrix[5][7] =  9*t;
    matrix[5][8] =  0;

    matrix[6][0] =  4*t;
    matrix[6][1] =  2*t;
    matrix[6][2] =  1*t;
    matrix[6][3] = -6*t;
    matrix[6][4] = -3*t;
    matrix[6][5] = -6*t;
    matrix[6][6] = -3*t;
    matrix[6][7] =  0;
    matrix[6][8] =  9*t;

    matrix[7][0] =  4*t;
    matrix[7][1] = -1*t;
    matrix[7][2] = -2*t;
    matrix[7][3] =  0;
    matrix[7][4] =  0;
    matrix[7][5] = -6*t;
    matrix[7][6] =  6*t;
    matrix[7][7] = -9*t;
    matrix[7][8] =  0;

    matrix[8][0] =  4*t;
    matrix[8][1] =  2*t;
    matrix[8][2] =  1*t;
    matrix[8][3] =  6*t;
    matrix[8][4] =  3*t;
    matrix[8][5] = -6*t;
    matrix[8][6] = -3*t;
    matrix[8][7] =  0;
    matrix[8][8] = -9*t;

    return matrix;
}
