#ifndef CONSTANTS_H
#define CONSTANTS_H

#include"Matrix.h"
#include"Definitions.h"

//=========================== CONSTANTS ===========================

/// Pi
const double PI = 3.14159265358979323846264338327950288419716939937510;

/// arbitrary  definition
const double MACH_MAX = 0.1; // maximal erlaubte Mach-Zahl


/// D2Q9 weights
const array2D WEIGHTS_2D = {{0.4444444444444444, 0.1111111111111111, 0.02777777777777778, 0.1111111111111111, 0.02777777777777778, 0.1111111111111111, 0.02777777777777778, 0.1111111111111111, 0.02777777777777778}};

/// D2Q13 directions
const Vector2D e0(0,0),e1(1,0),e2(1,1),e3(0,1),e4(-1,1),e5(-1,0),e6(-1,-1),e7(0,-1),e8(1,-1),e9(2,0),e10(0,2),e11(-2,0),e12(0,-2);
const boost::array<Vector2D,13> DIRECTION_2D = {{e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12}};

/// D2Q9 streaming indices
const boost::array<int,9> PULL_INDEX_2D = {{0,5,6,7,8,1,2,3,4}};

/// Transformation-Matrix
const boost::multi_array<double,2> define_trafo_matrix();
const boost::multi_array<double,2> define_inverse_trafo_matrix();

const Matrix TRAFO_MATRIX(define_trafo_matrix());
const Matrix INV_TRAFO_MATRIX(define_inverse_trafo_matrix());

#endif