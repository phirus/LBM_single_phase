#ifndef DEFINITIONS_BASIC_H
#define DEFINITIONS_BASIC_H

#include<boost/multi_array.hpp>
#include"Vector2D.h"

/// contains custom typedefs
//=========================== TYPES ===========================
//=========================== TYPES ===========================

typedef boost::array<Vector2D,13> direction2D ; /// < collection of 13 direction vectors (D2Q13)
typedef boost::array<int,2> SizeSet2D; 
typedef boost::array<double,9> array2D;       /// < used to describe single distributions

typedef boost::array<int,2> DimSet2D;         /// < simple 2d vector for lattice positions and size

struct RelaxationPar2D
{
    double s_2,s_3,s_5;
    RelaxationPar2D(double s2=1, double s3=1, double s5=1):s_2(s2),s_3(s3),s_5(s5){};
};

//=========================== FUNCTIONS ===========================

/// functions handling basic operations on arrays
const array2D array_diff_2D(const array2D &one, const array2D &two);
const array2D array_add_2D(const array2D &one, const array2D &two);
const array2D array_times_2D(const array2D &foo, double factor);

#endif