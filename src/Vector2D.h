/** TheVector class handles typical 2D vectors
* provides functionfor absolute value and angles between two vectors 
*/

#ifndef VECTOR2D_H
#define VECTOR2D_H

#include<cmath>

class Vector2D
{
    public:
    double x,y;

    /// Lifecycle
    Vector2D(double xn = 0, double yn = 0):x(xn),y(yn){};

    /// operators
    inline const Vector2D operator+(const Vector2D& other)const{return Vector2D(x+other.x, y+other.y);};  /// < addition-operator
    inline const Vector2D operator-(const Vector2D& other)const{return Vector2D(x-other.x, y-other.y);};  /// < subtraction-operator
    inline const double operator*(const Vector2D&other)const{return x*other.x + y*other.y;};   /// < scalar product
    inline const Vector2D operator*(double c)const{return Vector2D(x*c,y*c);};                   /// < multiplication with a number
    inline const bool operator==(const Vector2D& other)const{return (x == other.x && y == other.y);}; 

    /// operations
    inline const double Abs()const{return sqrt(x*x+y*y);};   /// < absolute value
    inline const double Sum()const{return x+y;};             /// < sum over all components
};

#endif
