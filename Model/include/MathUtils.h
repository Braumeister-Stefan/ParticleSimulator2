//This file contains geometric and other mathematical structures

//Internal libraries
#include "Particles.h"

#include <boost/multiprecision/cpp_dec_float.hpp>


//namespaces
using namespace boost::multiprecision;
using high_prec = boost::multiprecision::cpp_dec_float_50;


//Boosterized structures

struct Vector2D {
    high_prec x;
    high_prec y;

};

struct Vector2D_double {
    double x;
    double y;

};

// Define the multiplication operator for high_prec and Vector2D
inline Vector2D operator*(const high_prec& scalar, const Vector2D& vec) {
    return { vec.x * scalar, vec.y * scalar };
}

//define the Vector & Scalar operators

inline Vector2D operator+(const high_prec& scalar, const Vector2D& vec) {
    return { vec.x + scalar, vec.y + scalar };}
inline Vector2D operator-(const high_prec& scalar, const Vector2D& vec) {
    return { vec.x - scalar, vec.y - scalar };}
inline Vector2D operator/(const Vector2D& vec, const high_prec& scalar) {
    return { vec.x / scalar, vec.y / scalar };}
//define the Vector & Vector operators
inline Vector2D operator+(const Vector2D& vec1, const Vector2D& vec2) {
    return { vec1.x + vec2.x, vec1.y + vec2.y };}

inline Vector2D operator-(const Vector2D& vec1, const Vector2D& vec2) {
    return { vec1.x - vec2.x, vec1.y - vec2.y };}

inline Vector2D operator*(const Vector2D& vec1, const Vector2D& vec2) {
    return { vec1.x * vec2.x, vec1.y * vec2.y };}

inline Vector2D operator/(const Vector2D& vec1, const Vector2D& vec2) {
    return { vec1.x / vec2.x, vec1.y / vec2.y };}




// Dot product function
inline high_prec dot(const Vector2D& v1, const Vector2D& v2) {
    return v1.x * v2.x + v1.y * v2.y;
}



struct Circle {
    Vector2D center;
    double radius;
};

struct Line {
    Vector2D start;
    Vector2D end;
};

high_prec calc_distance(Vector2D point1, Vector2D point2);
Vector2D calc_midpoint(Vector2D point1, Vector2D point2);
high_prec calc_magnitude(Vector2D vector);
Vector2D find_intersection(Line line1, Circle circle, Vector2D pos);
Vector2D sample_in_circle(Vector2D center, double radius);


