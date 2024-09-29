//This file contains geometric and other mathematical structures

#include "Particles.h"

#ifndef MATH_UTILS_H
#define MATH_UTILS_H

//Statistics structures

//Geometry structures

struct Vector2D {
    double x;
    double y;

    // Addition operator
    Vector2D operator+(const Vector2D& other) const {
        return { x + other.x, y + other.y };
    }

    // Subtraction operator
    Vector2D operator-(const Vector2D& other) const {
        return { x - other.x, y - other.y };
    }

    // Multiplication by scalar (Vector * Scalar)
    Vector2D operator*(double scalar) const {
        return { x * scalar, y * scalar };
    }

    // Multiplication by scalar (Scalar * Vector)
    friend Vector2D operator*(double scalar, const Vector2D& vector) {
        return { vector.x * scalar, vector.y * scalar };
    }

    // Division by scalar
    Vector2D operator/(double scalar) const {
        return { x / scalar, y / scalar };
    }
};

// Dot product function
inline double dot(const Vector2D& v1, const Vector2D& v2) {
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

double calc_distance(Vector2D point1, Vector2D point2);

Vector2D calc_midpoint(Vector2D point1, Vector2D point2);

double calc_magnitude(Vector2D vector);

Vector2D find_intersection(Line line1, Circle circle, Vector2D pos);

Vector2D sample_in_circle(Vector2D center, double radius);


















#endif // MATH_UTILS_H