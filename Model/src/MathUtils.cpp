//implements the functions declared in MathUtils.h

//Standard libraries
#include <iostream>
#include <cmath>
#include <vector>

//Internal libraries
#include "../include/MathUtils.h"
#include "../include/Particles.h"

//defining global constants
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif




//namespaces
using namespace std;

//Statistics functions

Vector2D sample_in_circle(Vector2D center, double radius) {
    //This function will sample a random point within a circle of a given radius
    
    //sample from a uniform distribution a random angle between 0 and 2*pi and a random radius between 0 and the circle radius
    double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
    double r = radius * sqrt(rand() / (double)RAND_MAX);

    //print random number
    cout << "Random number in sample_in_circle(): " << rand() << endl;

    //initialize the point struct

    Vector2D point;

    //calculate the x and y coordinates of the sampled point
    point.x = center.x + r * cos(angle);
    point.y = center.y + r * sin(angle);
    
    return point;
}


//Geometry functions - Boosterized

// Function to calculate the distance between two points
high_prec calc_distance(Vector2D point1, Vector2D point2) {
    return hypot(point1.x - point2.x, point1.y - point2.y);
}

// Function to calculate the midpoint between two points
Vector2D calc_midpoint(Vector2D point1, Vector2D point2) {
    Vector2D midpoint;
    midpoint.x = (point1.x + point2.x) / 2.0;
    midpoint.y = (point1.y + point2.y) / 2.0;
    return midpoint;
}

// Function to calculate the magnitude of a vector
high_prec calc_magnitude(Vector2D vector) {
    return hypot(vector.x, vector.y);
}

// Function to find the intersection of two lines
Vector2D find_intersection(Line line1, Circle circle, Vector2D particle_coords) {
    // 1. Calculate the slope of the line
    high_prec slope = (line1.end.y - line1.start.y) / (line1.end.x - line1.start.x);

    // 2. Calculate the y-intercept of the line
    high_prec y_intercept = line1.start.y - slope * line1.start.x;

    // 3. Calculate the distance between the center of the circle and the line
    high_prec distance = abs(slope * circle.center.x - circle.center.y + y_intercept) / sqrt(pow(slope, 2) + 1);

    // 4. Calculate the intersection points
    high_prec x1 = (circle.center.x + slope * (circle.center.y - y_intercept) + sqrt(pow(circle.radius, 2) * pow(slope, 2) + pow(circle.radius, 2))) / (pow(slope, 2) + 1);
    high_prec x2 = (circle.center.x + slope * (circle.center.y - y_intercept) - sqrt(pow(circle.radius, 2) * pow(slope, 2) + pow(circle.radius, 2))) / (pow(slope, 2) + 1);

    high_prec y1 = slope * x1 + y_intercept;
    high_prec y2 = slope * x2 + y_intercept;

    // 5. Select the intersection point closest to the particle
    Vector2D intersection1;
    Vector2D intersection2;

    intersection1.x = x1;
    intersection1.y = y1;

    intersection2.x = x2;
    intersection2.y = y2;

    if (calc_distance(intersection1, particle_coords) < calc_distance(intersection2, particle_coords)) {
        return intersection1;
    } else {
        return intersection2;
    }
}