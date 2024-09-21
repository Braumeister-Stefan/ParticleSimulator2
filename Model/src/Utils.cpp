

//helper functions that dont belong elsewhere
#include "../include/Utils.h"

using namespace std;


int intensity_to_rgb(double r, double g, double b) {
    //this function will convert the rgb values to a hex code

    //1. convert the intensity values to 255 base

    int r255 = r * 255;
    int g255 = g * 255;
    int b255 = b * 255;

    //2. combine the rgb  

    int rgb = (r255 << 16) | (g255 << 8) | (b255);

    return rgb;

}

// particle.r = get<float>(rgb_parameters[0]);
//                 particle.g = get<float>(rgb_parameters[1]);
//                 particle.b = get<float>(rgb_parameters[2]);

//                 int rgb = ((int)(get<float>(rgb_parameters[0]) * 255) << 16) | ((int)(get<float>(rgb_parameters[1]) * 255) << 8) | (int)(get<float>(rgb_parameters[2]) * 255);