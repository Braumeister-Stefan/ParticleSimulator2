#include "../include/ParticlePlotter.h"
#include "../include/InitStructs.h"
#include <memory> // Add this line to include the necessary header for shared_ptr
#include <iostream>
#include <windows.h>
#include <fstream>
#include <sstream>
#include <stdio.h>

extern "C" FILE *popen(const char *command, const char *mode); //these 2 lines are added as popen and pclose were not recognized. Note that these are C functions, but worked without issue on ParticleSimulator1...likely need to change something about the compiler settings
extern "C" void pclose(FILE *pipe);


//namespaces
using namespace std;

//global variables
FILE *gnuplotPipe = nullptr;

// Constructor
Plotter::Plotter() {
    cout << "Plotting Engine initialized." << endl;
}

// Destructor
Plotter::~Plotter() {
    cout << "Plotting Engine destroyed." << endl;
}

// Define the plot_run method
void Plotter::plot_run(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states, shared_ptr<test_metrics> metrics) { 
    
    //1. convert r,g,b to hex
    for (int i = 0; i < particle_states->snaps.size(); i++) {
        particle_states->snaps[i] = convert_intensity_to_rgb(particle_states->snaps[i]);
    }
    cout << "RGB values have been calculated." << endl;

    //2. Initialize the plot using GNUplot

    init_GNU(scenario);

    //3. Populate the plot with the first snapshot


    //4. Confirm that the user is ready to start the simulation
    cout << "ready to plot scenario:" << scenario->name << endl;
    cout << "Press enter to continue." << endl;
    cin.ignore();
    cin.get();


    //5. Loop through the snapshots and plot each one on the plot
    for (int i = 0; i < particle_states->snaps.size(); i++) {
        cout << "plotting snapshot " << i << endl;

        plot_GNU(particle_states->snaps[i], metrics->metrics[i]); 

        fprintf(gnuplotPipe, "set label 1 'Step: %d' at screen 0.01,0.01\n", i);
        fflush(gnuplotPipe);
        
    }

    //6. Close the plot when the user is done

    cout << "Simulation completed. Close the plot window or press enter to exit." << endl;
    cin.ignore();
    cin.get();
    close_GNU();

}

//this function will initialize the plot using GNUplot. It is static as it only will be used by functions in this file
void Plotter::init_GNU(shared_ptr<scenario> scenario) {
    //initialize the plot using GNUplot (open a pipe)
    cout << "Initializing GNU";

    gnuplotPipe = popen("gnuplot -persistent", "w");

    //set some basic parameters for the plot, based on the scenario
    fprintf(gnuplotPipe, "set title '%s' font 'Arial Bold,16'\n", scenario->name.c_str());
    

    //set some basic parameters that are independent of the scenario
    fprintf(gnuplotPipe, "set xrange [-100:100]\n");
    fprintf(gnuplotPipe, "set yrange [-100:100]\n");
    //fprintf(gnuplotPipe, "set zrange [-100:100]\n");

    //set aspect ratio
    fprintf(gnuplotPipe, "set size ratio -1\n");

    //ensure objects are filled
    fprintf(gnuplotPipe, "set style fill solid 1.0 noborder\n");


    

    fflush(gnuplotPipe);

    

    
}

void Plotter::plot_GNU(shared_ptr<Particles> particles, shared_ptr<test_metrics_t> metrics_t) {
    //this function will plot the particles to the GNUplot window

    //1. reorder the particles, so they are visualized correctly for distance to the viewer (not required for 2D). Prior, assess viewers location.
    //to implement


    //2. plot the particles

    fprintf(gnuplotPipe, "plot '-' with circles lc rgb variable\n");
    for (int i = 0; i < particles->particle_list.size(); i++) {

        //set radius of the point using the rad field

        float point_size = static_cast<float>(particles->particle_list[i]->rad);

        //plot location, radius and rgb for the particle
        fprintf(gnuplotPipe, "%f %f %f %d\n", particles->particle_list[i]->x, particles->particle_list[i]->y, point_size, particles->particle_list[i]->rgb);

    }

    //tell gnuplot that the data input for this snapshot has ended
    fprintf(gnuplotPipe, "e\n");

    //3. update the metrics on the plot 
    //fprintf(gnuplotPipe, "set label 'Time: %d' at 0,0,0\n", metrics_t->time); //does not show up correctly, but placeholder for metrics to be implemented
    //fprintf(gnuplotPipe, "set label 'Memory: %d' at 0,0,0\n", metrics_t->memory);

    //4. refresh the plot
    
    fflush(gnuplotPipe);


}

void Plotter::close_GNU() {
    if (gnuplotPipe != nullptr) {
        pclose(gnuplotPipe);
        gnuplotPipe = nullptr;
    }
}


int Plotter::intensity_to_rgb(double r, double g, double b) {
    //this function will convert the rgb values to a hex code

    //1. convert the intensity values to 255 base

    int r255 = r * 255;
    int g255 = g * 255;
    int b255 = b * 255;

    //2. combine the rgb  

    int rgb = (r255 << 16) | (g255 << 8) | (b255);

    return rgb;

}


shared_ptr<Particles> Plotter::convert_intensity_to_rgb(shared_ptr<Particles> particles) {
    //this function will convert the rgb values of the particles to hex values

    for (int i = 0; i < particles->particle_list.size(); i++) {
        //convert the rgb values to hex using the rgb_to_hex function from Utils.cpp
        int rgb = intensity_to_rgb(particles->particle_list[i]->r, particles->particle_list[i]->g, particles->particle_list[i]->b);

        //store the hex value in the rgb field of the particle

        particles->particle_list[i]->rgb = rgb;

    }

    return particles;
}