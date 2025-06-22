//Standard libraries
#include <memory> // Add this line to include the necessary header for shared_ptr
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <windows.h>  

//Internal libraries
#include "../include/ParticlePlotter.h"
#include "../include/InitStructs.h"

//namespaces
using namespace std;

//defining the gnuplot pipe
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
void Plotter::plot_run(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) { 
    
    //1. brighten the rgb values for particles with higher temperature
    particle_states = heat_to_rgb(particle_states);

    //2. convert r,g,b to hex
    for (int i = 0; i < particle_states->snaps.size(); i++) {

        
        

        particle_states->snaps[i] = convert_intensity_to_rgb(particle_states->snaps[i]);
    }

    


    cout << "RGB values have been calculated." << endl;

    //3. Initialize the plot using GNUplot

    init_GNU(scenario);

    plot_GNU(particle_states->snaps[0], particle_states->metrics[0]);




    fprintf(gnuplotPipe, "set label 1 'Step: 0' at screen 0.01,0.01\n");
    fflush(gnuplotPipe);

    //4. Confirm that the user is ready to start the simulation
    cout << "............................................" << endl;
    cout << "Ready to plot scenario:" << scenario->name << endl;
    cout << "............................................" << endl;
    cout << "Press enter to start." << endl;
    cin.ignore();
    cin.get();

    
    int frame_speed = static_cast<int>(1.0 / static_cast<double>(scenario->dt));
    int step = static_cast<int>(1.0 / static_cast<double>(scenario->dt));
    if (step < 1) step = 1; // Ensure at least step of 1

    //5. Loop through the snapshots and plot each n-th one on the plot
    for (int i = 0; i < static_cast<int>(particle_states->snaps.size()); i += step) {
        // Print progress every 10% of the simulation
        if (i % (particle_states->snaps.size() / 10) == 0) {
            cout << "Simulation " << (i * 100) / particle_states->snaps.size() << "% complete." << endl;
        }

        plot_GNU(particle_states->snaps[i], particle_states->metrics[i]);

        fprintf(gnuplotPipe, "set label 1 'Step: %d' at screen 0.01,0.95 textcolor rgb 'white'\n", i + 1);
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
    cout << "Initializing GNU" << endl;

    gnuplotPipe = popen("gnuplot -persistent", "w");

    //set some basic parameters for the plot, based on the scenario
    
    fprintf(gnuplotPipe, "set title '%s' font 'Arial Bold,16' textcolor rgb 'white'\n", scenario->name.c_str());
    

    //set some basic parameters that are independent of the scenario
    fprintf(gnuplotPipe, "set xrange [-100:100]\n");
    fprintf(gnuplotPipe, "set yrange [-100:100]\n");
    //fprintf(gnuplotPipe, "set zrange [-100:100]\n");

    //set aspect ratio
    fprintf(gnuplotPipe, "set size ratio -1\n");

    //ensure objects are filled
    fprintf(gnuplotPipe, "set style fill solid 1.0 noborder\n");

    //to have the plot background darkmode
    fprintf(gnuplotPipe, "set terminal wxt background '#000000'\n");

    //to hide the legend
    fprintf(gnuplotPipe, "unset key\n");

    fflush(gnuplotPipe);

    
}

void Plotter::plot_GNU(shared_ptr<Particles> particles, shared_ptr<test_metrics_t> metrics_t) {
    //this function will plot the particles to the GNUplot window

    //1. reorder the particles, so they are visualized correctly for distance to the viewer (not required for 2D). Prior, assess viewers location.
    //to implement


    //2. plot the particles

    fprintf(gnuplotPipe, "plot '-' with circles lc rgb variable\n");
    //fprintf(gnuplotPipe, "plot '-' with circles lc rgb variable fs empty\n");//for debugging you can use hollow circles
    for (int i = 0; i < particles->particle_list.size(); i++) {

        //set radius of the point using the rad field

        float point_size = static_cast<float>(particles->particle_list[i]->rad);
             
        //fprintf(gnuplotPipe, "%f %f %f %d\n", particles->particle_list[i]->x, particles->particle_list[i]->y, point_size, particles->particle_list[i]->rgb);
        //note that the above line will not work if x, y are high_prec, so we assume they are float or double
        //instead we can use the following line to convert them to float
        fprintf(gnuplotPipe, "%f %f %f %d\n", static_cast<float>(particles->particle_list[i]->x), static_cast<float>(particles->particle_list[i]->y), point_size, particles->particle_list[i]->rgb);
    }

    //tell gnuplot that the data input for this snapshot has ended
    fprintf(gnuplotPipe, "e\n");

    //3. update the metrics on the plot

    //print number of particles
    
    fprintf(gnuplotPipe, "set label 2 'N= %d' at screen 0.01,0.90 textcolor rgb 'white'\n", particles->particle_list.size());
    fprintf(gnuplotPipe, "set label 3 'FPS= %f' at screen 0.01,0.85 textcolor rgb 'white'\n", static_cast<float>(metrics_t->fps));
    fprintf(gnuplotPipe, "set label 4 'KE= %f K' at screen 0.01,0.80 textcolor rgb 'white'\n", static_cast<float>(metrics_t->KE/1000));
    fprintf(gnuplotPipe, "set label 5 'PE= %f K' at screen 0.01,0.75 textcolor rgb 'white'\n", static_cast<float>(metrics_t->PE/1000));
    fprintf(gnuplotPipe, "set label 6 'HE= %f K' at screen 0.01,0.70 textcolor rgb 'white'\n", static_cast<float>(metrics_t->HE/1000));
    fprintf(gnuplotPipe, "set label 7 'TE= %f K' at screen 0.01,0.65 textcolor rgb 'white'\n", static_cast<float>(metrics_t->TE/1000));
    fprintf(gnuplotPipe, "set label 8 'Mom x= %f K' at screen 0.01,0.60 textcolor rgb 'white'\n", static_cast<float>(metrics_t->mom_x/1000));
    fprintf(gnuplotPipe, "set label 9 'Mom y= %f K' at screen 0.01,0.55 textcolor rgb 'white'\n", static_cast<float>(metrics_t->mom_y/1000));
    //fprintf(gnuplotPipe, "set label 8 'Mom Change x= %f' at screen 0.01,0.60 textcolor rgb 'white'\n", metrics_t->mom_x_change);
    //fprintf(gnuplotPipe, "set label 9 'Mom Change y= %f' at screen 0.01,0.55 textcolor rgb 'white'\n", metrics_t->mom_y_change);
    fprintf(gnuplotPipe, "set label 10 'Relative TE Error= %f' at screen 0.01,0.50 textcolor rgb 'white'\n", static_cast<float>(metrics_t->relative_error));


    //format by dividing by /1000 and calling it K
    //fprintf(gnuplotPipe, "set label 8 'Aggregate TE Error: %f K' at screen 0.01,0.75\n", metrics_t->TE_error / 1000);

    //fprintf(gnuplotPipe, "set label 8 ' TE Error: %f %' at screen 0.01,0.70\n", metrics_t->relative_error);



    

    //4. refresh the plot
    
    fflush(gnuplotPipe);


}

void Plotter::close_GNU() {
    if (gnuplotPipe != nullptr) {
        pclose(gnuplotPipe);
        gnuplotPipe = nullptr;
    }
}


int Plotter::intensity_to_rgb(high_prec r, high_prec g, high_prec b) {
    //this function will convert the rgb values to a hex code

    //1. convert the intensity values to 255 base

    int r255 = static_cast<int>(r * 255);
    int g255 = static_cast<int>(g * 255);
    int b255 = static_cast<int>(b * 255);

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


shared_ptr<snapshots> Plotter::heat_to_rgb(shared_ptr<snapshots> snapshots) 
{
    // 1. Find the highest temperature of all particles in all snapshots.
    high_prec max_temp = 0.0;
    for (int i = 0; i < snapshots->snaps.size(); i++) {
        for (int j = 0; j < snapshots->snaps[i]->particle_list.size(); j++) {
            high_prec t = snapshots->snaps[i]->particle_list[j]->temp;
            if (t > max_temp) {
                max_temp = t;
            }
        }
    }

    high_prec min_temp = 0.0;
    high_prec temp_range = max_temp - min_temp;

    // If max_temp == 0 (or < 0), there's no brightness increase to do.
    if (temp_range <= 0.0) {
        return snapshots;
    }

    // 2. For each snapshot & each particle, linearly bring (r,g,b) closer to 1
    //    according to (temp - min_temp) / temp_range, clamping them at 1.0.
    for (int i = 0; i < snapshots->snaps.size(); i++) {
        for (int j = 0; j < snapshots->snaps[i]->particle_list.size(); j++) {

            high_prec t = snapshots->snaps[i]->particle_list[j]->temp;
            high_prec fraction = (t - min_temp) / temp_range;
            // Clamp the fraction between 0 and 1
            if (fraction < 0.0) fraction = 0.0;
            if (fraction > 1.0) fraction = 1.0;

            // Current color
            high_prec r = snapshots->snaps[i]->particle_list[j]->r;
            high_prec g = snapshots->snaps[i]->particle_list[j]->g;
            high_prec b = snapshots->snaps[i]->particle_list[j]->b;

            // Move each channel from its current value toward 1 by 'fraction'
            r = r + (1.0 - r) * fraction;
            g = g + (1.0 - g) * fraction;
            //b = b + (1.0 - b) * fraction;

            // Ensure each channel does not exceed 1.0
            if (r > 1.0) r = 1.0;
            if (g > 1.0) g = 1.0;
            if (b > 1.0) b = 1.0;

            // Write back the updated values
            snapshots->snaps[i]->particle_list[j]->r = r;
            snapshots->snaps[i]->particle_list[j]->g = g;
            snapshots->snaps[i]->particle_list[j]->b = b;
        }
    }

    return snapshots;
}