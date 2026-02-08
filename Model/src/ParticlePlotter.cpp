#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <windows.h>
#include <algorithm>
#include <cmath>

#include "../include/ParticlePlotter.h"
#include "../include/InitStructs.h"
#include "../include/PhysEngineWrapper.h"

using namespace std;

static const high_prec HEAT_GAMMA = 0.7;
static const int kPlotSpeedMultiplier = 20; // 1 = normal, 2 = 2x, etc.

FILE *gnuplotPipe = nullptr;

Plotter::Plotter() {
    cout << "Plotting Engine initialized." << endl;
}

Plotter::~Plotter() {
    cout << "Plotting Engine destroyed." << endl;
}

void Plotter::plot_run(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {

    apply_heat_brightness_and_pack_rgb(particle_states);

    cout << "RGB values have been calculated." << endl;

    init_GNU(scenario);
    plot_GNU(particle_states->snaps[0], particle_states->metrics[0]);

    // Set N label once (particle count never changes)
    fprintf(gnuplotPipe, "set label 2 'N= %d' at screen 0.01,0.90 textcolor rgb 'white'\n",
            static_cast<int>(particle_states->snaps[0]->particle_list.size()));
    fprintf(gnuplotPipe, "set label 1 'Step: 0' at screen 0.01,0.95 textcolor rgb 'white'\n");
    fflush(gnuplotPipe);

    cout << "............................................" << endl;
    cout << "Ready to plot scenario:" << scenario->name << endl;
    cout << "............................................" << endl;
    cout << "Press enter to start." << endl;
    cin.ignore();
    cin.get();

    int step = static_cast<int>(1.0 / static_cast<double>(scenario->dt)) * kPlotSpeedMultiplier;
    if (step < 1) step = 1;

    int base_step = static_cast<int>(1.0 / static_cast<double>(scenario->dt));
    if (base_step < 1) base_step = 1;
    cout << "Playback: " << kPlotSpeedMultiplier << "x speed (step " << step << " vs base " << base_step << ")" << endl;

    const int n = static_cast<int>(particle_states->snaps.size());
    const int totalIters = (n + step - 1) / step; // how many loop iterations will run

    for (int i = 0; i < n; i += step) {
        int iter = i / step; // 0..totalIters-1

        int progressStep = totalIters / 10;
        if (progressStep > 0 && (iter % progressStep == 0)) {
            cout << "Simulation " << (iter * 100) / totalIters << "% complete." << endl;
        }

        plot_GNU(particle_states->snaps[i], particle_states->metrics[i]);

        fprintf(gnuplotPipe, "set label 1 'Step: %d' at screen 0.01,0.95 textcolor rgb 'white'\n", i + 1);
        fflush(gnuplotPipe);
    }

    cout << "Simulation completed. Close the plot window or press enter to exit." << endl;
    cin.ignore();
    cin.get();
    close_GNU();
}
void Plotter::apply_heat_brightness_and_pack_rgb(shared_ptr<snapshots> snaps) {

    // Config: cutoff, gamma
    static const high_prec HEAT_CUTOFF_FRAC   = 0.00; // no cutoff — full temp range
    static const high_prec HEAT_GAMMA_HP      = HEAT_GAMMA;

    high_prec max_temp = 0.0;
    high_prec min_temp = 0.0;
    for (int i = 0; i < static_cast<int>(snaps->snaps.size()); i++) {
        for (int j = 0; j < static_cast<int>(snaps->snaps[i]->particle_list.size()); j++) {
            high_prec t = snaps->snaps[i]->particle_list[j]->temp;
            if (t > max_temp) max_temp = t;
        }
    }
    high_prec temp_range = max_temp - min_temp;
    if (temp_range <= 0.0) return;

    cout << "Max temp: " << max_temp << ", Min temp: " << min_temp << endl;

    high_prec cutoff_temp = min_temp + HEAT_CUTOFF_FRAC * temp_range;
    high_prec denom = (max_temp - cutoff_temp);
    if (denom <= 0.0) return;

    for (int i = 0; i < static_cast<int>(snaps->snaps.size()); i++) {
        for (int j = 0; j < static_cast<int>(snaps->snaps[i]->particle_list.size()); j++) {

            auto &p = snaps->snaps[i]->particle_list[j];

            high_prec fraction = 0.0;
            if (p->temp > cutoff_temp) {
                fraction = (p->temp - cutoff_temp) / denom;
                if (fraction < 0.0) fraction = 0.0;
                if (fraction > 1.0) fraction = 1.0;

                fraction = std::pow(fraction, HEAT_GAMMA_HP);
            }

            // Hue shift toward red + 1/4 as brightening
            high_prec brighten = fraction * 0.25;
            p->r = p->r + (1.0 - p->r) * fraction + (1.0 - p->r) * brighten;
            p->g = p->g * (1.0 - fraction) + (1.0 - p->g) * brighten;
            p->b = p->b * (1.0 - fraction) + (1.0 - p->b) * brighten;

            if (p->r > 1.0) p->r = 1.0;
            if (p->g > 1.0) p->g = 1.0;
            if (p->b > 1.0) p->b = 1.0;

            int r255 = static_cast<int>(p->r * 255);
            int g255 = static_cast<int>(p->g * 255);
            int b255 = static_cast<int>(p->b * 255);
            p->rgb = (r255 << 16) | (g255 << 8) | (b255);
        }
    }
}



void Plotter::init_GNU(shared_ptr<scenario> scenario) {
    cout << "Initializing GNU" << endl;

    gnuplotPipe = popen("gnuplot -persistent", "w");

    string scenario_name = scenario->name;
    replace(scenario_name.begin(), scenario_name.end(), '_', ' ');

    fprintf(gnuplotPipe, "set title '%s' font 'Arial Bold,16' textcolor rgb 'white'\n", scenario_name.c_str());
    fprintf(gnuplotPipe, "set xrange [-100:100]\n");
    fprintf(gnuplotPipe, "set yrange [-100:100]\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set style fill solid 1.0 noborder\n");
    fprintf(gnuplotPipe, "set terminal wxt noraise background '#000000'\n");
    fprintf(gnuplotPipe, "unset key\n");

    fflush(gnuplotPipe);
}

void Plotter::plot_GNU(shared_ptr<Particles> particles, shared_ptr<test_metrics_t> metrics_t) {

    // Batch all gnuplot commands into a single string to minimize pipe write syscalls
    std::string buf;
    buf.reserve(particles->particle_list.size() * 48 + 512);

    buf += "plot '-' with circles lc rgb variable\n";

    char line[128];
    for (int i = 0; i < static_cast<int>(particles->particle_list.size()); i++) {
        int len = snprintf(line, sizeof(line), "%f %f %f %d\n",
                static_cast<float>(particles->particle_list[i]->x),
                static_cast<float>(particles->particle_list[i]->y),
                static_cast<float>(particles->particle_list[i]->rad),
                particles->particle_list[i]->rgb);
        buf.append(line, len);
    }

    buf += "e\n";

    // N label is set once in plot_run (never changes)
    // Energy/momentum labels only in debug mode
    if (Engine::debug_mode()) {
        snprintf(line, sizeof(line), "set label 3 'FPS= %f' at screen 0.01,0.85 textcolor rgb 'white'\n",
                static_cast<float>(metrics_t->fps));
        buf += line;
        snprintf(line, sizeof(line), "set label 4 'KE= %f K' at screen 0.01,0.80 textcolor rgb 'white'\n",
                static_cast<float>(metrics_t->KE / 1000));
        buf += line;
        snprintf(line, sizeof(line), "set label 5 'PE= %f K' at screen 0.01,0.75 textcolor rgb 'white'\n",
                static_cast<float>(metrics_t->PE / 1000));
        buf += line;
        snprintf(line, sizeof(line), "set label 6 'HE= %f K' at screen 0.01,0.70 textcolor rgb 'white'\n",
                static_cast<float>(metrics_t->HE / 1000));
        buf += line;
        snprintf(line, sizeof(line), "set label 7 'TE= %f K' at screen 0.01,0.65 textcolor rgb 'white'\n",
                static_cast<float>(metrics_t->TE / 1000));
        buf += line;
        snprintf(line, sizeof(line), "set label 8 'Mom x= %f K' at screen 0.01,0.60 textcolor rgb 'white'\n",
                static_cast<float>(metrics_t->mom_x / 1000));
        buf += line;
        snprintf(line, sizeof(line), "set label 9 'Mom y= %f K' at screen 0.01,0.55 textcolor rgb 'white'\n",
                static_cast<float>(metrics_t->mom_y / 1000));
        buf += line;
        snprintf(line, sizeof(line), "set label 10 'Relative TE Error= %f' at screen 0.01,0.50 textcolor rgb 'white'\n",
                static_cast<float>(metrics_t->relative_error));
        buf += line;
    }

    fwrite(buf.data(), 1, buf.size(), gnuplotPipe);
    fflush(gnuplotPipe);
}

void Plotter::close_GNU() {
    if (gnuplotPipe != nullptr) {
        pclose(gnuplotPipe);
        gnuplotPipe = nullptr;
    }
}

int Plotter::intensity_to_rgb(high_prec r, high_prec g, high_prec b) {
    int r255 = static_cast<int>(r * 255);
    int g255 = static_cast<int>(g * 255);
    int b255 = static_cast<int>(b * 255);
    return (r255 << 16) | (g255 << 8) | (b255);
}

shared_ptr<Particles> Plotter::convert_intensity_to_rgb(shared_ptr<Particles> particles) {
    for (int i = 0; i < static_cast<int>(particles->particle_list.size()); i++) {
        int rgb = intensity_to_rgb(particles->particle_list[i]->r,
                                  particles->particle_list[i]->g,
                                  particles->particle_list[i]->b);
        particles->particle_list[i]->rgb = rgb;
    }
    return particles;
}

shared_ptr<snapshots> Plotter::heat_to_rgb(shared_ptr<snapshots> snapshots) {

    high_prec max_temp = 0.0;
    for (int i = 0; i < snapshots->snaps.size(); i++) {
        for (int j = 0; j < snapshots->snaps[i]->particle_list.size(); j++) {
            high_prec t = snapshots->snaps[i]->particle_list[j]->temp;
            if (t > max_temp) max_temp = t;
        }
    }

    high_prec min_temp = 0.0;
    high_prec temp_range = max_temp - min_temp;

    if (temp_range <= 0.0) {
        return snapshots;
    }

    for (int i = 0; i < snapshots->snaps.size(); i++) {
        for (int j = 0; j < snapshots->snaps[i]->particle_list.size(); j++) {

            high_prec t = snapshots->snaps[i]->particle_list[j]->temp;
            high_prec fraction = (t - min_temp) / temp_range;
            if (fraction < 0.0) fraction = 0.0;
            if (fraction > 1.0) fraction = 1.0;

            high_prec r = snapshots->snaps[i]->particle_list[j]->r;
            high_prec g = snapshots->snaps[i]->particle_list[j]->g;
            high_prec b = snapshots->snaps[i]->particle_list[j]->b;

            r = r + (1.0 - r) * fraction;
            g = g + (1.0 - g) * fraction;
            b = b + (1.0 - b) * fraction;

            if (r > 1.0) r = 1.0;
            if (g > 1.0) g = 1.0;
            if (b > 1.0) b = 1.0;

            snapshots->snaps[i]->particle_list[j]->r = r;
            snapshots->snaps[i]->particle_list[j]->g = g;
            snapshots->snaps[i]->particle_list[j]->b = b;
        }
    }

    return snapshots;
}
