#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <windows.h>
#include <algorithm>
#include <cmath>
#include <limits>

#include "../include/ParticlePlotter.h"
#include "../include/InitStructs.h"

using namespace std;

FILE *gnuplotPipe = nullptr;

Plotter::Plotter() {
    cout << "Plotting Engine initialized." << endl;
}

Plotter::~Plotter() {
    cout << "Plotting Engine destroyed." << endl;
}

void Plotter::plot_run(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {

    debug_mode_ = scenario->debug_mode;
    

    apply_heat_brightness_and_pack_rgb(particle_states, scenario);

    init_GNU(scenario);

    cout <<"Plotting the initial state..." << endl;

    //check if particle_states->snaps[0] and particle_states->metrics[0] are valid before plotting
    
    bool ps_state_exists = (particle_states->snaps.size() > 0 && particle_states->snaps[0] != nullptr);
    bool ps_metrics_exists = (particle_states->metrics.size() > 0 && particle_states->metrics[0] != nullptr);
    plot_GNU(particle_states->snaps[0], particle_states->metrics[0]);
    cout <<"Initial state plotted..." << endl;

    

    // Set N label once (particle count never changes)
    fprintf(gnuplotPipe, "set label 2 'N= %d' at screen 0.01,0.90 textcolor rgb 'white'\n",
            static_cast<int>(particle_states->snaps[0]->particle_list.size()));
    fprintf(gnuplotPipe, "set label 1 'Step: 0' at screen 0.01,0.95 textcolor rgb 'white'\n");
    fflush(gnuplotPipe);

    cout << "pipe flushed " << endl;

    cout << "............................................" << endl;
    cout << "Ready to plot scenario:" << scenario->name << endl;
    cout << "............................................" << endl;
    cout << "Press enter to start." << endl;
    cin.ignore();
    cin.get();

    int step = static_cast<int>(1.0 / static_cast<double>(scenario->dt)) * scenario->plot_speed_multiplier;
    if (step < 1) step = 1;

    int base_step = static_cast<int>(1.0 / static_cast<double>(scenario->dt));
    if (base_step < 1) base_step = 1;
    cout << "Playback: " << scenario->plot_speed_multiplier << "x speed (step " << step << " vs base " << base_step << ")" << endl;

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
void Plotter::apply_heat_brightness_and_pack_rgb(shared_ptr<snapshots> snaps, shared_ptr<scenario> scenario) {

    if (!snaps || snaps->snaps.empty() || !snaps->snaps[0] || snaps->snaps[0]->particle_list.empty()) {
        return;
    }

    // Config: cutoff, gamma
    const high_prec HEAT_CUTOFF_FRAC = scenario->heat_cutoff_frac;
    const high_prec HEAT_GAMMA_HP    = scenario->heat_gamma;

    // Use snapshot 0 as the authoritative "base RGB" reference
    const auto& base_frame = snaps->snaps[0];
    const int base_count   = static_cast<int>(base_frame->particle_list.size());

    // Find global min/max temperature across all frames (as before)
    high_prec max_temp = 0.0;
    high_prec min_temp = std::numeric_limits<high_prec>::max();

    for (int i = 0; i < static_cast<int>(snaps->snaps.size()); i++) {
        if (!snaps->snaps[i]) continue;
        for (int j = 0; j < static_cast<int>(snaps->snaps[i]->particle_list.size()); j++) {
            if (!snaps->snaps[i]->particle_list[j]) continue;
            high_prec t = snaps->snaps[i]->particle_list[j]->temp;
            if (t > max_temp) max_temp = t;
            if (t < min_temp) min_temp = t;
        }
    }

    if (min_temp == std::numeric_limits<high_prec>::max()) min_temp = 0.0;
    high_prec temp_range = max_temp - min_temp;

    auto pack_rgb = [](high_prec r, high_prec g, high_prec b) -> int {
        if (r < 0.0) r = 0.0; if (r > 1.0) r = 1.0;
        if (g < 0.0) g = 0.0; if (g > 1.0) g = 1.0;
        if (b < 0.0) b = 0.0; if (b > 1.0) b = 1.0;

        int r255 = static_cast<int>(r * 255);
        int g255 = static_cast<int>(g * 255);
        int b255 = static_cast<int>(b * 255);
        return (r255 << 16) | (g255 << 8) | (b255);
    };

    // If heating is effectively absent or unusable, just pack BASE RGB for all frames (and always for frame 0)
    if (max_temp < 0.00001 || temp_range <= 0.0) {
        for (int i = 0; i < static_cast<int>(snaps->snaps.size()); i++) {
            if (!snaps->snaps[i]) continue;
            for (int j = 0; j < static_cast<int>(snaps->snaps[i]->particle_list.size()); j++) {
                auto& p = snaps->snaps[i]->particle_list[j];
                if (!p) continue;

                // base RGB from frame 0 if possible
                high_prec br = (j < base_count && base_frame->particle_list[j]) ? base_frame->particle_list[j]->r : p->r;
                high_prec bg = (j < base_count && base_frame->particle_list[j]) ? base_frame->particle_list[j]->g : p->g;
                high_prec bb = (j < base_count && base_frame->particle_list[j]) ? base_frame->particle_list[j]->b : p->b;

                p->rgb = pack_rgb(br, bg, bb);
            }
        }
        return;
    }

    high_prec cutoff_temp = min_temp + HEAT_CUTOFF_FRAC * temp_range;
    high_prec denom = (max_temp - cutoff_temp);
    if (denom <= 0.0) denom = 1.0; // safety

    for (int i = 0; i < static_cast<int>(snaps->snaps.size()); i++) {
        if (!snaps->snaps[i]) continue;

        for (int j = 0; j < static_cast<int>(snaps->snaps[i]->particle_list.size()); j++) {
            auto& p = snaps->snaps[i]->particle_list[j];
            if (!p) continue;

            // --- ALWAYS use base RGB for frame 0 ---
            // And for all other frames, start from base RGB too (prevents any accidental "baking")
            high_prec base_r = (j < base_count && base_frame->particle_list[j]) ? base_frame->particle_list[j]->r : p->r;
            high_prec base_g = (j < base_count && base_frame->particle_list[j]) ? base_frame->particle_list[j]->g : p->g;
            high_prec base_b = (j < base_count && base_frame->particle_list[j]) ? base_frame->particle_list[j]->b : p->b;

            // Frame 0: no heat tinting, period.
            if (i == 0) {
                p->rgb = pack_rgb(base_r, base_g, base_b);
                continue;
            }

            // Heat fraction (0..1) above cutoff with gamma shaping
            high_prec fraction = 0.0;
            if (p->temp > cutoff_temp) {
                fraction = (p->temp - cutoff_temp) / denom;
                if (fraction < 0.0) fraction = 0.0;
                if (fraction > 1.0) fraction = 1.0;
                fraction = std::pow(fraction, HEAT_GAMMA_HP);
            }

            // Apply the existing "shift to red + brighten" look, but WITHOUT mutating p->r/g/b
            high_prec brighten = fraction * 0.25;

            high_prec out_r = base_r + (1.0 - base_r) * fraction + (1.0 - base_r) * brighten;
            high_prec out_g = base_g * (1.0 - fraction) + (1.0 - base_g) * brighten;
            high_prec out_b = base_b * (1.0 - fraction) + (1.0 - base_b) * brighten;

            p->rgb = pack_rgb(out_r, out_g, out_b);
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

        // Debugging: check pointers and sizes
        if (!particles) {
        std::cerr << "[plot_GNU] particles is nullptr!" << std::endl;
        return;
        }
        size_t pcount = particles->particle_list.size();
        std::cerr << "[plot_GNU] particle count: " << pcount << std::endl;
        for (size_t i = 0; i < pcount; ++i) {
        if (!particles->particle_list[i]) {
            std::cerr << "[plot_GNU] particle_list[" << i << "] is nullptr!" << std::endl;
        }
        }

        if (!metrics_t) {
        std::cerr << "[plot_GNU] metrics_t is nullptr!" << std::endl;
        }

        // Batch all gnuplot commands into a single string to minimize pipe write syscalls
        std::string buf;
        buf.reserve(pcount * 48 + 512);

        buf += "plot '-' with circles lc rgb variable\n";

        char line[128];
        for (int i = 0; i < static_cast<int>(pcount); i++) {
        if (!particles->particle_list[i]) {
            snprintf(line, sizeof(line), "# NULL PARTICLE %d\n", i);
            buf.append(line);
            continue;
        }
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
        if (debug_mode_ && metrics_t) {
        snprintf(line, sizeof(line), "set label 3 'KE  = %.4e' at screen 0.01,0.85 font 'Consolas,9' textcolor rgb 'white'\n",
            static_cast<double>(metrics_t->KE));
        buf += line;
        snprintf(line, sizeof(line), "set label 4 'PE  = %.4e' at screen 0.01,0.80 font 'Consolas,9' textcolor rgb 'white'\n",
            static_cast<double>(metrics_t->PE));
        buf += line;
        snprintf(line, sizeof(line), "set label 5 'HE  = %.4e' at screen 0.01,0.75 font 'Consolas,9' textcolor rgb 'white'\n",
            static_cast<double>(metrics_t->HE));
        buf += line;
        snprintf(line, sizeof(line), "set label 6 'TE  = %.4e' at screen 0.01,0.70 font 'Consolas,9' textcolor rgb 'white'\n",
            static_cast<double>(metrics_t->TE));
        buf += line;
        snprintf(line, sizeof(line), "set label 7 'Px  = %.4e' at screen 0.01,0.65 font 'Consolas,9' textcolor rgb 'white'\n",
            static_cast<double>(metrics_t->mom_x));
        buf += line;
        snprintf(line, sizeof(line), "set label 8 'Py  = %.4e' at screen 0.01,0.60 font 'Consolas,9' textcolor rgb 'white'\n",
            static_cast<double>(metrics_t->mom_y));
        buf += line;

        // Show relative TE error with color coding
        double rel_err = static_cast<double>(metrics_t->relative_error);
        const char* err_color = (rel_err < 0.001) ? "#00FF00" : (rel_err < 0.01) ? "#FFFF00" : "#FF4444";
        snprintf(line, sizeof(line), "set label 9 'TE err = %.6e' at screen 0.01,0.54 font 'Consolas,9' textcolor rgb '%s'\n",
            rel_err, err_color);
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
