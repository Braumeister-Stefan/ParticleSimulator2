// ParticlePlotter.cpp

#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <windows.h>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <cstdlib>
#include <string>
#include <vector>
#include <type_traits>
#include <utility>

#include "../include/PhysEngineCore.h"

constexpr bool plot_debugmode = false; // Set to false to disable timing log

#include "../include/ParticlePlotter.h"
#include "../include/InitStructs.h"

using namespace std;

FILE *gnuplotPipe = nullptr;

// ------------------------------------
// Small internal helpers (mass detect)
// ------------------------------------
namespace {

template <typename T, typename = void>
struct has_member_m : std::false_type {};
template <typename T>
struct has_member_m<T, std::void_t<decltype(std::declval<T>().m)>> : std::true_type {};

template <typename T, typename = void>
struct has_member_mass : std::false_type {};
template <typename T>
struct has_member_mass<T, std::void_t<decltype(std::declval<T>().mass)>> : std::true_type {};

template <typename T, typename = void>
struct has_member_M : std::false_type {};
template <typename T>
struct has_member_M<T, std::void_t<decltype(std::declval<T>().M)>> : std::true_type {};

template <typename ParticleT>
static inline double get_mass_d(const ParticleT& p) {
    if constexpr (has_member_m<ParticleT>::value) {
        return static_cast<double>(p.m);
    } else if constexpr (has_member_mass<ParticleT>::value) {
        return static_cast<double>(p.mass);
    } else if constexpr (has_member_M<ParticleT>::value) {
        return static_cast<double>(p.M);
    } else {
        return 1.0; // fallback: equal-mass average
    }
}

static inline double abs_d(double v) { return (v < 0.0) ? -v : v; }

} // namespace

Plotter::Plotter() {

    // Small initial reserve; plot_GNU will grow this as needed
    frame_buf_.reserve(1024);
}

Plotter::~Plotter() {
    cout << "Plotting Engine destroyed." << endl;
}

void Plotter::plot_run(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {

    debug_mode_ = scenario->debug_mode;

    // ---- Plot bounds: COM + 98th-percentile from frame 0 ----
    bounds_from_com0(particle_states, plot_padding_frac_, plot_xmin_, plot_xmax_, plot_ymin_, plot_ymax_);
    plot_bounds_ready_ = true;

    // ---- OFFLINE render setup ----
    init_GNU(scenario);
    if (!gnuplotPipe) {
        cout << "ERROR: GNUplot pipe not available. Aborting render." << endl;
        return;
    }

    // Static plot settings that don't change per-frame should be sent ONCE.
    // (Saves pipe bytes + gnuplot parse work per frame.)
    {
        std::string once;
        once.reserve(256);
        append_ranges(once);
        if (!once.empty()) {
            fwrite(once.data(), 1, once.size(), gnuplotPipe);
            fflush(gnuplotPipe);
        }
    }

    // Set N label once (particle count never changes).
    fprintf(gnuplotPipe, "set label 2 'N= %d' at screen 0.01,0.90 font 'Arial,80' textcolor rgb 'white'\n",
            static_cast<int>(particle_states->snaps[0]->particle_list.size()));
    fflush(gnuplotPipe);

    cout << "............................................" << endl;
    cout << "Offline rendering scenario: " << scenario->name << endl;
    cout << "Output: " << offline_output_path_ << endl;
    cout << "............................................" << endl;

    cout << "Plot bounds (COM0+square): "
         << "x[" << plot_xmin_ << ", " << plot_xmax_ << "] "
         << "y[" << plot_ymin_ << ", " << plot_ymax_ << "]"
         << endl;

    int step = static_cast<int>(1.0 / static_cast<double>(scenario->dt)) * scenario->plot_speed_multiplier;
    if (step < 1) step = 1;

    int base_step = static_cast<int>(1.0 / static_cast<double>(scenario->dt));
    if (base_step < 1) base_step = 1;
    cout << "Playback: " << scenario->plot_speed_multiplier << "x speed (step " << step << " vs base " << base_step << ")" << endl;

    const int n = static_cast<int>(particle_states->snaps.size());
    const int totalIters = (n + step - 1) / step;

    // Progress logging setup (compute once)
    const int progressStep = totalIters / 20; // same behavior: if 0, no progress printing

    // Render frames (single pass; avoids rendering frame 0 twice)
    for (int i = 0; i < n; i += step) {
        const int iter = i / step; // 0..totalIters-1

        if (progressStep > 0 && (iter % progressStep == 0)) {
            cout << "Rendering " << (iter * 100) / totalIters << "% complete." << endl;
        }

        // Step label matches the snapshot index
        current_step_label_ = i;

        plot_GNU(particle_states->snaps[i], particle_states->metrics[i]);
    }

    // Close the output file to finalize the GIF, then quit gnuplot
    fprintf(gnuplotPipe, "set output\n");
    fprintf(gnuplotPipe, "quit\n");
    fflush(gnuplotPipe);

    close_GNU();

    cout << "Offline rendering completed." << endl;

    //print number of frames in gif and particles per frame
    cout << "Frames rendered: " << (n + step - 1) / step << " with " << particle_states->snaps[0]->particle_list.size() << " particles per frame." << endl;

    // ---- Autoplay after render completes ----
    playback_offline(offline_output_path_);
}

void Plotter::init_GNU(shared_ptr<scenario> scenario) {
    cout << "Initializing GNU (offline GIF renderer)" << endl;

    offline_output_path_ = "Inputs/rendered_scenarios/" + scenario->name + ".gif";

    gnuplotPipe = popen("gnuplot 2>nul", "w");
    if (!gnuplotPipe) {
        cout << "ERROR: Failed to start gnuplot process." << endl;
        return;
    }

    string scenario_name = scenario->name;
    replace(scenario_name.begin(), scenario_name.end(), '_', ' ');

    fprintf(gnuplotPipe, "set terminal gif animate delay 4 loop 0 size 6144,6144\n");
    fprintf(gnuplotPipe, "set output '%s'\n", offline_output_path_.c_str());

    fprintf(gnuplotPipe, "set object 1 rectangle from screen 0,0 to screen 1,1 behind fillcolor rgb '#000000' fillstyle solid 1.0\n");

    fprintf(gnuplotPipe, "set tmargin 14\n");
    fprintf(gnuplotPipe, "set bmargin 2\n");
    fprintf(gnuplotPipe, "set lmargin 2\n");
    fprintf(gnuplotPipe, "set rmargin 2\n");
    
    fprintf(gnuplotPipe, "set title '%s' font 'Arial,104' textcolor rgb 'white' offset 0,-7\n", scenario_name.c_str());
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set style fill solid 1.0 noborder\n");
    fprintf(gnuplotPipe, "unset key\n");

    fflush(gnuplotPipe);
}

void Plotter::plot_GNU(shared_ptr<Particles> particles, shared_ptr<test_metrics_t> metrics_t) {

    EngineCore::clock_t::time_point t_total_start, t_build_end, t_io_end;
    if (plot_debugmode) {
        t_total_start = EngineCore::clock_t::now();
    }

    const auto& plist = particles->particle_list;
    const int np = static_cast<int>(plist.size());

    // Reuse one buffer to reduce alloc churn.
    // Rough estimate: ~50-70 bytes per particle line + labels.
    const size_t need = static_cast<size_t>(np) * 64u + 1024u;
    if (frame_buf_.capacity() < need) frame_buf_.reserve(need);
    frame_buf_.clear();

    // Labels must be set BEFORE the plot command so they appear in the same frame.
    {
        char line[256];

        snprintf(line, sizeof(line),
                 "set label 1 'Step: %d' at screen 0.01,0.95 font 'Arial,80' textcolor rgb 'white'\n",
                 current_step_label_);
        frame_buf_ += line;

        if (debug_mode_ && metrics_t) {
            snprintf(line, sizeof(line),
                     "set label 3 'KE  = %.4e' at screen 0.01,0.85 font 'Consolas,60' textcolor rgb 'white'\n",
                     static_cast<double>(metrics_t->KE));
            frame_buf_ += line;

            snprintf(line, sizeof(line),
                     "set label 4 'PE  = %.4e' at screen 0.01,0.80 font 'Consolas,60' textcolor rgb 'white'\n",
                     static_cast<double>(metrics_t->PE));
            frame_buf_ += line;

            snprintf(line, sizeof(line),
                     "set label 5 'HE  = %.4e' at screen 0.01,0.75 font 'Consolas,60' textcolor rgb 'white'\n",
                     static_cast<double>(metrics_t->HE));
            frame_buf_ += line;

            snprintf(line, sizeof(line),
                     "set label 6 'TE  = %.4e' at screen 0.01,0.70 font 'Consolas,60' textcolor rgb 'white'\n",
                     static_cast<double>(metrics_t->TE));
            frame_buf_ += line;

            snprintf(line, sizeof(line),
                     "set label 7 'Px  = %.4e' at screen 0.01,0.65 font 'Consolas,60' textcolor rgb 'white'\n",
                     static_cast<double>(metrics_t->mom_x));
            frame_buf_ += line;

            snprintf(line, sizeof(line),
                     "set label 8 'Py  = %.4e' at screen 0.01,0.60 font 'Consolas,60' textcolor rgb 'white'\n",
                     static_cast<double>(metrics_t->mom_y));
            frame_buf_ += line;

            const double rel_err = static_cast<double>(metrics_t->relative_error);
            const char* err_color = (rel_err < 0.001) ? "#00FF00" : (rel_err < 0.01) ? "#FFFF00" : "#FF4444";
            snprintf(line, sizeof(line),
                     "set label 9 'TE err = %.6e' at screen 0.01,0.54 font 'Consolas,60' textcolor rgb '%s'\n",
                     rel_err, err_color);
            frame_buf_ += line;
        }
    }

    // Plot command + inline data (circles)
    frame_buf_ += "plot '-' with circles lc rgb variable\n";

    // Build particle lines
    char line[192];
    for (int i = 0; i < np; ++i) {
        const auto& p = plist[i];

        const double x   = static_cast<double>(p->x);
        const double y   = static_cast<double>(p->y);
        const double rad = static_cast<double>(p->rad);

        const int len = snprintf(line, sizeof(line), "%.10g %.10g %.10g %d\n",
                                 x, y, rad, p->rgb);
        frame_buf_.append(line, static_cast<size_t>(len));
    }
    frame_buf_ += "e\n";

    if (plot_debugmode) {
        t_build_end = EngineCore::clock_t::now();
    }

    if (gnuplotPipe) {
        fwrite(frame_buf_.data(), 1, frame_buf_.size(), gnuplotPipe);
        fflush(gnuplotPipe);
    }

    if (plot_debugmode) {
        t_io_end = EngineCore::clock_t::now();

        const double total_s = EngineCore::seconds_between(t_total_start, t_io_end);
        const double io_s    = EngineCore::seconds_between(t_build_end,   t_io_end);
        const double io_pct  = (total_s > 0.0) ? (io_s / total_s) * 100.0 : 0.0;

        std::cout << "plot_GNU took " << total_s
                  << " seconds (" << io_pct << "% due to GNUplot interaction)"
                  << std::endl;
    }
}

void Plotter::playback_offline(const std::string& rendered_file) {
    cout << "Press enter to start GIF " << rendered_file << endl;
    cin.get();

#ifdef _WIN32
    std::string cmd = "start \"\" \"" + rendered_file + "\"";
    system(cmd.c_str());
#elif __APPLE__
    std::string cmd = "open \"" + rendered_file + "\"";
    system(cmd.c_str());
#else
    std::string cmd = "xdg-open \"" + rendered_file + "\"";
    system(cmd.c_str());
#endif

    cout << "Press enter to close programme" << endl;
    cin.get();
}

void Plotter::close_GNU() {
    if (gnuplotPipe != nullptr) {
        pclose(gnuplotPipe);
        gnuplotPipe = nullptr;
    }
}

int Plotter::intensity_to_rgb(double r, double g, double b) {
    const int r255 = static_cast<int>(r * 255.0);
    const int g255 = static_cast<int>(g * 255.0);
    const int b255 = static_cast<int>(b * 255.0);
    return (r255 << 16) | (g255 << 8) | (b255);
}

shared_ptr<Particles> Plotter::convert_intensity_to_rgb(shared_ptr<Particles> particles) {
    auto& plist = particles->particle_list;
    const int np = static_cast<int>(plist.size());
    for (int i = 0; i < np; ++i) {
        auto &p = plist[i];
        const double r = static_cast<double>(p->r);
        const double g = static_cast<double>(p->g);
        const double b = static_cast<double>(p->b);
        p->rgb = intensity_to_rgb(r, g, b);
    }
    return particles;
}

shared_ptr<snapshots> Plotter::heat_to_rgb(shared_ptr<snapshots> snapshots) {

    double max_temp = 0.0;
    auto& frames = snapshots->snaps;
    const int nf = static_cast<int>(frames.size());

    for (int i = 0; i < nf; ++i) {
        auto& plist = frames[i]->particle_list;
        const int np = static_cast<int>(plist.size());
        for (int j = 0; j < np; ++j) {
            const auto& p = plist[j];
            const double t = static_cast<double>(p->temp);
            if (t > max_temp) max_temp = t;
        }
    }

    const double min_temp = 0.0;
    const double temp_range = max_temp - min_temp;

    if (temp_range <= 0.0) {
        return snapshots;
    }

    for (int i = 0; i < nf; ++i) {
        auto& plist = frames[i]->particle_list;
        const int np = static_cast<int>(plist.size());
        for (int j = 0; j < np; ++j) {
            auto &p = plist[j];

            const double t = static_cast<double>(p->temp);
            double fraction = (t - min_temp) / temp_range;
            if (fraction < 0.0) fraction = 0.0;
            if (fraction > 1.0) fraction = 1.0;

            double r = static_cast<double>(p->r);
            double g = static_cast<double>(p->g);
            double b = static_cast<double>(p->b);

            r = r + (1.0 - r) * fraction;
            g = g + (1.0 - g) * fraction;
            b = b + (1.0 - b) * fraction;

            if (r > 1.0) r = 1.0;
            if (g > 1.0) g = 1.0;
            if (b > 1.0) b = 1.0;

            p->r = r;
            p->g = g;
            p->b = b;
        }
    }

    return snapshots;
}

// ---------------------------------------------------------
// Plot bounds helpers – COM + 98th-percentile from frame 0
// ---------------------------------------------------------

void Plotter::bounds_from_com0(shared_ptr<snapshots> snaps,
                               double pad_frac,
                               double& xmin,
                               double& xmax,
                               double& ymin,
                               double& ymax) const {
    xmin = -1.0; xmax = 1.0;
    ymin = -1.0; ymax = 1.0;

    if (!snaps || snaps->snaps.empty() || !snaps->snaps[0]) return;

    const auto& frame0 = snaps->snaps[0];
    if (frame0->particle_list.empty()) return;

    const auto& plist = frame0->particle_list;
    const int np = static_cast<int>(plist.size());

    double sum_m  = 0.0;
    double sum_mx = 0.0;
    double sum_my = 0.0;

    for (int i = 0; i < np; ++i) {
        const auto& p = plist[i];
        if (!p) continue;

        const double m = get_mass_d(*p);
        const double x = static_cast<double>(p->x);
        const double y = static_cast<double>(p->y);

        if (!std::isfinite(m) || !std::isfinite(x) || !std::isfinite(y)) continue;

        sum_m  += m;
        sum_mx += m * x;
        sum_my += m * y;
    }

    double cx = 0.0;
    double cy = 0.0;

    if (sum_m != 0.0) {
        cx = sum_mx / sum_m;
        cy = sum_my / sum_m;
    } else {
        double sx = 0.0, sy = 0.0, cnt = 0.0;
        for (int i = 0; i < np; ++i) {
            const auto& p = plist[i];
            if (!p) continue;
            const double x = static_cast<double>(p->x);
            const double y = static_cast<double>(p->y);
            if (!std::isfinite(x) || !std::isfinite(y)) continue;
            sx += x;
            sy += y;
            cnt += 1.0;
        }
        if (cnt != 0.0) {
            cx = sx / cnt;
            cy = sy / cnt;
        }
    }

    // --------------
    // Radial percentile (98%) of square-radius-from-center
    // d_i = max(|x-cx|+rad, |y-cy|+rad)
    // --------------
    std::vector<double> dists;
    dists.reserve(static_cast<size_t>(np));

    for (int i = 0; i < np; ++i) {
        const auto& p = plist[i];
        if (!p) continue;

        const double x = static_cast<double>(p->x);
        const double y = static_cast<double>(p->y);
        const double r = static_cast<double>((p->rad > 0) ? p->rad : 0);

        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(r)) continue;

        const double dx = abs_d(x - cx) + r;
        const double dy = abs_d(y - cy) + r;
        const double d  = (dx > dy) ? dx : dy;

        if (!std::isfinite(d)) continue;
        dists.push_back(d);
    }

    double half = 0.0;
    if (!dists.empty()) {
        const double q = 0.95; // requested 98% radial percentile
        const size_t n = dists.size();

        // Use a conservative index so that at least ~98% are within bounds.
        size_t k = static_cast<size_t>(std::ceil(q * static_cast<double>(n - 1)));
        if (k >= n) k = n - 1;

        std::nth_element(dists.begin(), dists.begin() + static_cast<std::ptrdiff_t>(k), dists.end());
        half = dists[k];
    }

    if (half <= 0.0) half = 1.0;

    if (pad_frac < 0.0) pad_frac = 0.0;
    half *= (1.0 + 10.0 * pad_frac);

    xmin = cx - half;
    xmax = cx + half;
    ymin = cy - half;
    ymax = cy + half;

    // ---- Width ratio: last frame vs first frame ----
    const double width_first = std::abs(xmax - xmin);

    if (snaps->snaps.size() > 1) {
        const auto& last = snaps->snaps.back();
        if (last && !last->particle_list.empty()) {
            double mn_x =  std::numeric_limits<double>::max();
            double mx_x = -std::numeric_limits<double>::max();
            double mn_y =  std::numeric_limits<double>::max();
            double mx_y = -std::numeric_limits<double>::max();

            for (const auto& p : last->particle_list) {
                if (!p) continue;
                const double x = static_cast<double>(p->x);
                const double y = static_cast<double>(p->y);
                const double r = static_cast<double>((p->rad > 0) ? p->rad : 0);
                if (!std::isfinite(x) || !std::isfinite(y)) continue;
                if (x - r < mn_x) mn_x = x - r;
                if (x + r > mx_x) mx_x = x + r;
                if (y - r < mn_y) mn_y = y - r;
                if (y + r > mx_y) mx_y = y + r;
            }

            if (mn_x < mx_x && mn_y < mx_y) {
                const double span_x_last = mx_x - mn_x;
                const double span_y_last = mx_y - mn_y;
                const double width_last  = (span_x_last > span_y_last) ? span_x_last : span_y_last;

                if (width_first > 0.0) {
                    cout << "Last/First frame width ratio: "
                         << std::abs(width_last) / std::abs(width_first)
                         << endl << endl;
                }
            }
        }
    }
}

void Plotter::append_ranges(std::string& buf) const {
    if (!plot_bounds_ready_) return;

    char line[256];

    snprintf(line, sizeof(line), "set xrange [%.17g:%.17g]\n", plot_xmin_, plot_xmax_);
    buf += line;

    snprintf(line, sizeof(line), "set yrange [%.17g:%.17g]\n", plot_ymin_, plot_ymax_);
    buf += line;
}