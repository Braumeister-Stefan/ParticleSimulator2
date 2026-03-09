#include "../include/Interfacer.h"
#define CSV_IO_NO_THREAD
#include "../include/3party/csv.h"

#include <iostream>
#include <fstream>
#include <memory>



//namespaces
using namespace std;

// Constructor
Interfacer::Interfacer() {

   cout << "======================================================" << endl;
   cout << "  ParticleSimulator2 - By Braumeister Stefan (2026)   " << endl;
   cout << "======================================================" << endl << endl;

   cout << "What this simulator can do for you:" << endl;
   cout << "------------------------------------" << endl;
   cout << "  * Simulate gravitational attraction and kinetic collisions" << endl;
   cout << "    between n spherical particles in 2D space." << endl;
   cout << "  * Support both elastic and inelastic (energy-absorbing)" << endl;
   cout << "    collisions with configurable restitution parameters." << endl;
   cout << "  * Generate complex objects (e.g. particle circles / clouds)" << endl;
   cout << "    and save them to a cache for reuse across scenarios." << endl;
   cout << "  * Load previously computed object states as building blocks" << endl;
   cout << "    for new, more elaborate simulations." << endl;
   cout << "  * Visualize every timestep live via GNUplot, with color-coded" << endl;
   cout << "    temperature (heat glow) and adjustable playback speed." << endl;
   cout << "  * Track physics metrics: kinetic energy, momentum conservation," << endl;
   cout << "    total energy error, and frames-per-second performance." << endl;
   cout << "  * Compare simulation runs against saved benchmarks to verify" << endl;
   cout << "    correctness and measure computational performance." << endl;
   cout << endl;
   cout << "To get started, select a scenario from the list below." << endl;
   cout << endl;
} 

// Destructor
Interfacer::~Interfacer() {
    //cout << "Interfacer destroyed." << endl;
}

void Interfacer::setup_console_window() {
    //This function will set up the console window to be displayed on the left third of the screen. WIP


}



// Methods for interfacing with the code
shared_ptr<scenario> Interfacer::select_scenario() {
    //This function will load the list of scenarios and allow the user to select one. It will report the details of the selected scenario.
    setup_console_window();

    cout << "The following scenarios are available:" << endl;

    //1. Retrieve the scenario_inputs from csv file
    string scenario_input_path = "Inputs/scenario_inputs.csv";


    //2. Create csv reader object
    typedef io::trim_chars<' ', '\t'> TrimPolicy;
    typedef io::double_quote_escape<',', '\"'> QuotePolicy;

    const int column_count = 19;

    io::CSVReader<column_count, TrimPolicy, QuotePolicy> in(scenario_input_path);
    
    
    //3. Read the contents of the csv file and store each row in a scenario object. Store all scenario objects in a vector. 
    //scenario new_scenario;
    scenarios scenario_list;

    string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19;


    int scenarios_loaded = 1;

    //discard the header row. This list is a guide to the columns in the csv file.
    in.read_header(io::ignore_extra_column, "SCENARIO_NAME", "OBJ_LIST", "TIME", "DT", "BETA",
        "SAVED_OBJ", "COLLISION_DIST_TOL",
        "HEAT_GAMMA", "HEAT_CUTOFF_FRAC", "PLOT_SPEED_MULTIPLIER",
        "BENCHMARK_TE_ERROR_PCT", "BENCHMARK_SIM_TIME_SEC",
        "TRY_CACHE", "REFRESH_OBJ", "DEBUG_MODE", "REPORT_ENERGY_PER_STEP", "SAVE_SCENARIO", "GLOW_MODE", "RELAX_OVERLAPS");

    while(in.read_row(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19)) {

        unique_ptr<scenario> new_scenario(new scenario);


        new_scenario -> scenario_id  = scenarios_loaded;
        new_scenario-> name = col1;
        new_scenario-> obj_list = col2;

        new_scenario-> time = stod(col3);
        new_scenario-> dt = stod(col4);
        new_scenario-> beta = stod(col5);

        if (col6 != "FALSE") {
            new_scenario-> save_obj = col6;
        } else {
            new_scenario-> save_obj = "FALSE";
        }

        new_scenario-> collision_distance_tolerance = stod(col7);

        new_scenario-> heat_gamma = stod(col8);
        new_scenario-> heat_cutoff_frac = stod(col9);
        new_scenario-> plot_speed_multiplier = stoi(col10);

        new_scenario-> benchmark_te_error_pct = stod(col11);
        new_scenario-> benchmark_sim_time_sec = stod(col12);

        new_scenario-> try_cache = (col13 == "TRUE");
        new_scenario-> refresh_obj = (col14 == "TRUE");
        new_scenario-> debug_mode = (col15 == "TRUE");
        new_scenario-> report_energy_per_step = (col16 == "TRUE");
        new_scenario-> save_scenario = (col17 == "TRUE");
        new_scenario-> glow_mode = (col18 == "TRUE");
        new_scenario-> relax_overlaps = (col19 != "FALSE");


        
        
        //add the new scenario to the list of scenarios
        scenario_list.scenario_list.push_back(move(new_scenario));

        scenarios_loaded++;

    }

    //4. Print for each scenario in the list, the scenario_id and the name on a new line in the console.
    for (int i = 0; i < scenario_list.scenario_list.size(); i++) {
        cout << scenario_list.scenario_list[i] -> scenario_id << ". " << scenario_list.scenario_list[i] -> name << endl;
    }
    cout << endl;

    //5. Ask the user to select a scenario by entering the scenario_id. Store the selected scenario_id in a variable.

    int selected_scenario_id;

    cout << "Please enter the scenario ID of the scenario you would like to select: ";

    cin >> selected_scenario_id;
    cout << endl;

    //6. Print the details of the selected scenario to the console.

    for (int i = 0; i < scenario_list.scenario_list.size(); i++) {
        if (scenario_list.scenario_list[i] -> scenario_id == selected_scenario_id) {
            cout << "Selected Scenario Details: " << endl;
            cout << "Scenario ID: " << scenario_list.scenario_list[i]->scenario_id << endl;
            cout << "Name: " << scenario_list.scenario_list[i]->name << endl;
            cout << "Object List: " << scenario_list.scenario_list[i]->obj_list << endl;
            cout << "Time (s): " << scenario_list.scenario_list[i]->time << endl;
            cout << "Time Step Length: " << scenario_list.scenario_list[i]->dt << endl;
            cout << "Contact Bias Beta: " << scenario_list.scenario_list[i]->beta << endl;
            cout << "Saved OBJ (FALSE if not saved as object): " << scenario_list.scenario_list[i]->save_obj << endl;
            cout << "Collision Distance Tolerance: " << scenario_list.scenario_list[i]->collision_distance_tolerance << endl;
            cout << "Heat Gamma: " << scenario_list.scenario_list[i]->heat_gamma << endl;
            cout << "Heat Cutoff Frac: " << scenario_list.scenario_list[i]->heat_cutoff_frac << endl;
            cout << "Plot Speed Multiplier: " << scenario_list.scenario_list[i]->plot_speed_multiplier << endl;
            cout << "Benchmark TE Error (%): " << scenario_list.scenario_list[i]->benchmark_te_error_pct << endl;
            cout << "Benchmark Sim Time (s): " << scenario_list.scenario_list[i]->benchmark_sim_time_sec << endl;
            cout << "Try Cache: " << scenario_list.scenario_list[i]->try_cache << endl;
            cout << "Refresh Object: " << scenario_list.scenario_list[i]->refresh_obj << endl;
            cout << "Debug Mode: " << scenario_list.scenario_list[i]->debug_mode << endl;
            cout << "Report Energy Per Step: " << scenario_list.scenario_list[i]->report_energy_per_step << endl;
            cout << "Save Scenario: " << scenario_list.scenario_list[i]->save_scenario << endl;
            cout << "Glow Mode: " << scenario_list.scenario_list[i]->glow_mode << endl;
            cout << "Relax Overlaps: " << scenario_list.scenario_list[i]->relax_overlaps << endl;
        }
    }
    cout << endl;

    

    //7. Return the selected scenario_id to the main function by storing it in a public variable. in the Interfacer class.
    
    shared_ptr<scenario> selected_scenario = move(scenario_list.scenario_list[selected_scenario_id - 1]);


    //ask the user to press enter to continue
    //cout << "Press enter to continue." << endl;
    //cin.ignore();
    //cin.get();

    return selected_scenario;

}

