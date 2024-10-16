
//Standard libraries
#include <iostream>
#include <fstream>
#include <memory>

//External libraries
#define CSV_IO_NO_THREAD
#include "../include/3party/csv.h"

//Internal libraries
#include "../include/Interfacer.h"

//namespaces
using namespace std;

// Constructor
Interfacer::Interfacer() {

    cout << "Interfacer initialized." << endl;
} 

// Destructor
Interfacer::~Interfacer() {
    cout << "Interfacer destroyed." << endl;
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

    const int column_count = 8;

    io::CSVReader<column_count, TrimPolicy, QuotePolicy> in(scenario_input_path);
    
    
    //3. Read the contents of the csv file and store each row in a scenario object. Store all scenario objects in a vector. 
    //scenario new_scenario;
    scenarios scenario_list;

    string col1, col2, col3, col4, col5, col6, col7, col8;


    int scenarios_loaded = 1;

    //discard the header row. This list is a guide to the columns in the csv file.
    in.read_header(io::ignore_extra_column, "SCENARIO_NAME", "OBJ_LIST", "TIME", "INTERACTION_FUNC", "TRY_CACHE", "REFRESH_OBJ", "TIMESTEP", "3D");

    while(in.read_row(col1, col2, col3, col4, col5, col6, col7, col8)) {

        unique_ptr<scenario> new_scenario(new scenario);


        new_scenario -> scenario_id  = scenarios_loaded;
        new_scenario-> name = col1;
        new_scenario-> obj_list = col2;

        new_scenario-> time = stod(col3);
        
        new_scenario-> interaction_func = col4;

        if (col5 == "TRUE") {
            new_scenario-> try_cache = true;
        } else {
            new_scenario-> try_cache = false;
        }
        
        if (col6 == "TRUE") {
            new_scenario-> refresh_obj = true;
        } else {
            new_scenario-> refresh_obj = false;
        }

        if (col8 == "TRUE") {
            new_scenario-> three_d = true;
        } else {
            new_scenario-> three_d = false;
        }

        new_scenario-> dt = stod(col7);
        
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
            cout << "Interaction Function: " << scenario_list.scenario_list[i]->interaction_func << endl;
            cout << "Try Cache: " << scenario_list.scenario_list[i]->try_cache << endl;
            cout << "Refresh Object: " << scenario_list.scenario_list[i]->refresh_obj << endl;
            cout << "Total time (s): " << scenario_list.scenario_list[i]->time << endl;
            cout << "Time Step length (s): " << scenario_list.scenario_list[i]->dt << endl;
            cout << "Steps: " << scenario_list.scenario_list[i]->time / scenario_list.scenario_list[i]->dt << endl;
            cout << "3D: " << scenario_list.scenario_list[i]->three_d << endl;
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

