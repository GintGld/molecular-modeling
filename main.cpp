#include "Gas/gas.h"
#include <iostream>
#include <chrono>

#define delimiter "----------------------\n"

using std::cout;

struct options_struct {
    fs::path input_file, output_extension;
    MF time, dt, target_temperature;
    bool last_state = false, ovito = false;

    bool f_time = false, f_dt = false, f_input = false, f_output = false;

    options_struct() {};

    bool is_full() const {
        return f_time && f_dt && f_input;
    }

    void update_output() {
        if (!f_output)
            output_extension = fs::path(".csv");
        return;
    }
};

bool compare(char* s1, std::string s2) {
    bool f = true;
    for (int i = 0; i < s2.size(); ++i)
        f &= (s1[i] == s2[i]);
    return f;
}

options_struct get_parameters(int argc, char** argv) {
    options_struct options;

    for (int i = 1; i < argc; ++i) {
        if ( compare(argv[i], "-t") ) {
            try {
                options.time = std::stod(argv[++i]);
                options.f_time = true;
            } catch(std::exception& e) {throw e;}
            continue;
        }
        if ( compare(argv[i], "-dt") ) {
            try {
                options.dt = std::stod(argv[++i]);
                options.f_dt = true;
            } catch(std::exception& e) {throw e;}
            continue;
        }
        if ( compare(argv[i], "--input") ) {
            options.input_file = argv[++i];
            options.f_input = true;
            continue;
        }
        if ( compare(argv[i], "--output") ) {
            options.output_extension = argv[++i];
            options.f_output = true;
            continue;
        }
        if ( compare(argv[i], "--last-state") ) {
            options.last_state = true;
            continue;
        }
        if ( compare(argv[i], "--without-centering-CM") ) {
            model::without_centering_CM = true;
            continue;
        }
        if ( compare(argv[i], "--scale") ) {
            try {
                options.target_temperature = std::stod(argv[++i]);
                model::scale = true;
            } catch(std::exception& e) {throw e;}
            continue;
        }
        if ( compare(argv[i], "--ovito") ) {
            options.ovito = true;
            continue;
        }
        throw std::runtime_error("Unrecognized option "+std::string(argv[i]));
    }

    if (!options.is_full()) {
        throw std::runtime_error ("Not enought parameters.\nRequired: -t, -dt, --input");
    }
    
    options.update_output();

    return options;
}

int main(int argc, char** argv) {
    parser p;

    options_struct options;

    try{
        options = get_parameters(argc, argv);
    } catch(std::exception& e) {
        cout << e.what() << "\n";
        return 1;
    }

    // defined directories for saving data
    fs::path input  = fs::path("./configs")/options.input_file;
    fs::path output = fs::path("./outputs")/options.input_file.stem().concat(options.output_extension.string());

    cout << "input  file: " << input.string() << "\n"
         << "output file: " << output.string() << "\n"
         << "simulation time: " << options.time << "\n"
         << "time step:       " << options.dt << "\n"
         << delimiter;

    try{
        p.read(input);
    } catch(std::exception& e) {
        cout << e.what() << "\n";
        return 1;
    }

    cout << "Read input file\n";

    model gas;

    try {
        gas = p.get_model();
    } catch(std::exception& e) {
        cout << e.what() << "\n";
        return 1;
    }

    cout << "Created gas model\n"
         << delimiter;

    cout << "Start simulation...\n";

    auto begin = std::chrono::steady_clock::now();

    gas.simulate(options.time, options.dt, options.target_temperature);

    auto end = std::chrono::steady_clock::now();

    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    std::cout << "Simulation time: " << elapsed_ms.count() / 1000 << " s\n";


    cout << "Finish simulation\n"
         << delimiter;

    cout << "Writing data...\n";

    gas.write(output);

    cout << "Data written to " << output.string() << "\n";

    cout << "Writing energies...\n";

    fs::path energy_file = fs::path(
        "./outputs/"+input.stem().string()+"-energy"+output.extension().string()
    );

    gas.write_energy(energy_file);

    cout << "Data written to " << energy_file.string() << "\n";

    if (options.ovito) {
        fs::path ovito_path = fs::path("./ovito")/input.stem().concat(".xyz");

        cout << "Writing .xyz file for ovito...\n";

        gas.write_ovito(ovito_path);

        cout << "XYZ file written to " << ovito_path.string() << "\n";
    }

    if (options.last_state) {
        fs::path last_state_path = input.parent_path()/
            fs::path(input.stem().string() + "-last-state" +
            input.extension().string());

        gas.write_last_state(last_state_path);

        cout << "Saved last state to " << last_state_path.string() << "\n";
    }
}