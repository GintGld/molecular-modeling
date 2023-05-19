#include "Gas/gas.h"
#include <iostream>

#define delimiter "----------------------\n"

using std::cout;

bool compare(char* s1, std::string s2) {
    bool f = true;
    for (int i = 0; i < s2.size(); ++i)
        f &= (s1[i] == s2[i]);
    return f;
}

std::tuple<fs::path, fs::path, MF, MF, bool> // {input_file, output_file, time, dt, last_state}
get_parameters(int argc, char** argv) {
    fs::path input_file, input_prefix, output_file, output_prefix;
    MF time, dt;
    bool f_time = false, f_dt = false, f_input = false, f_output = false;
    bool last_state = false;

    for (int i = 1; i < argc; ++i) {
        if ( compare(argv[i], "-t") ) {
            try {
                time = std::stod(argv[++i]);
                f_time = true;
            } catch(std::exception& e) {throw e;}
            continue;
        }
        if ( compare(argv[i], "-dt") ) {
            try {
                dt = std::stod(argv[++i]);
                f_dt = true;
            } catch(std::exception& e) {throw e;}
            continue;
        }
        if ( compare(argv[i], "--input") ) {
            input_file = argv[++i];
            f_input = true;
            continue;
        }
        if ( compare(argv[i], "--output") ) {
            output_file = argv[++i];
            f_output = true;
            continue;
        }
        if ( compare(argv[i], "--last-state") ) {
            last_state = true;
            continue;
        }
        throw std::runtime_error("Unrecognized option "+std::string(argv[i]));
    }

    if (!f_time || !f_dt || !f_input) {
        throw std::runtime_error ("Not enought parameters.\nRequired: -t, -dt, --input");
    }
    if (!f_output)
        output_file = fs::path(input_file.stem().string()+".csv");

    return {input_file, output_file, time, dt, last_state};
}

int main(int argc, char** argv) {
    parser p;

    fs::path input, output;
    MF time, dt;
    bool last_state = false;

    try{
        std::tie(input, output, time, dt, last_state) = get_parameters(argc, argv);
    } catch(std::exception& e) {
        cout << e.what() << "\n";
        return 1;
    }

    // defined directories for saving data
    input = fs::path("./configs")/input;
    output = fs::path("./outputs")/output;

    cout << "input file: " << input.string() << "\n"
         << "output file: " << output.string() << "\n"
         << "simulation time: " << time << "\n"
         << "time step: " << dt << "\n"
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

    gas.simulate(time, dt);

    cout << "Finish simulation\n"
         << delimiter;

    cout << "Writing data...\n";

    gas.write(output);

    cout << "Data written to " << output.string() << "\n";

    if (last_state) {
        fs::path last_state_path = input.parent_path()/
            fs::path(input.stem().string() + "-last-state" +
            input.extension().string());

        gas.write_last_state(last_state_path);

        cout << "Saved last state to " << last_state_path.string() << "\n";
    }
}