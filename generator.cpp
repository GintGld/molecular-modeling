#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <filesystem>

namespace fs=std::filesystem;

using std::cout;

/*
    ToDo
    add flag "-h"
*/

bool compare(char* s1, const std::string s2) {
    bool f = true;
    for (int i = 0; i < s2.size(); ++i)
        f &= (s1[i] == s2[i]);
    return f;
}

int main(int argc, char** argv) {

    fs::path file("file.mdl"), prefix("./configs");

    std::string sep = " ";

    double density, temperature = 1.;
    int cells, seed = time(NULL);

    bool f_density = false, f_cells = false, f_file = false;

    for (int i = 1; i < argc; ++i) {
        if ( compare(argv[i], "-d") ) {
            try {
                density = std::stod(argv[++i]);
                f_density = true;
            } catch(...) {
                cout << "incorrect parameter for density\n";
                return 1;
            }
            continue;
        }
        if ( compare(argv[i], "-c") ) {
            try {
                cells = std::stoi(argv[++i]);
                f_cells = true;
            } catch(...) {
                cout << "incorrect parameter for number of cells\n";
                return 1;
            }
            continue;
        }
        if ( compare(argv[i], "-f") ) {
            file = fs::path(argv[++i]);
            f_file = true;
            continue;
        }
        if ( compare(argv[i], "-t") ) {
            try {
                temperature = std::stod(argv[++i]);
            } catch(...) {
                cout << "Incorrect parameter for temperature\n";
                return 1;
            }
        }
        if ( compare(argv[i], "--prefix") ) {
            prefix = fs::path(argv[++i]);
            continue;
        }
        if ( compare(argv[i], "--seed") == 1 ) {
            try {
                seed = std::stoi(argv[++i]);
            } catch(...) {
                cout << "incorrect parameter for random seed\n";
                return 1;
            }
            continue;
        }
        if ( compare(argv[i], "--sep") == 1 ) {
            sep = argv[++i];
            continue;
        }
    }

    if (!f_density || !f_cells || !f_file) {
        cout << "Not enough parameters. Needed \"-d\", \"-c\", and \"-f\".\n";
        return 1;
    }

    std::default_random_engine gen(seed);
    std::normal_distribution<double> theta(0, M_PI), phi(0, 2 * M_PI);

    double box_size = (cells - 1) / cbrt(density);
    double distance = box_size / cells;

    cout << "Generating parameters:\n"
         << "Density:                    " << density << "\n"
         << "Cells:                      " << cells << "\n"
         << "Box size:                   " << box_size << "\n"
         << "Distance between particles: " << distance << "\n"
         << "Number of particles:        " << (cells - 1) * (cells - 1) * (cells - 1) << "\n"
         << "Temperature:                " << temperature << "\n"
         << "Mean square velocity:       " << sqrt(3 * temperature) << "\n";

    std::ofstream out(prefix/file);

    if (!out.good()) {
        cout << "Can't open file";
    }

    out << box_size << sep << (cells - 1) * (cells - 1) * (cells - 1) << "\n";

    //for (int i = 1; i < cells; ++i) for (int j = 1; j < cells; ++j) for (int k = 1; k < cells; ++k)
    //    out << i * distance << sep << j * distance << sep << k * distance << sep 
    //        << dist(gen) << sep << dist(gen) << sep << dist(gen) << "\n";

    double velocity = sqrt(3 * temperature);

    for (int i = 1; i < cells; ++i) for (int j = 1; j < cells; ++j) for (int k = 1; k < cells; ++k)
        out << i * distance << sep << j * distance << sep << k * distance << sep
            << velocity * sin(theta(gen)) * cos(phi(gen)) << sep
            << velocity * sin(theta(gen)) * sin(phi(gen)) << sep
            << velocity * cos(theta(gen)) << "\n";
    
    // technical info, model.cpp will not read this
    out << "DENS.: " << density << "\n"
        << "CELLS: " << cells << "\n"
        << "TEMP.: " << temperature << "\n"
        << "SEED:  " << seed << "\n";

    out.close();
}