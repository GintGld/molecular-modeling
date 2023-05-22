#include <iostream>
#include <fstream>
#include <exception>
#include "parser.h"

using std::cout;

parser::parser(): particle_number(0) {}

parser::parser(const fs::path& file): particle_number(0) {
    read(file);
}

void parser::read(const fs::path& _file) {
    file = fs::path(_file);

    std::ifstream reader(file);

    if (!reader.good()) {
        reader.close();
        throw std::runtime_error("Can't read file "+file.string());
    }

    try {
        reader >> model::size_of_box >> particle_number;
    } catch(std::exception& e) {throw e;}

    MF x, y, z, vx, vy, vz;

    for (int i = 0; i < particle_number; ++i) {
        try{
            reader >> x >> y >> z >> vx >> vy >> vz;
        } catch(std::exception& e) {throw e;}
        Particles.push_back(particle(x, y, z, vx, vy, vz));
    }

    reader.close();

    return;
}

model parser::get_model() const {
    model mod;

    for (const particle& p : Particles)
        mod.add_particle(p);

    //mod.particle_number = particle_number;
    //mod.Particles = Particles;

    return mod;
}