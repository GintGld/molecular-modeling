#include "model.h"
#include "progressbar.hpp"

#include <iomanip> // std::setprecision
#include <cmath>
#include <tuple>
#include <fstream>

MF model::size_of_box = 10;
bool model::without_centering_CM = false;
MF model::scale = false;

model::model() : 
    particle_number(0), 
    potential_energy_tmp(0), 
    kinetic_energy_tmp(0) {
}

void model::add_particle(const particle& p) {
    Particles.push_back(p);
    ++particle_number;

    // recalculate potential energy
    update_acceleration();
    kinetic_energy_tmp += (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz) / 2;

}

void model::correct_CM() {
    /*
        scale velocities
        to make center of mass motionless
    */
    MF Px = 0, Py = 0, Pz = 0;
    for (const particle& p : Particles) {
        Px += p.vx;
        Py += p.vy;
        Pz += p.vz;
    }

    // substract center mass kinetic energy
    kinetic_energy_tmp -= (Px * Px + Py * Py + Pz * Pz) / 2;

    Px /= particle_number;
    Py /= particle_number;
    Pz /= particle_number;
    for (particle& p : Particles) {
        p.vx -= Px;
        p.vy -= Py;
        p.vz -= Pz;
    }

    return;
}

std::tuple<MF, MF, MF, MF>
model::nearest_reflection(const particle& p, const particle& q) {
    MF L = size_of_box, d;
    int kx, ky, kz;

    for (int i = -1; i <= 1; ++i) for (int j = -1; j <= 1; ++j) for (int k = -1; k <= 1; ++k) {
        d = particle::dist(p, {q.x + i * size_of_box, q.y + j * size_of_box, q.z + k * size_of_box});
        if (L > d) {
            L = d;
            kx = i;
            ky = j;
            kz = k;
        }
    }

    return {
        L, 
        p.x - q.x - kx * size_of_box, 
        p.y - q.y - ky * size_of_box,
        p.z - q.z - kz * size_of_box
    };
}

void model::update_acceleration() {
    /*
        calculate acceleration
        and potential energy
        for all particles
    */

    // dump accelerations
    for (particle& p : Particles) {
        p.wx = 0;
        p.wy = 0;
        p.wz = 0;
    }

    potential_energy_tmp = 0;

    for (int i = 0; i < particle_number; ++i) {
        for (int j = i + 1; j < particle_number; ++j) {
            // Get tuple {L, x, y, z}, unpuck it
            // (x, y, z) -- vector between 2 particles
            const auto [L, x, y, z] = nearest_reflection(Particles[i], Particles[j]);

            MF k = 24 * ( 2 * pow(L, -14) - pow(L, -8) );
            MF wx = x * k;
            MF wy = y * k;
            MF wz = z * k;
            Particles[i].update_w( wx,  wy,  wz);
            Particles[j].update_w(-wx, -wy, -wz);

            potential_energy_tmp += 4 * ( pow(L, -12) - pow(L, -6) );
        }
    }
}

void model::update_coordinates(MF dt) {
    /*
        move particles and
        calculate kinetic energy
    */

    // First integration half-step
    for (particle& p : Particles) {
        p.vx += .5 * p.wx * dt;
        p.vy += .5 * p.wy * dt;
        p.vz += .5 * p.wz * dt;
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;

        // handle borders
        if (p.x > size_of_box) p.x -= size_of_box;
        if (p.y > size_of_box) p.y -= size_of_box;
        if (p.z > size_of_box) p.z -= size_of_box;
        if (p.x < 0) p.x += size_of_box;
        if (p.y < 0) p.y += size_of_box;
        if (p.z < 0) p.z += size_of_box;
    }

    update_acceleration();

    kinetic_energy_tmp = 0;

    // Second integration half-step
    for (particle& p : Particles) {
        p.vx += .5 * p.wx * dt;
        p.vy += .5 * p.wy * dt;
        p.vz += .5 * p.wz * dt;

        kinetic_energy_tmp += (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz) / 2;
    }
}

void model::scale_particles(MF target_temp) {
    /*
        Scales particles' velocities
    */
    MF coeff = sqrt(target_temp * 1.5 * particle_number / kinetic_energy_tmp);

    for (particle& p : Particles) {
        p.vx *= coeff;
        p.vy *= coeff;
        p.vz *= coeff;
    }

    return;
}

void model::init_step(MF dt) {
    /*
        initial step when 
        there is no previous state
    */

    // acceleration is actual because of model::add_particle(const particle&)


    MF kitetic_energy_tmp = 0;

    for (particle& p : Particles) {
        p.x += p.vx * dt + 0.5 * p.wx * dt * dt;
        p.y += p.vy * dt + 0.5 * p.wy * dt * dt;
        p.z += p.vz * dt + 0.5 * p.wz * dt * dt;
        p.vx += p.wx * dt;
        p.vy += p.wy * dt;
        p.vz += p.wz * dt;

        kitetic_energy_tmp += (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz) / 2;
    }
}

void model::make_step(MF dt) {
    update_coordinates(dt);
}

void model::simulate(MF time, MF dt, MF target_temp) {
    /*
        simulation loop
    */

    if (!without_centering_CM)
        correct_CM();

    // save initial state
    commit();

    // make first step
    init_step(dt);

    commit();

    progressbar pBar(int (time / dt) - 2);

    // main algorithm
    for (size_t i = 2; dt * i < time; ++i) {
        make_step(dt);

        // scaling
        if (scale)
            scale_particles(target_temp);

        commit();

        pBar.update();
    }
    std::cout << "\n";
}

void model::commit() {
    /*
        save current system state as
        (x, y, z, vx, vy, vz), (...), ...
    */
    for (const particle& p : Particles) {
        hist_tmp.push_back(p.x);
        hist_tmp.push_back(p.y);
        hist_tmp.push_back(p.z);
        hist_tmp.push_back(p.vx);
        hist_tmp.push_back(p.vy);
        hist_tmp.push_back(p.vz);
    }
    history.push_back(hist_tmp);
    hist_tmp.clear();
    kinetic_energy.push_back(kinetic_energy_tmp);
    potential_energy.push_back(potential_energy_tmp);
}

std::string model::header() const {
    std::string s, res = "";
    for (int i = 0; i < particle_number; ++i) {
        s = std::to_string(i+1);
        res += "x"+s+",y"+s+",z"+s+",vx"+s+",vy"+s+",vz"+s;
        if (i != particle_number - 1)
            res += ",";
        else
            res += "\n";
    }
    return res;
}

void model::write(const fs::path& file) const {
    /*
        write history to a file
    */
    std::ofstream out;

    fs::path ext = file.extension();

    progressbar pBar(history.size());

    // writing into .csv format (for pandas)
    if (ext == ".csv") {
        out.open(file);

        if (!out.good()) {
            out.close();
            throw std::runtime_error("Can't open file "+file.string());
        }

        out << header();

        for (const std::vector<MF>& v : history) {
            for (int i = 0; i < v.size() - 1; ++i) 
                out << v[i] << ',';
            out << v.back() << "\n";

            pBar.update();
        }

        out.close();

        std::cout << "\n";
        return;
    }

    // writing into binary format (read by numpy.fromfile in python)
    if (ext == ".dat" || ext == ".data" || ext == ".binary" || ext == "") {
        out.open(file, std::ios::binary);

        if (!out.good()) {
            out.close();
            throw std::runtime_error("Can't open file"+file.string());
        }

        for (const std::vector<MF>& v : history) {
            for (int i = 0; i < v.size(); ++i)
                out.write((char*)(&v[i]), sizeof(MF));

            pBar.update();
        }

        out.close();

        std::cout << "\n";
        return;
    }
    
    out.close();
}

void model::write_last_state(const fs::path& file, const std::string& sep) const {
    std::ofstream out(file);

    if (!out.good()) {
        out.close();
        throw std::runtime_error ("Can't open file "+file.string());
    }

    out << size_of_box << sep << particle_number << "\n";

    for (const particle& p : Particles) {
        out << p.x << sep 
            << p.y << sep 
            << p.z << sep 
            << p.vx << sep 
            << p.vy << sep 
            << p.vz << "\n";
    }

    out.close();
}

void model::write_energy(const fs::path& file, const std::string& sep) const {
    /*
        write Energies to a file
    */
    std::ofstream out;

    fs::path ext = file.extension();

    progressbar pBar(kinetic_energy.size());

    // writing into .csv format (for pandas)
    if (ext == ".csv") {
        out.open(file);

        if (!out.good()) {
            out.close();
            throw std::runtime_error("Can't open file "+file.string());
        }

        out << "K" << sep << "P\n";

        for (int i = 0; i < kinetic_energy.size(); ++i) {
            out << kinetic_energy[i] << sep << potential_energy[i] << "\n";

            pBar.update();
        }

        out.close();

        std::cout << "\n";
        return;
    }

    // writing into binary format (read by numpy.fromfile in python)
    if (ext == ".dat" || ext == ".data" || ext == ".binary" || ext == "") {
        out.open(file, std::ios::binary);

        if (!out.good()) {
            out.close();
            throw std::runtime_error("Can't open file"+file.string());
        }

        for (int i = 0; i < kinetic_energy.size(); ++i) {
            out.write((char*)(&kinetic_energy[i]), sizeof(MF));
            out.write((char*)(&potential_energy[i]), sizeof(MF));

            pBar.update();
        }

        out.close();

        std::cout << "\n";
        return;
    }
    
    out.close();
}

void model::write_ovito(const fs::path& file) const {
    std::ofstream out(file);

    if (!out.good()) {
        out.close();
        throw std::runtime_error("Can't open "+file.string());
    }

    progressbar pBar(history.size());

    for (size_t i = 0; i < history.size(); ++i) {
        out << std::fixed << std::setprecision(8) << particle_number << "\n\n";

        for (int k = 0; k < particle_number; ++k)
            out << "1 " << history[i][6*k] << ' ' << history[i][6*k+1] << ' ' << history[i][6*k+2] << "\n";

        pBar.update();
    }
    std::cout << "\n";

    out.close();
    return;
}

const std::vector<particle>& model::get_particles() const {
    return Particles;
}

std::vector<particle>& model::get_particles() {
    return Particles;
}

const std::vector< std::vector<MF> >& model::get_history() {
    return history;
}