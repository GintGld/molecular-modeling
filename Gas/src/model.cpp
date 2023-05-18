#include "model.h"
#include <cmath>
#include <tuple>
#include <fstream>

MF model::size_of_box = 10, model::relaxation_time = 1, model::record_time = 1;

model::model(): particle_number(0) {
}

void model::add_particle(const particle& p) {
    Particles.push_back(p);
    ++particle_number;
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
        for all particles
    */
    MF k, L, wx, wy, wz;

    // dump accelearations
    for (particle& p : Particles) {
        p.wx_prev = p.wx;
        p.wy_prev = p.wy;
        p.wz_prev = p.wz;
        p.wx = 0;
        p.wy = 0;
        p.wz = 0;
    }

    for (int i = 0; i < particle_number; ++i) {
        for (int j = i + 1; j < particle_number; ++j) {
            // Get tuple {L, x, y, z}, unpuck it
            // (x, y, z) -- vector between 2 particles
            const auto [L, x, y, z] = nearest_reflection(Particles[i], Particles[j]);
            k = 24 * (2 * pow(L, -14) - pow(L, -8));
            wx = x * k;
            wy = y * k;
            wz = z * k;
            Particles[i].update_w( wx,  wy,  wz);
            Particles[j].update_w(-wx, -wy, -wz);
        }
    }
}

void model::update_coordinates(MF dt) {
    /*
        move particles
    */
    MF x, y, z;
    for (particle& p : Particles) {
        // make Verlet scheme step
        x = p.x; y = p.y; z = p.z;
        p.x = 2 * p.x - p.x_prev + p.wx * dt * dt;
        p.y = 2 * p.y - p.y_prev + p.wy * dt * dt;
        p.z = 2 * p.z - p.z_prev + p.wz * dt * dt;
        p.x_prev = x; p.y_prev = y; p.z_prev = z;

        // Calculate new velocity
        p.vx += 0.5 * dt * (p.wx_prev + p.wx);
        p.vy += 0.5 * dt * (p.wy_prev + p.wy);
        p.vz += 0.5 * dt * (p.wz_prev + p.wz);

        // handle borders
        if (p.x > size_of_box) p.x -= size_of_box;
        if (p.y > size_of_box) p.y -= size_of_box;
        if (p.z > size_of_box) p.z -= size_of_box;
        if (p.x < 0) p.x += size_of_box;
        if (p.y < 0) p.y += size_of_box;
        if (p.z < 0) p.z += size_of_box;

        // ToDo: diffusion
    }
}

void model::init_step(MF dt) {
    /*
        initial step when 
        there is no previous state
    */
    update_acceleration();

    for (particle& p : Particles) {
        p.x_prev = p.x;
        p.y_prev = p.y;
        p.z_prev = p.z;
        p.x += p.vx * dt + 0.5 * p.wx * dt * dt;
        p.y += p.vy * dt + 0.5 * p.wy * dt * dt;
        p.z += p.vz * dt + 0.5 * p.wz * dt * dt;
        p.vx += p.wx * dt;
        p.vy += p.wy * dt;
        p.vz += p.wz * dt;
    }
}

void model::make_step(MF dt) {
    update_acceleration();
    update_coordinates(dt);
}

void model::simulate(MF time, MF dt) {
    /*
        simulation loop
    */

    // save initial coordinates
    commit();

    // make first step
    init_step(dt);
    commit();

    // main algorithm
    for (size_t i = 2; dt * i < time; ++i) {
        make_step(dt);
        commit();
    }
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

void model::write() const {
    /*
        write history to a file
    */
    std::ofstream out("output.csv");

    out << header();

    for (const std::vector<MF>& v : history) {
        for (int i = 0; i < particle_number - 1; ++i) 
            out << v[i] << ',';
        out << v.back() << "\n";
    }

    out.close();
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