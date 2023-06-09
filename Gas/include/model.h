#pragma once

#include "particle.h"

#include <vector>
#include <string>
#include <filesystem>

using std::size_t;
namespace fs=std::filesystem;

class parser;

class model {
private:
    std::vector<particle> Particles;
    std::vector< std::vector<MF> >history;
    std::vector<MF> kinetic_energy, potential_energy;
    MF kinetic_energy_tmp, potential_energy_tmp;
    std::vector<MF> hist_tmp;
    size_t particle_number;

public:
    model();

    void add_particle(const particle&);
    void correct_CM();

private:
    void update_acceleration();
    void update_coordinates(MF dt);

    static std::tuple<MF, MF, MF> // x, y, z
    nearest_reflection(const particle&, const particle&);

    void init_step(MF dt);
    void make_step(MF dt);

    void save_point(const particle&);
    void commit();
    void scale_particles(MF target_temp, MF dt = 0);

    std::string header() const;

public:
    void simulate(MF time, MF dt, MF target_temp = 1);
    void write(const fs::path&) const;
    void write_last_state(const fs::path&, const std::string& = " ") const;
    void write_energy(const fs::path&, const std::string& = ",") const;
    void write_ovito(const fs::path&, const int = 1) const;

    const std::vector<particle>& get_particles() const;
    std::vector<particle>& get_particles();
    const std::vector< std::vector<MF> >& get_history();

    static MF size_of_box;
    static bool without_centering_CM;
    static bool scale, diffusion;
    static MF relaxation_time;
};