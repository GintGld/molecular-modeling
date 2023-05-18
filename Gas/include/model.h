#pragma once
#include "particle.h"
#include <vector>
#include <string>

using std::size_t;

class model {
private:
    std::vector<particle> Particles;
    std::vector< std::vector<MF> >history;
    std::vector<MF> hist_tmp;
    size_t particle_number;

public:
    model();

    void add_particle(const particle&);
    void correct_CM();

private:
    void update_acceleration();
    void update_coordinates(MF dt);

    static std::tuple<MF, MF, MF, MF> // L, x, y, z
    nearest_reflection(const particle&, const particle&);

    void init_step(MF dt);
    void make_step(MF dt);

    void save_point(const particle&);
    void commit();

    std::string header() const;

public:
    void simulate(MF time, MF dt);
    void write() const;

    const std::vector<particle>& get_particles() const;
    std::vector<particle>& get_particles();
    const std::vector< std::vector<MF> >& get_history();

    static MF size_of_box, relaxation_time, record_time;
};