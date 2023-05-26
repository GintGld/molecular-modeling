#pragma once
#include "MF.h"
#include <tuple>

struct particle {
    MF 
        x, y, z,
        vx, vy, vz,
        wx, wy, wz,
        step_x, step_y, step_z;

    particle();
    particle(MF x, MF y, MF z, MF vx, MF vy, MF vz);

    void update_w(MF wx, MF wy, MF wz);

    static MF dist(const particle&, const particle&);
    static MF dist(const particle&, std::tuple<MF, MF, MF>);
};