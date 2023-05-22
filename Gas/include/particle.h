#pragma once
#include "MF.h"
#include <tuple>

struct particle {
    MF 
        x, y, z,
        //x_prev, y_prev, z_prev, deprecated
        vx, vy, vz,
        wx, wy, wz;
        //wx_prev, wy_prev, wz_prev; deprecated

    particle();
    particle(MF x, MF y, MF z, MF vx, MF vy, MF vz);

    void update_w(MF wx, MF wy, MF wz);

    static MF dist(const particle&, const particle&);
    static MF dist(const particle&, std::tuple<MF, MF, MF>);
};