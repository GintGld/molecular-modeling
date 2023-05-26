#include <cmath>
#include "particle.h"

particle::particle() {}

particle::particle(MF x, MF y, MF z, MF vx, MF vy, MF vz):
    x(x), y(y), z(z),
    vx(vx), vy(vy), vz(vz),
    wx(0), wy(0), wz(0),
    step_x(0), step_y(0), step_z(0) {}

void particle::update_w(MF wx_, MF wy_, MF wz_) {
    wx += wx_;
    wy += wy_;
    wz += wz_;
}

MF particle::dist(const particle& p1, const particle& p2) {
    return sqrt(
        (p1.x - p2.x) * (p1.x - p2.x) +
        (p1.y - p2.y) * (p1.y - p2.y) +
        (p1.z - p2.z) * (p1.z - p2.z)
    );
}

MF particle::dist(const particle& p, std::tuple<MF, MF, MF> t) {
    const auto [x, y, z] = t;
    return sqrt(
        (p.x - x) * (p.x - x) +
        (p.y - y) * (p.y - y) +
        (p.z - z) * (p.z - z)
    );
}