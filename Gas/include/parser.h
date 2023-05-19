#pragma once

#include <filesystem>

#include "model.h"

namespace fs=std::filesystem;

class parser {
private:
    fs::path file;
    int particle_number;

public:
    std::vector<particle> Particles;

public:
    parser();
    parser(const fs::path&);
    void read(const fs::path&);
    model get_model() const;
};