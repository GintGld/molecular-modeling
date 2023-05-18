include_directories("include")

add_library(Gas OBJECT
    "include/MF.h"
    "include/particle.h"
    "include/model.h"

    "src/particle.cpp"
    "src/model.cpp"
)