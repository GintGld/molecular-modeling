cmake_minimum_required(VERSION 3.14)
project(model)
add_executable("${PROJECT_NAME}" main.cpp)
add_executable("generator" generator.cpp)

add_subdirectory(Gas)

target_link_libraries("${PROJECT_NAME}" PUBLIC Gas)
target_include_directories("${PROJECT_NAME}" PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        "${PROJECT_SOURCE_DIR}/src")

# install binary
install(TARGETS "${PROJECT_NAME}" DESTINATION "./bin")
install(TARGETS "generator" DESTINATION "./bin")