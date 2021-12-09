cmake_minimum_required(VERSION 3.10.0)
project(SphericalOMT)

# == Deps
find_package(CGAL REQUIRED)# this seems to must be before add_executable to compile
add_subdirectory(deps/polyscope)

# == Srcs
set(SRCS src/main.cpp
)
add_executable(SphericalOMT "${SRCS}")

# == link & include
target_link_libraries(SphericalOMT PUBLIC CGAL::CGAL PUBLIC polyscope)
target_include_directories(SphericalOMT PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(SphericalOMT PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/Eigen/src")