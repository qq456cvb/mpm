#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <armadillo>
using namespace arma;
using namespace std;

class Particle {
public:
    vec2 _x;
    vec2 _v;
    mat22 _c;
    float _mass = 1;
};

class Cell {
public:
    vec2 _v;
    float _mass = 1;
};

class Simulator {
    Cell **_cells;
    Particle *_particles;
    float _grid_spacing;
    int _grid_size;
    int _num_particles;
    float _dt = 1.f, _gravity = -0.05f;

public:
    Simulator(int grid_size, int num_particles, float grid_spacing=1.f);

    void step();
};

#endif