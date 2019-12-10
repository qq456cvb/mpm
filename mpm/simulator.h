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
    mat22 _C;
    mat22 _F;
    float _mass = 0;
    float _vol0 = 0;
};

class Cell {
public:
    vec2 _v;
    float _mass = 0;
};

class Simulator {
    Cell **_cells;
    Particle *_particles;
    float _grid_spacing;
    int _grid_size;
    int _num_particles;
    float _dt = 0.1f, _gravity = -0.3f;
    const float _lambda = 10.f;
    const float _mu = 20.f;
    int _cnt = 0;

public:
    Simulator(int grid_size, int num_particles, float grid_spacing=1.f, float particle_spacing=1.f);

    void step();
    void clearG();
    void P2G();
    void updateG();
    void G2P();
    void render(Mat<unsigned char>&);
};

#endif