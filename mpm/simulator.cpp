#include "simulator.h"

Simulator::Simulator(int grid_size, int num_particles, float grid_spacing) 
    : _grid_size(grid_size), _num_particles(num_particles), _grid_spacing(grid_spacing) {
    srand(time(NULL));

    _cells = new Cell *[grid_size];
    _cells[0] = new Cell[grid_size * grid_size];
    for (size_t i = 1; i < grid_size; i++) {
        _cells[i] = &_cells[0][i * grid_size];
    }

    _particles = new Particle[num_particles];

    // init particles
    int num_x = int(sqrt(num_particles));
    int num_y = num_particles / num_x;
    
    for (size_t j = 0; j < num_y; j++) {
        for (size_t i = 0; i < num_x; i++) {
            _particles[j * num_x + i]._x = { grid_size / 2. - num_x / 2 + i, grid_size / 2. - num_y / 2 + j };

            //_particles[j * num_x + i]._v = { 0, 0 };
             _particles[j * num_x + i]._v = { double(rand()) / RAND_MAX - 0.5,  double(rand()) / RAND_MAX + 3.5 };
             _particles[j * num_x + i]._v *= 0.5;
            
            _particles[j * num_x + i]._c.zeros();
            _particles[j * num_x + i]._mass = 1.;
        }
    }
}

void Simulator::step() {
    for (size_t i = 0; i < _grid_size; i++)
    {
        for (size_t j = 0; j < _grid_size; j++)
        {
            _cells[i][j]._v.zeros();
            _cells[i][j]._mass = 0;
        }
    }

    // P2G
    for (size_t i = 0; i < _num_particles; i++)
    {
        const auto &p = _particles[i];
        ivec2 center_cell_idx = conv_to<ivec>::from(p._x);
        auto t = p._x - center_cell_idx - 0.5;  // 0.5 offset between grid center and index

        vec2 weights[3] = { 0.5 * square(0.5 - t), 0.75 - square(t), 0.5 * square(0.5 + t) };

        for (int gx = 0; gx < 3; gx++)
        {
            for (int gy = 0; gy < 3; gy++)
            {
                auto weight = weights[gx][0] * weights[gy][1];

                ivec2 cell_idx = center_cell_idx + ivec2({ gx - 1, gy - 1 });
                auto mass_contrib = weight * p._mass;
                _cells[cell_idx[0]][cell_idx[1]]._mass += mass_contrib;

                vec2 delta_x = center_cell_idx - p._x + 0.5;
                _cells[cell_idx[0]][cell_idx[1]]._v += mass_contrib * (p._v + p._c * delta_x);  // momentum
            }
        }
    }

    // Grid Dynamics
    for (size_t i = 0; i < _grid_size; i++)
    {
        for (size_t j = 0; j < _grid_size; j++)
        {
            auto &cell = _cells[i][j];
            if (cell._mass > 0)
            {
                cell._v /= cell._mass;
                cell._v += _dt * vec2({ 0, _gravity });

                if (i < 2 || i > _grid_size - 3)
                {
                    cell._v[0] = 0;
                }
                if (j < 2 || j > _grid_size - 3)
                {
                    cell._v[1] = 0;
                }
            }
        }
    }

    // G2P
    for (size_t i = 0; i < _num_particles; i++)
    {
        auto &p = _particles[i];

        p._v.zeros();

        ivec2 center_cell_idx = conv_to<ivec>::from(p._x);
        auto t = p._x - center_cell_idx - 0.5;  // 0.5 offset between grid center and index

        vec2 weights[3] = { 0.5 * square(0.5 - t), 0.75 - square(t), 0.5 * square(0.5 + t) };

        mat22 B;
        B.zeros();
        for (int gx = 0; gx < 3; gx++) {
            for (int gy = 0; gy < 3; gy++) {
                auto weight = weights[gx][0] * weights[gy][1];

                ivec2 cell_idx = center_cell_idx + ivec2({ gx - 1, gy - 1 });

                vec2 delta_x = center_cell_idx - p._x + 0.5;
                vec2 weighted_velocity = _cells[cell_idx[0]][cell_idx[1]]._v * weight;

                B += weighted_velocity * delta_x.t();

                p._v += weighted_velocity;
            }
        }
        p._c = B * 4;

        // advect particles
        p._x += p._v * _dt;

        // safety clamp to ensure particles don't exit simulation domain
        p._x = clamp(p._x, 1., _grid_size - 2.);
    }
}



void Simulator::render(Mat<unsigned char> &image) {
    double scale = min(image.n_rows / 3, image.n_cols) / double(_grid_size);
    int test = 10000;
    for (size_t i = 0; i < _num_particles; i++) 
    {
        const auto &p = _particles[i];
        image(int(p._x[0] * scale) * 3, int(p._x[1] * scale)) = 255;
        image(int(p._x[0] * scale) * 3 + 1, int(p._x[1] * scale)) = 255;
        image(int(p._x[0] * scale) * 3 + 2, int(p._x[1] * scale)) = 255;
        test = min(test, int(p._x[1] * scale));
    }
    cout << test << endl;
}