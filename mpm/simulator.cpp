#include "simulator.h"

Simulator::Simulator(int grid_size, int num_particles, float grid_spacing, float particle_spacing)
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
            Particle &p = _particles[j * num_x + i];
            p._x = { grid_size / 2. + (i - num_x / 2.) * particle_spacing,
                grid_size / 2. + (j - num_y / 2.) * particle_spacing };

            //_particles[j * num_x + i]._v = { double(rand()) / RAND_MAX - 0.5,  double(rand()) / RAND_MAX + 3.5 };
            // _particles[j * num_x + i]._v *= 0.5;
            p._v = { 0, 0 };
            p._F.eye();
            
            
            p._C.zeros();
            p._mass = 1.;
        }
    }

    P2G();
    for (size_t i = 0; i < _num_particles; i++)
    {
        Particle &p = _particles[i];

        ivec2 center_cell_idx = conv_to<ivec>::from(p._x);
        vec2 t = p._x - center_cell_idx - 0.5;  // 0.5 offset between grid center and index

        vec2 weights[3] = { 0.5 * square(0.5 - t), 0.75 - square(t), 0.5 * square(0.5 + t) };
        double density = 0;
        for (int gx = 0; gx < 3; gx++)
        {
            for (int gy = 0; gy < 3; gy++)
            {
                auto weight = weights[gx][0] * weights[gy][1];

                ivec2 cell_idx = center_cell_idx + ivec2({ gx - 1, gy - 1 });
                density += _cells[cell_idx[0]][cell_idx[1]]._mass * weight;
                /*if (i == 0) {
                    cout << gx << ", " << gy << ", " << weight << ", " << _cells[cell_idx[0]][cell_idx[1]]._mass << endl;
                    fflush(stdout);
                }*/
            }
        }
        //cout << density << endl;
        p._vol0 = p._mass / density;
    }
}

void Simulator::clearG() {
    for (size_t i = 0; i < _grid_size; i++)
    {
        for (size_t j = 0; j < _grid_size; j++)
        {
            _cells[i][j]._v.zeros();
            _cells[i][j]._mass = 0;
        }
    }
}

void Simulator::P2G() {
    // P2G
    for (size_t i = 0; i < _num_particles; i++)
    {
        const Particle &p = _particles[i];

        double J = det(p._F);
        double vol = J * p._vol0;

        mat22 F_inv_t = inv(p._F.t());
        mat22 P = _mu * (p._F - F_inv_t) + _lambda * log(J) * F_inv_t;
        mat22 stress = 1. / J * P * p._F.t();
        mat22 force_update = -vol * 4 * stress * _dt;

        ivec2 center_cell_idx = conv_to<ivec>::from(p._x);
        vec2 t = p._x - center_cell_idx - 0.5;  // 0.5 offset between grid center and index

        vec2 weights[3] = { 0.5 * square(0.5 - t), 0.75 - square(t), 0.5 * square(0.5 + t) };

        for (int gx = 0; gx < 3; gx++)
        {
            for (int gy = 0; gy < 3; gy++)
            {
                auto weight = weights[gx][0] * weights[gy][1];

                ivec2 cell_idx = center_cell_idx + ivec2({ gx - 1, gy - 1 });
                auto mass_contrib = weight * p._mass;
                _cells[cell_idx[0]][cell_idx[1]]._mass += mass_contrib;

                /*if (i == 0) {
                    cout << gx << ", " << gy << ", " << weight << ", " << p._mass << endl;
                    fflush(stdout);
                }*/

                vec2 delta_x = cell_idx - p._x + 0.5;
                _cells[cell_idx[0]][cell_idx[1]]._v += mass_contrib * (p._v + p._C * delta_x);  // momentum

                // force is doing work
                _cells[cell_idx[0]][cell_idx[1]]._v += force_update * weight * delta_x;  // momentum

            }
        }
    }
}

void Simulator::updateG() {
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
}

void Simulator::G2P() {
    // G2P
    for (size_t i = 0; i < _num_particles; i++)
    {
        Particle &p = _particles[i];

        p._v.zeros();

        ivec2 center_cell_idx = conv_to<ivec>::from(p._x);
        if (i == 0 && _cnt == 0) {
            cout << center_cell_idx << endl;
            cout << p._x << endl;
        }
        vec2 t = p._x - center_cell_idx - 0.5;  // 0.5 offset between grid center and index

        vec2 weights[3] = { 0.5 * square(0.5 - t), 0.75 - square(t), 0.5 * square(0.5 + t) };

        mat22 B;
        B.zeros();
        for (int gx = 0; gx < 3; gx++) {
            for (int gy = 0; gy < 3; gy++) {
                auto weight = weights[gx][0] * weights[gy][1];

                ivec2 cell_idx = center_cell_idx + ivec2({ gx - 1, gy - 1 });

                vec2 delta_x = cell_idx - p._x + 0.5;
                vec2 weighted_velocity = _cells[cell_idx[0]][cell_idx[1]]._v * weight;

                B += weighted_velocity * delta_x.t();
                /*if (_cnt == 0 && i == 0)
                {
                    cout << _cells[cell_idx[0]][cell_idx[1]]._v << ", " << weight << ", " << delta_x << endl;
                }*/

                p._v += weighted_velocity;
            }
        }
        p._C = B * 4;

        // advect particles
        p._x += p._v * _dt;

        // safety clamp to ensure particles don't exit simulation domain
        p._x = clamp(p._x, 1., _grid_size - 2.);

        // deformation gradient update
        p._F += _dt * p._C * p._F;
    }
}

void Simulator::step() {
    clearG();
    P2G();
    updateG();
    G2P();
    _cnt++;
}



void Simulator::render(Mat<unsigned char> &image) {
    double scale = min(image.n_rows / 3, image.n_cols) / double(_grid_size);
    //int test = 10000;
    for (size_t i = 0; i < _num_particles; i++) 
    {
        auto &p = _particles[i];
        int x = int(p._x[0] * scale), y = int(p._x[1] * scale);
        image(x * 3, y) = 255;
        image(x * 3 + 1, y) = 255;
        image(x * 3 + 2, y) = 255;
        //test = min(test, int(p._x[1] * scale));
    }
    //printf("test %d\n", test);
}