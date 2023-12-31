#include "SOR.h"

SOR_Method::SOR_Method(double _w, double _e, int _gridSize, std::pair<int, int> _domain, double BC_V, double SRC_V, double BCD_V) 
    : w(_w), epsilon(_e), grid_size(_gridSize), domain(_domain), BC_value(BC_V), SRC_value(SRC_V), BC_Derivative(BCD_V) {
    init();
}

void SOR_Method::init() {
    T = std::vector< std::vector<double> >(grid_size, std::vector<double>(grid_size, 0));
    set_Dirichlet_BC();
}

void SOR_Method::set_para(double _w, double _e, int _gridSize) {
    w = _w;
    epsilon = _e;
    grid_size = _gridSize;
    init();
}

void SOR_Method::set_Dirichlet_BC() {
    for(int i = 0; i < grid_size; ++i) {
        T[i][0] = BC_value; // bottom
        T[0][i] = BC_value; // left
        T[grid_size - 1][i] = BC_value; // right
    }
}

void SOR_Method::set_Neuman_BC() {
    double h = (1.0 * domain.first) / (grid_size - 1);
    for(int i = 0; i < grid_size; ++i) {
        T[i][grid_size - 1] = T[i][grid_size - 1] + w * (T[i][grid_size - 2] + h * BC_Derivative - T[i][grid_size - 1]);
    }
}

double SOR_Method::get_source(int i, int j) {
    if(i == (grid_size - 1) / 2 && j == i)
        return SRC_value;
    else
        return (0.0);
}

double SOR_Method::SOR_iter() {
    double err = 0.0;
    double h = (1.0 * domain.first) / (grid_size - 1);
    for(int i = 1; i < grid_size - 1; ++i) {
        for(int j = 1; j < grid_size - 1; ++j) {
            double old = T[i][j];
            double s = get_source(i, j);
            T[i][j] = T[i][j] + (w / 4.0) * (-s * h * h + (T[i - 1][j] + T[i][j - 1] + T[i + 1][j] + T[i][j + 1]) - 4 * T[i][j]);
            double temp = fabs(T[i][j] - old);
            err = (temp > err) ? temp : err;
        }
    }

    set_Neuman_BC();
    return err;
}

void SOR_Method::run() {
    iter_num = 0;
    double err = 999;
    while(err > epsilon) {
        err = SOR_iter();
        iter_num += 1;
    }
    std::cout << "total iteration = " << iter_num << std::endl;
    // std::cout << std::setprecision(7);
    // for(int i = 0; i < grid_size; ++i) {
    //     for(int j = 0; j < grid_size; ++j) {
    //         std::cout << T[i][j];
    //         if(j != grid_size - 1)
    //             std::cout << std::setw(10);
    //     }
    //     std::cout << std::endl;
    // }
}