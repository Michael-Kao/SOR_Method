#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <utility>

struct SOR_Method {
    SOR_Method(double _w, double _e, int _gridSize, std::pair<int, int> _domain, double BC_V, double SRC_V, double BCD_V);
    void init();
    void set_Dirichlet_BC();
    void set_Neuman_BC();
    void set_para(double _w, double _e, int _gridSize);
    double get_source(int i, int j);
    double SOR_iter();
    void run();

    double epsilon;
    int grid_size;
    std::pair<int, int> domain;
    double w;
    int iter_num;
    double BC_Derivative, BC_value, SRC_value;
    std::vector< std::vector<double> > T;
};