#ifndef _BASE_H
#define _BASE_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <utility>
#include <ctime>

#define THREADS_NUM 48

#define M_PI 3.1415926535898
#define M_PI_2 3.1415926535898/2

namespace py = pybind11;
using std::cout;

enum Direction {
    GO_UP,
    GO_DOWN,
    GO_LEFT,
    GO_RIGHT
};

inline double sqr(double x) {
    return x * x;
}

template <typename T>
inline T min(const T&a, const T& b) {
    return a < b ? a : b;
}

template <typename T>
inline T max(const T&a, const T& b) {
    return a > b ? a : b;
}

double getDist(double x1, double y1, double x2, double y2){
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

int getLabel(
const std::vector<std::vector<int>> &grid_label,
const int &x, const int &y) {
    int n = grid_label.size();
    if (x < 0 || y < 0 || x >= n || y >= n) return -1;
    return grid_label[x][y];
}

#endif