#ifndef _MEASURE_BOUNDARY_H
#define _MEASURE_BOUNDARY_H

#include "../utils/convexHull.h"
#include "../utils/base.h"
#include <cmath>
#include <iostream>
#include <set>
#include <float.h>

double measureOfRotate(
const int grids[][2],
const int grids_n,
const double angle) {
    double (*rotate_grids)[2] = new double[grids_n][2];
    double x_min = DBL_MAX;
    double y_min = DBL_MAX;
    double x_max = -DBL_MAX;
    double y_max = -DBL_MAX;
    for (int i = 0; i < grids_n; i++) {
        rotate_grids[i][0] = grids[i][0] * std::cos(angle) - grids[i][1] * std::sin(angle);
        rotate_grids[i][1] = grids[i][0] * std::sin(angle) + grids[i][1] * std::cos(angle);

        double x = rotate_grids[i][0];
        double y = rotate_grids[i][1];
        if (x < x_min) x_min = x;
        if (x > x_max) x_max = x;
        if (y < y_min) y_min = y;
        if (y > y_max) y_max = y;
    }
    // for (int i = 0; i < grids_n; i++) {
    //     std::cout << rotate_grids[i][0] << ", " << rotate_grids[i][1] << '\n';
    // }
    double per2 = 2 * (x_max - x_min + y_max - y_min);

    double per1 = 0;
    for (int i = 0; i < grids_n - 1; i++) {
        per1 += std::abs(rotate_grids[i][0] - rotate_grids[i+1][0]);
        per1 += std::abs(rotate_grids[i][1] - rotate_grids[i+1][1]);
    }
    per1 += std::abs(rotate_grids[0][0] - rotate_grids[grids_n-1][0]);
    per1 += std::abs(rotate_grids[0][1] - rotate_grids[grids_n-1][1]);
    delete [] rotate_grids;

    return per2 / per1;
}

double convexityMeasure(
const int grids[][2],
const int grids_n) {
    double (*convex_hull)[2] = new double[grids_n][2];
    for (int i = 0; i < grids_n; ++i) {
        convex_hull[i][0] = grids[i][0];
        convex_hull[i][1] = grids[i][1];
    }
    int vertex_num = getConvexHull(grids_n, convex_hull); // @vertex_num: num of vertices in the convex hull
    
    // compute all the rotation angles
    std::set<double> angles;
    for (int i = 0; i < vertex_num; i++) {
        double x1 = convex_hull[i][0];
        double y1 = convex_hull[i][1];
        size_t j = (i+1) % vertex_num;
        double x2 = convex_hull[j][0];
        double y2 = convex_hull[j][1];
        double angle = std::atan2(y2 - y1, x2 - x1);
        angles.insert(-angle);
        angles.insert(M_PI_2 - angle);
    }

    double measure = 1.0;
    for (std::set<double>::const_iterator iter = angles.begin(); iter != angles.end(); ++iter) {
        measure = std::min(measure, measureOfRotate(grids, grids_n, *iter));
        // std::cout << "angle: " << *(iter) / M_PI * 180 << '\n';
        // std::cout << "measure: " << measureOfRotate(grids, grids_n, *iter) << '\n';
    }

    delete [] convex_hull;
    return measure;
}

std::vector<double> checkConvexForBArray(
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    std::vector<std::vector<int>> grid_label(square_len);
    for (int i = 0; i < square_len; i++) {
        grid_label[i].resize(square_len);
    }
    int *head = new int[maxLabel];
    int *last = new int[num];
    int *element_asses = new int[num];
    int *label_count = new int[maxLabel];
    for (int i = 0; i < maxLabel; i++) {
        head[i] = -1;
        label_count[i] = 0;
    }
    for (int gid = 0; gid < N; gid++) {
        int x = gid % square_len;
        int y = gid / square_len;
        if (grid_asses[gid] < num) {
            int id = grid_asses[gid];
            int lb = cluster_labels[id];
            element_asses[id] = gid;
            last[id] = head[lb]; head[lb] = id;
            grid_label[x][y] = lb;
            label_count[lb]++;
        } else {
            grid_label[x][y] = -1;
        }
    }

    double S0 = 0, S1 = 0;
    int (*nodes)[2] = new int[num*4][2];
    for (int li = 0; li < maxLabel; li++) {
        //std::cout << "label: " << li << '\n';
        int cnt = 0;
        double tmp_S0 = 0, tmp_S1 = 0;
        int px = square_len; int py = square_len;
        for (int id = head[li]; id >= 0; id = last[id]) {
            int gid = element_asses[id];
            int x = gid % square_len;
            int y = gid / square_len;
            if ((y < py) || (y == py && x < px)) {
                px = x;
                py = y;
            }
        }
        nodes[0][0] = px; nodes[0][1] = py; cnt = 1;
        int nx = px+1; int ny = py; int dir = GO_RIGHT;
        // 规定走向为逆时针
        while (true) {
            //std::cout << "nx: " << nx << ", ny: " << ny << '\n';
            if (nx == px && ny == py) break;
            switch (dir) {
                case GO_UP: {
                    int label = getLabel(grid_label, nx, ny);
                    if (label == li) {
                        nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                        nx++;
                        dir = GO_RIGHT;
                    } else {
                        label = getLabel(grid_label, nx-1, ny);
                        if (label != li) {
                            nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                            nx--;
                            dir = GO_LEFT;
                        } else {
                            ny++;
                            dir = GO_UP;
                        }
                    }
                    break;
                }
                case GO_DOWN: {
                    int label = getLabel(grid_label, nx-1, ny-1);
                    if (label == li) {
                        nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                        nx--;
                        dir = GO_LEFT;
                    } else {
                        label = getLabel(grid_label, nx, ny-1);
                        if (label != li) {
                            nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                            nx++;
                            dir = GO_RIGHT;
                        } else {
                            ny--;
                            dir = GO_DOWN;
                        }
                    }
                    break;
                }
                case GO_LEFT: {
                    int label = getLabel(grid_label, nx-1, ny);
                    if (label == li) {
                        nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                        ny++;
                        dir = GO_UP;
                    } else {
                        label = getLabel(grid_label, nx-1, ny-1);
                        if (label != li) {
                            nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                            ny--;
                            dir = GO_DOWN;
                        } else {
                            nx--;
                            dir = GO_LEFT;
                        }
                    }
                    break;
                }
                case GO_RIGHT: {
                    int label = getLabel(grid_label, nx, ny-1);
                    if (label == li) {
                        nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                        ny--;
                        dir = GO_DOWN;
                    } else {
                        label = getLabel(grid_label, nx, ny);
                        if (label != li) {
                            nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                            ny++;
                            dir = GO_UP;
                        } else {
                            nx++;
                            dir = GO_RIGHT;
                        }
                    }
                    break;
                }
            }
        }
        tmp_S1 = label_count[li];

        // debug
        // std::cout << "polygon:\n";
        // for (int i = 0; i < cnt; i++) {
        //     std::cout << nodes[i][0] << ", " << nodes[i][1] << '\n';
        // }

        tmp_S0 = convexityMeasure(nodes, cnt) * tmp_S1;
        // std::cout << "tmp_S0: " << tmp_S0 << ", tmp_S1: " << tmp_S1 << '\n';

        S1 += tmp_S1;
        S0 += tmp_S0;
    }
    std::vector<double> S_pair(2, 0);
    S_pair[0] = S1 - S0;
    S_pair[1] = S1;

    delete [] head;
    delete [] last;
    delete [] element_asses;
    delete [] nodes;
    delete [] label_count;
    return S_pair;
}

std::vector<double> checkConvexForB(
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels) {
    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = std::max(maxLabel, _cluster_labels[i]+1);

    int *grid_asses = new int[N];
    int *cluster_labels = new int[num];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];

    std::vector<double> ret = checkConvexForBArray(grid_asses, cluster_labels, N, num, square_len, maxLabel);

    delete[] grid_asses;
    delete[] cluster_labels;
    
    return ret;
}

std::vector<double> checkCostForB(
    const double Similar_cost_matrix[],
    const double Compact_cost_matrix[],
    const int grid_asses[], const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    const double &alpha, const double &beta) {

    std::vector<double> S_pair = checkConvexForBArray(grid_asses, cluster_labels, N, num, square_len, maxLabel);
    double correct=0, full=0;
    full = S_pair[1];
    correct = S_pair[1]-S_pair[0];

    double Convex_cost = (full-correct)/full;
    double Similar_cost = 0;
    double Compact_cost = 0;
    double cost = 0;

    int *element_asses = new int[num];
    for(int i=0;i<N;i++)if(grid_asses[i]<num)element_asses[grid_asses[i]] = i;

    for(int i=0;i<num;i++){
        Similar_cost += Similar_cost_matrix[element_asses[i]*N+i];
        Compact_cost += Compact_cost_matrix[element_asses[i]*N+i];
        cost += Similar_cost_matrix[element_asses[i]*N+i] * (1-beta-alpha);
        cost += Compact_cost_matrix[element_asses[i]*N+i] * beta;
    }
    cost += Convex_cost * N * alpha;

//    printf("cost %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n", Similar_cost, Compact_cost, Convex_cost*N, beta, alpha, cost);

    delete[] element_asses;
    std::vector<double> ret(4, 0);
    ret[0] = cost;
    ret[1] = Similar_cost;
    ret[2] = Compact_cost;
    ret[3] = Convex_cost*N;
    return ret;
}

#endif