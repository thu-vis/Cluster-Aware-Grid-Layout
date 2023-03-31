#ifndef _MEASURE_DEVIATION_H
#define _MEASURE_DEVIATION_H

#include "../utils/convexHull.h"
#include "../utils/base.h"
#include <iostream>
#include <cmath>
#include <float.h>

// constants
namespace Deviation {
    const double epsilon = 0.01;
    const double alpha = 0.2;   // convexity = exp(-alpha*deviation)
    const double beta = 0.4;    // convexity = 1/(1+beta*deviation)
}

// 将 deviation 转化为 convexity measure 的策略
enum ConvertStrategy {
    EXPONENTIAL, // 指数函数 exp(-kx)
    INVERSE // 反比例函数 1/(1+kx)
};

// deviation 的计算方法
enum DeviationStrategy {
    SIMPLE_DEVIATION,   // sinple index of nonconvexity
    TOTAL_DEVIATION     // total index of nonconvexity
};

// 计算点 (x, y) 到线段 { (x1, y1), (x2, y2) } 的距离
double pointDist2Line(
const double &x, const double &y,
const double &x1, const double &y1,
const double &x2, const double &y2) {
    double p = std::abs((y2-y1)*x - (x2-x1)*y + x2*y1 - x1*y2);
    double q = std::sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1));
    return p / q;
}

// scale_rate: sqrt(多边形围成区域的面积)
// 引入 scale_rate 的原因：使凸性度量在放缩变换下保持不变
// (论文里似乎没有这么做)
double getSimpleDeviation(
const int grids[][2],
const int grids_n,
const double scale_rate) {
    double (*convex_hull)[2] = new double[grids_n][2];
    for (int i = 0; i < grids_n; i++) {
        convex_hull[i][0] = grids[i][0];
        convex_hull[i][1] = grids[i][1];
    }
    int vertex_num = getConvexHull(grids_n, convex_hull);

    // 找到 grids[] 里与 convex_hull[0] 相同的顶点
    int index = 0;
    for (int i = 0; i < grids_n; i++) {
        if ((double)grids[i][0] == convex_hull[0][0]
            && (double)grids[i][1] == convex_hull[0][1]) {
            index = i;
            break;
        }
    }

    double result = 0;
    for (int pos = 0; pos < vertex_num; pos++) {
        int new_pos = (pos + 1) % vertex_num;
        double x1 = convex_hull[pos][0];
        double y1 = convex_hull[pos][1];
        double x2 = convex_hull[new_pos][0];
        double y2 = convex_hull[new_pos][1];
        while (true) {
            index = (index + 1) % grids_n;
            if ((double)grids[index][0] == x2
            && (double)grids[index][1] == y2) {
                break;
            }
            // else
            double x = (double)grids[index][0];
            double y = (double)grids[index][1];
            double deviation = pointDist2Line(x, y, x1, y1, x2, y2);
            result = std::max(result, deviation);
        }
    }

    delete [] convex_hull;
    return result / scale_rate;
}

double getTotalDeviation(
const int grids[][2],
const int grids_n,
const double scale_rate) {
    double (*convex_hull)[2] = new double[grids_n][2];
    for (int i = 0; i < grids_n; i++) {
        convex_hull[i][0] = grids[i][0];
        convex_hull[i][1] = grids[i][1];
    }
    int vertex_num = getConvexHull(grids_n, convex_hull);

    int index = 0;
    for (int i = 0; i < grids_n; i++) {
        if ((double)grids[i][0] == convex_hull[0][0]
        && (double)grids[i][1] == convex_hull[0][1]) {
            index = i;
            break;
        }
    }

    double result = 0;
    for (int pos = 0; pos < vertex_num; pos++) {
        int new_pos = (pos + 1) % vertex_num;
        double x1 = convex_hull[pos][0];
        double y1 = convex_hull[pos][1];
        double x2 = convex_hull[new_pos][0];
        double y2 = convex_hull[new_pos][1];
        double deviation = 0;
        while (true) {
            index = (index + 1) % grids_n;
            if ((double)grids[index][0] == x2
            && (double)grids[index][1] == y2) {
                break;
            }
            // else
            double x = (double)grids[index][0];
            double y = (double)grids[index][1];
            deviation = std::max(deviation, pointDist2Line(x, y, x1, y1, x2, y2) / scale_rate);
        }
        result = (deviation > Deviation::epsilon) ? result + deviation : result;
    }

    delete [] convex_hull;
    return result;
}

inline double deviation2Convexity_1(
const double deviation) {
    return std::exp(-Deviation::alpha*deviation);
}

inline double deviation2Convexity_2(
const double deviation) {
    return 1 / (1 + Deviation::beta * deviation);
}

double convexityMeasure(
const int grids[][2],
const int grids_n, const double scale_rate,
const int convert=ConvertStrategy::EXPONENTIAL,
const int deviation=DeviationStrategy::TOTAL_DEVIATION) {
    if (convert == ConvertStrategy::EXPONENTIAL
    && deviation == DeviationStrategy::SIMPLE_DEVIATION) {
        return deviation2Convexity_1(getSimpleDeviation(grids, grids_n, scale_rate));
    } else if (convert == ConvertStrategy::EXPONENTIAL
    && deviation == DeviationStrategy::TOTAL_DEVIATION) {
        return deviation2Convexity_1(getTotalDeviation(grids, grids_n, scale_rate));
    } else if (convert == ConvertStrategy::INVERSE
    && deviation == DeviationStrategy::SIMPLE_DEVIATION) {
        return deviation2Convexity_2(getSimpleDeviation(grids, grids_n, scale_rate));
    }  else if (convert == ConvertStrategy::INVERSE
    && deviation == DeviationStrategy::TOTAL_DEVIATION) {
        return deviation2Convexity_2(getTotalDeviation(grids, grids_n, scale_rate));
    }
    else return 0;
}

std::vector<double> checkConvexForDevArray(
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

        tmp_S0 = convexityMeasure(nodes, cnt, std::sqrt(label_count[li])) * tmp_S1;
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

std::vector<double> checkConvexForDev(
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

    std::vector<double> ret = checkConvexForDevArray(grid_asses, cluster_labels, N, num, square_len, maxLabel);

    delete[] grid_asses;
    delete[] cluster_labels;

    return ret;
}

std::vector<double> checkCostForDev(
    const double Similar_cost_matrix[],
    const double Compact_cost_matrix[],
    const int grid_asses[], const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    const double &alpha, const double &beta) {

    std::vector<double> S_pair = checkConvexForDevArray(grid_asses, cluster_labels, N, num, square_len, maxLabel);
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