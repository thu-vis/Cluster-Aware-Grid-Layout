#ifndef _MEASURE_ALPHA_TRIPLES_H
#define _MEASURE_ALPHA_TRIPLES_H

#include "../utils/base.h"
#include "../utils/geometry.h"

static std::vector<double> grid_pos(int point){
    int x = point / grid_square;
    int y = point % grid_square;
    return std::vector<double> {1.0 * x, 1.0 * y};
}

std::vector<double> checkConvexForAlphaTriplesArray(
//    const std::vector<int> &grid_asses,
//    const std::vector<int> &cluster_labels,
    const int grid_asses[], const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    double alpha = 0.33)
{
    // get position of each grid cell and calculate
    grid_square = square_len;
    std::vector<std::vector<double>> points;
    for(auto i = 0; i < N;i ++) {
        points.push_back(grid_pos(i));
    }
    int max_label = maxLabel;
    std::vector<std::vector<double>> label_pairs(max_label, std::vector<double>(2, 0));
    for(auto i = 0; i < N; i++) {
        for(auto j = 0; j < N; j++) {
            if(i == j) continue;
            if(grid_asses[i]>=num)continue;
            if(grid_asses[j]>=num)continue;
            if(cluster_labels[grid_asses[i]] != cluster_labels[grid_asses[j]]) continue;
            auto label = cluster_labels[grid_asses[i]];
            auto & pos1 = points[i];
            auto & pos2 = points[j];
            int alpha_grid = grid(round(pos1[0] * alpha + pos2[0] *(1 - alpha)), round(pos1[1] * alpha + pos2[1] * (1 - alpha)));
            if((grid_asses[alpha_grid]<num)&&(cluster_labels[grid_asses[alpha_grid]] == label)) {
                label_pairs[label][0] += 1;
            }
            label_pairs[label][1] += 1;
        }
    }
    double x = 0, y = 0;
    for(auto & pair: label_pairs) {
        x += pair[0];
        y += pair[1];
    }
    return {y-x, y};
}

std::vector<double> checkCostForAlphaT(
    const double Similar_cost_matrix[],
    const double Compact_cost_matrix[],
    const int grid_asses[], const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    const double &alpha, const double &beta) {

    std::vector<double> S_pair = checkConvexForAlphaTriplesArray(grid_asses, cluster_labels, N, num, square_len, maxLabel);
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