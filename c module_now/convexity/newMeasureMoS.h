#ifndef _MEASURE_MOS_H
#define _MEASURE_MOS_H

/************************************
 *
 * Author: Jiashu
 * Date: 2023-03-12
 * Description: Convexity measure of MoS from paper: Mesh Segmentation via Spectral Embedding and Contour Analysis
 *      MoS = Measure of Segmentation
 *
 * ***********************************/

#include "measureCHC.h"
#include "measureCHS.h"
#include "../utils/geometry.h"

// polygon form
std::vector<double> checkPolygonConvexByMoS(
    const PointList & polygon,
    const std::vector<PointList> & holes,
    double alpha = 4)
{
    // get polygon convex hull
    double (*polygon_cp)[2] = new double [polygon.size()][2];
    for(auto i = 0;i < polygon.size(); i++) {
        polygon_cp[i][0] = polygon[i][0];
        polygon_cp[i][1] = polygon[i][1];
    }
    int num_nodes = getConvexHull(polygon.size(), polygon_cp);
    PointList convexHull;
    for(auto i = 0;i < num_nodes;i++){
        convexHull.push_back({polygon_cp[i][0], polygon_cp[i][1]});
    }
    delete[] polygon_cp;

    // caculate convexhull perimeter and area
    double area_hull = getArea(convexHull);
    double perimeter_hull = getPerimeter(convexHull);
    double area = getArea(polygon);
    double perimeter = getPerimeter(polygon);
    for(auto & hole: holes) {
        area -= getArea(hole);
        perimeter += getPerimeter(hole);
    }
    double C = max(area / area_hull, pow(perimeter_hull / perimeter, alpha));
    return {C, 1};
}


// grid form based on measure CHC and CHS
std::vector<double> checkConvexForMoS(
    const std::vector<int> &_grid_asses,
    const std::vector<int> &_cluster_labels,
    double alpha = 4)
{
    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = max(maxLabel, _cluster_labels[i]+1);

    // check convexhull C
    int *grid_asses = new int[N];
    int *cluster_labels = new int[num];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];
    double (*res_c)[2] = new double [maxLabel][2];
    checkConvexForCArray(grid_asses, cluster_labels, N, num, square_len, maxLabel, true, false, res_c);

    // check convexhull S
    int *grid_asses2 = new int[N];
    int *cluster_labels2 = new int[num];
    for(int i=0;i<N;i++)grid_asses2[i] = _grid_asses[i];
    for(int i=0;i<num;i++)cluster_labels2[i] = _cluster_labels[i];
    double (*res_s)[2] = new double [maxLabel][2];
    checkConvexForSArray(grid_asses2, cluster_labels2, N, num, square_len, maxLabel, true, false, res_s);

    // caculate C = max(C1, C2_ALPHA)
    std::vector<std::vector<double>> C;
    double total_perimeter = 0, total_area = 0;
    for(int i=0;i<maxLabel;i++){
        double C1 = res_c[i][0] / res_c[i][1];
        double C2 = res_s[i][0] / res_s[i][1];
        C.push_back({max(C1, pow(C2, alpha)), res_c[i][1], res_s[i][0]});
        total_perimeter += res_c[i][1];
        total_area += res_s[i][0];
    }

    delete[] grid_asses;
    delete[] cluster_labels;
    delete[] res_c;
    delete[] grid_asses2;
    delete[] cluster_labels2;
    delete[] res_s;

    // summary: to be discussed
    double summary_C = 0;
    for(int i=0;i<maxLabel;i++){
        summary_C += C[i][0] * (C[i][1] / total_perimeter + C[i][2] / total_area) / 2;
    }
    return {1 - summary_C, 1};
}

std::vector<double> checkConvexForMoSArray(
    const int grid_asses[],
    const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel)
{
    return checkConvexForMoS(std::vector<int>(grid_asses, grid_asses+N), std::vector<int>(cluster_labels, cluster_labels+num));
}

std::vector<double> checkCostForMoS(
    const double Similar_cost_matrix[],
    const double Compact_cost_matrix[],
    const int grid_asses[], const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    const double &alpha, const double &beta) {

    std::vector<double> MoS_pair = checkConvexForMoSArray(grid_asses, cluster_labels, N, num, square_len, maxLabel);
    double correct=0, full=0;
    full = MoS_pair[1];
    correct = MoS_pair[1]-MoS_pair[0];

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

    // printf("cost %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n", Similar_cost, Compact_cost, Convex_cost*N, beta, alpha, cost);

    delete[] element_asses;

    std::vector<double> ret(4, 0);
    ret[0] = cost;
    ret[1] = Similar_cost;
    ret[2] = Compact_cost;
    ret[3] = Convex_cost*N;
    return ret;
}


#endif