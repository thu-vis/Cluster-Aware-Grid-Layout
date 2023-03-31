#ifndef _MEASURE_CHS_H
#define _MEASURE_CHS_H

#include "../utils/base.h"
#include "../utils/convexHull.h"
#include "../utils/util.h"

//calculate the convexity based on area(where S means Square) of convex hull
//save: if save convexity measure pairs of every cluster
//load: if load saved convexity measure pairs of every cluster and only re-calculate clusters of "mainLabel1" and "mainLabel2"
//mainLabel1, mainLabel2: re-calculate clusters of "mainLabel1" and "mainLabel2" when load saved convexity measure pairs
//return S[0..1]: convexity(in the form of loss/cost) pair
std::vector<double> checkConvexForSArray(
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
bool save=false, bool load=false, double label_pairs[][2]=nullptr, int mainLabel1=-1, int mainLabel2=-1) {

    int *head = new int[maxLabel];
    int *last = new int[num];
    int *element_asses = new int[num];
    int (*E_grid)[4] = new int[N][4];

    checkEdgeArray4(E_grid, grid_asses, cluster_labels, N, num, square_len, maxLabel);

    for(int i=0;i<maxLabel;i++)head[i] = -1;

    for(int gid=0;gid<N;gid++)
    if(grid_asses[gid]<num) {
        int id = grid_asses[gid];
        element_asses[id] = gid;
        int lb = cluster_labels[id];
        last[id] = head[lb]; head[lb] = id;
    }

    double (*nodes)[2] = new double[num*4][2];
    double S0 = 0, S1 = 0;
    for(int li=0;li<maxLabel;li++){
        if((load)&&(li!=mainLabel1)&&(li!=mainLabel2)){
            S0 += label_pairs[li][0];
            S1 += label_pairs[li][1];
            continue;
        }
        int cnt=0;

        double tmp_S0=0, tmp_S1=0;
        for(int id=head[li];id>=0;id=last[id]){
            int gid = element_asses[id];
            int x = gid/square_len;
            int y = gid%square_len;

            if((E_grid[gid][0]==1)&&(E_grid[gid][2]==1)) {
                nodes[cnt][0] = x; nodes[cnt][1] = y; cnt++;
            }
            if((E_grid[gid][1]==1)&&(E_grid[gid][2]==1)) {
                nodes[cnt][0] = x+1; nodes[cnt][1] = y; cnt++;
            }
            if((E_grid[gid][0]==1)&&(E_grid[gid][3]==1)) {
                nodes[cnt][0] = x; nodes[cnt][1] = y+1; cnt++;
            }
            if((E_grid[gid][1]==1)&&(E_grid[gid][3]==1)) {
                nodes[cnt][0] = x+1; nodes[cnt][1] = y+1; cnt++;
            }

            tmp_S0 += 1;
        }
        cnt = getConvexHull(cnt, nodes);
        tmp_S1 = getSofPoly(cnt, nodes);
        S1 += tmp_S1;
        S0 += tmp_S0;
        if(save){
            label_pairs[li][0] = tmp_S0;
            label_pairs[li][1] = tmp_S1;
        }
    }
    std::vector<double> S_pair(2, 0);
    S_pair[0] = (S1-S0);
    S_pair[1] = S1;

    delete[] head;
    delete[] last;
    delete[] element_asses;
    delete[] E_grid;
    delete[] nodes;
    return S_pair;
}

//calculate the full cost, using convexity based on area(where S means Square) of convex hull
//save: if save convexity measure pairs of every cluster
//load: if load saved convexity measure pairs of every cluster and only re-calculate clusters of "mainLabel1" and "mainLabel2"
//mainLabel1, mainLabel2: re-calculate clusters of "mainLabel1" and "mainLabel2" when load saved convexity measure pairs
//return cost[0..3], full/Similar/Compact/Convex cost
std::vector<double> checkCostForS(
    const double Similar_cost_matrix[],
    const double Compact_cost_matrix[],
    const int grid_asses[], const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    const double &alpha, const double &beta,
    bool save=false, bool load=false, double label_pairs[][2]=nullptr, int mainLabel1=-1, int mainLabel2=-1) {

    std::vector<double> S_pair = checkConvexForSArray(grid_asses, cluster_labels, N, num, square_len, maxLabel, save, load, label_pairs, mainLabel1, mainLabel2);
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

    delete[] element_asses;
    std::vector<double> ret(4, 0);
    ret[0] = cost;
    ret[1] = Similar_cost;
    ret[2] = Compact_cost;
    ret[3] = Convex_cost*N;
    return ret;
}
#endif