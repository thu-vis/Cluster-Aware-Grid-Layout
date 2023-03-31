#ifndef _MEASURE_DOUBLES_H
#define _MEASURE_DOUBLES_H

#include "measureTriples.h"

//std::vector<int> getDoubles(double a, double b) {
//    double tiny = 0.000001;
//    std::vector<int> T_pair(2, 0);
//    if(b>tiny) {
//        T_pair[1] = 1;
//        if(a>tiny)T_pair[0] = 1;
//    }
//    return T_pair;
//}

void getDoubles(double &T0, double &T1, double a, double b) {
    double tiny = 0.000001;
    T0 = 0, T1 = 0;
    if(b>tiny) {
        T1 = 1;
        if(a>tiny)T0 = 1;
    }
}

int getDoubleId(int x, int y, int N){
    if(x>y){ int tmp=x; x=y; y=tmp;}
    return (x*N+y)*2;
}

//check convexity of triples
//return T[0..1]: convexity(in the form of loss/cost) pair in layout now
std::vector<double> checkConvexForT2Array(
const int triples[][4],
const int &cnt,
bool save,
double save_D_pair[],  // N*N*2
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    double *D_pair = new double[N*N*2];
    for(int i=0;i<N*N*2;i++)D_pair[i] = 0;

    double T0 = 0;
    double T1 = 0;

    for(int starti=0;starti<N;starti++){
        for(int i=starti;i<cnt;i+=N){
            int gid1 = triples[i][0];
            int gid0 = triples[i][1];
            int gid2 = triples[i][2];
            int times = triples[i][3];

            // printf("triple %d %d %d %d\n", gid1, gid0, gid2, times);

            int did = getDoubleId(gid1, gid2, N);

            int id1 = grid_asses[gid1];
            int id0 = grid_asses[gid0];
            int id2 = grid_asses[gid2];
            if((id1>=num)||(id2>=num))continue;
            int lb1 = cluster_labels[id1];
            int lb2 = cluster_labels[id2];
            int lb0 = -1;
            if(id0<num)lb0 = cluster_labels[id0];
            if(lb1==lb2){
                D_pair[did+1] += times;
                if(lb1!=lb0)D_pair[did] += times;
            }
        }

    }

    for(int i=0;i<N;i++)
    for(int j=i+1;j<N;j++) {
        int did = getDoubleId(i, j, N);
        //std::vector<int> D_tmp = getDoubles(D_pair[did], D_pair[did+1]);
        double D_tmp[2];
        getDoubles(D_tmp[0], D_tmp[1], D_pair[did], D_pair[did+1]);
        // printf("%d %d %.2lf %.2lf %d %d\n", i, j, D_pair[did], D_pair[did+1], D_tmp[0], D_tmp[1]);
        T0 += D_tmp[0];
        T1 += D_tmp[1];
    }
    std::vector<double> T_pair(2, 0);
    T_pair[0] = T0;
    T_pair[1] = T1;

    if(save) {
        for(int i=0;i<N*N*2;i++)save_D_pair[i] = D_pair[i];
    }

    delete[] D_pair;
    return T_pair;
}

//check convexity of triples with saved information of old layout
//old_T_pair[0..1]: convexity(in the form of loss/cost) pair in old grid layout
//return T[0..1]: convexity(in the form of loss/cost) pair in layout now
std::vector<double> checkConvexForT2withOld(
int triples[][4],
int triples_head[],
int triples_list[][2],
bool save,
double old_T_pair[],
double old_D_pair[],  // N*N*2
const int grid_asses[],
const int old_grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    int *changed = new int[N];
    for(int gid=0;gid<N;gid++){
        int old_lb = -1;
        if(old_grid_asses[gid]<num)old_lb = cluster_labels[old_grid_asses[gid]];
        int lb = -1;
        if(grid_asses[gid]<num)lb = cluster_labels[grid_asses[gid]];
        if(old_lb==lb)changed[gid] = 0;
        else changed[gid] = 1;
    }

    double T0 = 0, T1 = 0;
    double dec_T0 = 0, dec_T1 = 0;

    for(int gid=0;gid<N;gid++)
    if(changed[gid]==1) {
        for(int ii=triples_head[gid];ii>=0;ii=triples_list[ii][1]){
            int i = triples_list[ii][0];
            int gid1 = triples[i][0];
            int gid0 = triples[i][1];
            int gid2 = triples[i][2];
            double times = triples[i][3];
            times = (1.0*times)/(changed[gid1]+changed[gid0]+changed[gid2]);

            int did = getDoubleId(gid1, gid2, N);
            double tmp_T0 = old_D_pair[did];
            double tmp_T1 = old_D_pair[did+1];

            int id1 = grid_asses[gid1];
            int id0 = grid_asses[gid0];
            int id2 = grid_asses[gid2];
            if((id1<num)&&(id2<num)) {
                int lb1 = cluster_labels[id1];
                int lb2 = cluster_labels[id2];
                int lb0 = -1;
                if(id0<num)lb0 = cluster_labels[id0];
                if(lb1==lb2){
                    tmp_T1 += times;
                    if(lb1!=lb0)tmp_T0 += times;
                }
            }

            id1 = old_grid_asses[gid1];
            id0 = old_grid_asses[gid0];
            id2 = old_grid_asses[gid2];
            if((id1<num)&&(id2<num)) {
                int lb1 = cluster_labels[id1];
                int lb2 = cluster_labels[id2];
                int lb0 = -1;
                if(id0<num)lb0 = cluster_labels[id0];
                if(lb1==lb2){
                    tmp_T1 -= times;
                    if(lb1!=lb0)tmp_T0 -= times;
                }
            }
//            std::vector<int> D_tmp1 = getDoubles(old_D_pair[did], old_D_pair[did+1]);
//            std::vector<int> D_tmp2 = getDoubles(tmp_T0, tmp_T1);
            double D_tmp1[2], D_tmp2[2];
            getDoubles(D_tmp1[0], D_tmp1[1], old_D_pair[did], old_D_pair[did+1]);
            getDoubles(D_tmp2[0], D_tmp2[1], tmp_T0, tmp_T1);
            T0 += D_tmp2[0];
            T1 += D_tmp2[1];
            dec_T0 += D_tmp1[0];
            dec_T1 += D_tmp1[1];
            old_D_pair[did] = tmp_T0;
            old_D_pair[did+1] = tmp_T1;
        }
    }

    std::vector<double> T_pair(2, 0);
    T_pair[0] = old_T_pair[0]+T0-dec_T0;
    T_pair[1] = old_T_pair[1]+T1-dec_T1;

    if(save) {
        delete[] changed;
        return T_pair;
    }

    for(int gid=0;gid<N;gid++)
    if(changed[gid]==1) {
        for(int ii=triples_head[gid];ii>=0;ii=triples_list[ii][1]){
            int i = triples_list[ii][0];
            int gid1 = triples[i][0];
            int gid0 = triples[i][1];
            int gid2 = triples[i][2];
            double times = triples[i][3];
            times = (1.0*times)/(changed[gid1]+changed[gid0]+changed[gid2]);

            int did = getDoubleId(gid1, gid2, N);
            double tmp_T0 = old_D_pair[did];
            double tmp_T1 = old_D_pair[did+1];

            int id1 = grid_asses[gid1];
            int id0 = grid_asses[gid0];
            int id2 = grid_asses[gid2];
            if((id1<num)&&(id2<num)) {
                int lb1 = cluster_labels[id1];
                int lb2 = cluster_labels[id2];
                int lb0 = -1;
                if(id0<num)lb0 = cluster_labels[id0];
                if(lb1==lb2){
                    tmp_T1 -= times;
                    if(lb1!=lb0)tmp_T0 -= times;
                }
            }

            id1 = old_grid_asses[gid1];
            id0 = old_grid_asses[gid0];
            id2 = old_grid_asses[gid2];
            if((id1<num)&&(id2<num)) {
                int lb1 = cluster_labels[id1];
                int lb2 = cluster_labels[id2];
                int lb0 = -1;
                if(id0<num)lb0 = cluster_labels[id0];
                if(lb1==lb2){
                    tmp_T1 += times;
                    if(lb1!=lb0)tmp_T0 += times;
                }
            }
            old_D_pair[did] = tmp_T0;
            old_D_pair[did+1] = tmp_T1;
        }
    }

    delete[] changed;
    return T_pair;
}

//calculate full cost, using convexity of triples
//alpha, beta: cost arguments
//return cost[0..3], full/Similar/Compact/Convex cost
std::vector<double> checkCostForT2(
const double Similar_cost_matrix[],
const double Compact_cost_matrix[],
const int grid_asses[], const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
const double &alpha, const double &beta,
bool save=false, bool load=false, int old_grid_asses[]=nullptr, double old_T_pair[]=nullptr, double old_D_pair[]=nullptr) {

    int (*triples)[4];
    int *triples_head;
    int (*triples_list)[2];
    int cnt;

    if(global_N==N){
        triples = global_triples;
        triples_head = global_triples_head;
        triples_list = global_triples_list;
        cnt = global_cnt;
    }else{
        if(global_cnt>=0){
            delete[] global_triples;
            delete[] global_triples_head;
            delete[] global_triples_list;
            global_cnt = -1;
            global_N = -1;
        }
        triples = new int[N*N*square_len][4];
        triples_head = new int[N];
        triples_list = new int[N*N*square_len*3][2];
        cnt = getTriples(triples, triples_head, triples_list, grid_asses, cluster_labels, N, num, square_len, maxLabel);
        global_triples = triples;
        global_triples_head = triples_head;
        global_triples_list = triples_list;
        global_cnt = cnt;
        global_N = N;
    }
    // printf("triples count: %d\n", cnt);

    std::vector<double> T_pair(2, 0);
    if(!load){
        T_pair = checkConvexForT2Array(triples, cnt, save, old_D_pair, grid_asses, cluster_labels, N, num, square_len, maxLabel);
        // printf("T_pair %.2lf %.2lf\n", T_pair[0], T_pair[1]);
    }else {
        T_pair = checkConvexForT2withOld(triples, triples_head, triples_list,
        save, old_T_pair, old_D_pair, grid_asses, old_grid_asses, cluster_labels, N, num, square_len, maxLabel);
    }

    if(save){
        old_T_pair[0] = T_pair[0];
        old_T_pair[1] = T_pair[1];
    }

    // printf("D_pair %d\n", save);
    // for(int i=0;i<N*N*2;i++)printf("%.2lf ", old_D_pair[i]);
    // printf("\n");
    // printf("T_pair %d\n", save);
    // for(int i=0;i<2;i++)printf("%.2lf ", old_T_pair[i]);
    // printf("\n");

    double correct=0, full=0;

    full = T_pair[1];
    correct = T_pair[1]-T_pair[0];

    double Convex_cost;
    if(full<0.000001)Convex_cost = 0;
    else Convex_cost = (full-correct)/full;

    // printf("Convex_cost %.2lf\n", Convex_cost);
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
        // printf("cost %.2lf\n", cost);
    }
    // printf("Convex_cost %.2lf N %d alpha %.2lf\n", Convex_cost, N, alpha);
    cost += Convex_cost * N * alpha;
    // printf("cost %.2lf\n", cost);
    // printf("cost %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n", Similar_cost, Compact_cost, Convex_cost*N, beta, alpha, cost);

    delete[] element_asses;

    std::vector<double> ret(4, 0);
    ret[0] = cost;
    ret[1] = Similar_cost;
    ret[2] = Compact_cost;
    ret[3] = Convex_cost*N;
    return ret;
}

# endif