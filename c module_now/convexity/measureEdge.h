#ifndef _MEASURE_EDGE_H
#define _MEASURE_EDGE_H

#include "../utils/base.h"

//grid_asses[0..N-1]: elements in every grid
//cluster_labels[0..num-1]: cluster label of every element
//N: number of grids
//num: number of elements
//square_len: sqrt(N)
//maxLabel: number of clusters


//calculate the convexity of edge-tangent-divide in grid "gid" and four adjacent grids
//x_pre,y_pre[0..maxLabel*(square_len+1)-1]: saved prefrocess information before, row/column prefix of count of elements in every clusters
//return E[0..1], convexity measure pairs
std::vector<double> checkConvexForESingle4(
const int grid_asses[],
const int cluster_labels[],
const int x_pre[], const int y_pre[], const int &gid, const int &lb,
const int &N, const int &num, const int &square_len, const int &maxLabel) {
    double E1 = 0;
    double E0 = 0;
    int i = gid/square_len;
    int j = gid%square_len;

    int bias = i*square_len;
    int bias0 = i*maxLabel;
    int bias1 = (i+1)*maxLabel;
    int bias2 = square_len*maxLabel;

//    if(lb<maxLabel) {
    if(true) {

        if(i-1>=0){
            int gid1 = gid-square_len;
            int id1 = grid_asses[gid1];
            if((id1>=num)||(cluster_labels[id1]!=lb)){
                double t = 1.0*x_pre[bias0+lb]/x_pre[bias2+lb];

                if((0<=lb)&&(lb<maxLabel)) {
                    E1 += 1;
                    E0 += t;
                }

                if(id1<num) {
                    int lb2 = cluster_labels[id1];
                    E1 += 1;
                    t = 1.0*(x_pre[bias2+lb2]-x_pre[bias0+lb2])/x_pre[bias2+lb2];
                    E0 += t;
                }
            }
        }else E1 += 1;

        if(i+1<square_len){
            int gid1 = gid+square_len;
            int id1 = grid_asses[gid1];
            if((id1>=num)||(cluster_labels[id1]!=lb)){
                double t = 1.0*(x_pre[bias2+lb]-x_pre[bias1+lb])/x_pre[bias2+lb];

                if((0<=lb)&&(lb<maxLabel)) {
                    E1 += 1;
                    E0 += t;
                }

                if(id1<num) {
                    int lb2 = cluster_labels[id1];
                    E1 += 1;
                    t = 1.0*x_pre[bias1+lb2]/x_pre[bias2+lb2];
                    E0 += t;
                }
            }
        }else E1 += 1;

        if(j-1>=0){
            int gid1 = gid-1;
            int id1 = grid_asses[gid1];
            if((id1>=num)||(cluster_labels[id1]!=lb)){
                double t = 1.0*y_pre[j*maxLabel+lb]/y_pre[bias2+lb];

                if((0<=lb)&&(lb<maxLabel)) {
                    E1 += 1;
                    E0 += t;
                }

                if(id1<num) {
                    int lb2 = cluster_labels[id1];
                    E1 += 1;
                    t = 1.0*(y_pre[bias2+lb2]-y_pre[j*maxLabel+lb2])/y_pre[bias2+lb2];
                    E0 += t;
                }
            }
        }else E1 += 1;

        if(j+1<square_len){
            int gid1 = gid+1;
            int id1 = grid_asses[gid1];
            if((id1>=num)||(cluster_labels[id1]!=lb)){
                double t = 1.0*(y_pre[bias2+lb]-y_pre[(j+1)*maxLabel+lb])/y_pre[bias2+lb];

                if((0<=lb)&&(lb<maxLabel)) {
                    E1 += 1;
                    E0 += t;
                }

                if(id1<num) {
                    int lb2 = cluster_labels[id1];
                    E1 += 1;
                    t = 1.0*y_pre[(j+1)*maxLabel+lb2]/y_pre[bias2+lb2];
                    E0 += t;
                }
            }
        }else E1 += 1;

    }

    std::vector<double> E_pair(2, 0);
    E_pair[0] = E0;
    E_pair[1] = E1;
    return E_pair;
}

//calculate the convexity of edge-tangent-divide
//save: if save prefix information
//save_x_pre, save_y_pre[0..maxLabel*(square_len+1)-1]: prefix
//return E[0..1]: convexity(in the form of loss/cost) pair
std::vector<double> checkConvexForEArray(
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
bool save=false, int save_x_pre[]=nullptr, int save_y_pre[]=nullptr) {

    int *x_pre = new int[(square_len+1)*maxLabel];
    int *y_pre = new int[(square_len+1)*maxLabel];
    for(int i=0;i<(square_len+1)*maxLabel;i++){
        x_pre[i] = y_pre[i] = 0;
    }

    // #pragma omp parallel for num_threads(THREADS_NUM)
    for(int i=0;i<square_len;i++){
        int bias1 = (i+1)*maxLabel;
        for(int j=0;j<square_len;j++){
            int gid1 = i*square_len+j;
            int id1 = grid_asses[gid1];
            int gid2 = j*square_len+i;
            int id2 = grid_asses[gid2];
            if(id1<num)x_pre[bias1+cluster_labels[id1]] += 1;
            if(id2<num)y_pre[bias1+cluster_labels[id2]] += 1;
        }
    }

    // #pragma omp parallel for num_threads(THREADS_NUM)
    for(int j=0;j<maxLabel;j++){
        for(int i=0;i<square_len;i++){
            int bias0 = i*maxLabel;
            int bias1 = (i+1)*maxLabel;
            x_pre[bias1+j] += x_pre[bias0+j];
            y_pre[bias1+j] += y_pre[bias0+j];
        }
    }

    double E1 = 0;
    double E0 = 0;

    // #pragma omp parallel for reduction(+:E0, E1) num_threads(THREADS_NUM)
    for(int i=0;i<square_len;i++){
        int bias = i*square_len;
        int bias0 = i*maxLabel;
        int bias1 = (i+1)*maxLabel;
        int bias2 = square_len*maxLabel;
        for(int j=0;j<square_len;j++){
            int gid = bias+j;
            int id = grid_asses[gid];
            if(id>=num)continue;
            int lb = cluster_labels[id];

            if(i-1>=0){
                int gid1 = gid-square_len;
                int id1 = grid_asses[gid1];
                if((id1>=num)||(cluster_labels[id1]!=lb)){
                    E1 += 1;
                    double t = 1.0*x_pre[bias0+lb]/x_pre[bias2+lb];
                    E0 += t;
                }
            }else E1 += 1;

            if(i+1<square_len){
                int gid1 = gid+square_len;
                int id1 = grid_asses[gid1];
                if((id1>=num)||(cluster_labels[id1]!=lb)){
                    E1 += 1;
                    double t = 1.0*(x_pre[bias2+lb]-x_pre[bias1+lb])/x_pre[bias2+lb];
                    E0 += t;
                }
            }else E1 += 1;

            if(j-1>=0){
                int gid1 = gid-1;
                int id1 = grid_asses[gid1];
                if((id1>=num)||(cluster_labels[id1]!=lb)){
                    E1 += 1;
                    double t = 1.0*y_pre[j*maxLabel+lb]/y_pre[bias2+lb];
                    E0 += t;
                }
            }else E1 += 1;

            if(j+1<square_len){
                int gid1 = gid+1;
                int id1 = grid_asses[gid1];
                if((id1>=num)||(cluster_labels[id1]!=lb)){
                    E1 += 1;
                    double t = 1.0*(y_pre[bias2+lb]-y_pre[(j+1)*maxLabel+lb])/y_pre[bias2+lb];
                    E0 += t;
                }
            }else E1 += 1;
        }
    }

    std::vector<double> E_pair(2, 0);
    E_pair[0] = E0;
    E_pair[1] = E1;

    if(save){
        for(int i=0;i<(square_len+1)*maxLabel;i++){
            save_x_pre[i] = x_pre[i];
            save_y_pre[i] = y_pre[i];
        }
    }
    delete[] x_pre;
    delete[] y_pre;
    return E_pair;
}

//calculate the convexity of edge-tangent-divide, python interface
std::vector<double> checkConvexForE(
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels) {
    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = max(maxLabel, _cluster_labels[i]+1);

    int *grid_asses = new int[N];
    int *cluster_labels = new int[num];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];

    std::vector<double> ret = checkConvexForEArray(grid_asses, cluster_labels, N, num, square_len, maxLabel);

    delete[] grid_asses;
    delete[] cluster_labels;
    return ret;
}

//calculate the convexity cost matrix of edge-tangent-divide
//cost_matrix_a[0..N*N-1]: cost matrix will be stored in it
void getCostMatrixForEArrayToArray(
int grid_asses[],
int cluster_labels[],
double cost_matrix_a[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    int *x_pre = new int[(square_len+1)*maxLabel];
    int *y_pre = new int[(square_len+1)*maxLabel];
    std::vector<double> E_pair1 = checkConvexForEArray(grid_asses, cluster_labels, N, num, square_len, maxLabel, true, x_pre, y_pre);
    double correct1=0, full1=0;
    full1 = E_pair1[1];
    correct1 = E_pair1[1]-E_pair1[0];

    double *correct_matrix1 = new double[N];    //当前layout上每个位置上的边分割
    double *full_matrix1 = new double[N];

    // #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid=0;gid<N;gid++){
        int lb = maxLabel;
        if(grid_asses[gid]<num)lb = cluster_labels[grid_asses[gid]];
        std::vector<double> pair = checkConvexForESingle4(grid_asses, cluster_labels, x_pre, y_pre, gid, lb, N, num, square_len, maxLabel);
        full_matrix1[gid] = pair[1];
        correct_matrix1[gid] = pair[1]-pair[0];
    }

    double *correct_matrix2 = new double[N*(maxLabel+1)];    //对于同label的元素放置在同个位置上，只需计算一次，因此进行记忆
    double *full_matrix2 = new double[N*(maxLabel+1)];
    for(int i=0;i<N*(maxLabel+1);i++){
        correct_matrix2[i] = full_matrix2[i] = -1;
    }

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid2=0;gid2<N;gid2++){
        int bias = gid2*N;
        int bias1 = gid2*(maxLabel+1);
        for(int gid1=0;gid1<N;gid1++){
            int x1 = gid1/square_len;
            int y1 = gid1%square_len;
            int x2 = gid2/square_len;
            int y2 = gid2%square_len;
            if((abs(x1-x2)>10)||(abs(y1-y2)>10)){
                cost_matrix_a[bias+grid_asses[gid1]] = 1;
                continue;
            }
            int lb = maxLabel;
            int id1 = grid_asses[gid1];
            if(id1<num)lb = cluster_labels[id1];
            if(full_matrix2[bias1+lb]<0){
                // int tmp = grid_asses[gid2];
                // grid_asses[gid2] = id1;
                std::vector<double> E_pair2 = checkConvexForESingle4(grid_asses, cluster_labels, x_pre, y_pre, gid2, lb, N, num, square_len, maxLabel);
                double correct2=0, full2=0;
                full2 = E_pair2[1];
                correct2 = E_pair2[1]-E_pair2[0];
                correct_matrix2[bias1+lb] = correct2;
                full_matrix2[bias1+lb] = full2;
                // grid_asses[gid2] = tmp;
            }
            double correct2 = correct1-correct_matrix1[gid2]+correct_matrix2[bias1+lb];
            double full2 = full1-full_matrix1[gid2]+full_matrix2[bias1+lb];
            cost_matrix_a[bias+id1] = (1.0*(correct1-correct2)/correct2+1.0*(full2-full1)/full1)*N;
        }
    }

    double mini = 0;

    #pragma omp parallel for reduction(min:mini) num_threads(THREADS_NUM)
    for(int gid2=0;gid2<N;gid2++){
        int bias = gid2*N;
        for(int id1=0;id1<N;id1++) mini = min(mini, cost_matrix_a[bias+id1]);
    }

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid2=0;gid2<N;gid2++){
        int bias = gid2*N;
        for(int id1=0;id1<N;id1++) cost_matrix_a[bias+id1] -= mini;
    }

    delete[] x_pre;
    delete[] y_pre;
    delete[] correct_matrix1;
    delete[] full_matrix1;
    delete[] correct_matrix2;
    delete[] full_matrix2;
    return;
}

//calculate the convexity cost matrix of edge-tangent-divide, python interface
std::vector<std::vector<double>> getCostMatrixForE(
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels) {
    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = max(maxLabel, _cluster_labels[i]+1);
    int *grid_asses = new int[N];
    int *cluster_labels = new int[num];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];

    std::vector<std::vector<double>> cost_matrix(N, std::vector<double>(N, 0));
    double *cost_matrix_a = new double [N*N];
    getCostMatrixForEArrayToArray(grid_asses, cluster_labels, cost_matrix_a, N, num, square_len, maxLabel);

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid=0;gid<N;gid++) {
        int bias = gid*N;
        for(int id=0;id<N;id++){
            cost_matrix[gid][id] = cost_matrix_a[bias+id];
        }
    }

    delete[] grid_asses;
    delete[] cluster_labels;
    delete[] cost_matrix_a;
    return cost_matrix;
}

//calculate full cost, using convexity of edge-tangent-divide
//alpha, beta: cost arguments
//return cost[0..3], full/Similar/Compact/Convex cost
std::vector<double> checkCostForE(
    const double Similar_cost_matrix[],
    const double Compact_cost_matrix[],
    const int grid_asses[], const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    const double &alpha, const double &beta) {

    std::vector<double> E_pair = checkConvexForEArray(grid_asses, cluster_labels, N, num, square_len, maxLabel);
    double correct=0, full=0;
    full = E_pair[1];
    correct = E_pair[1]-E_pair[0];

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