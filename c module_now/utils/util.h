#ifndef _MEASURE_OTHER_H
#define _MEASURE_OTHER_H

#include "base.h"

//DFS check connectivity
void checkConnectDFS(
    const int grid_asses[],
    const int cluster_labels[],
    int checked[], const int &check_cnt,
    const int &x, const int &y, const int &startgid,
    const int &lb, const int &num, const int &square_len, const int &connect) {
    int gid = x*square_len+y;
    int id = grid_asses[gid];
    if((gid!=startgid)&&((id>=num)||(cluster_labels[id]!=lb)))return;
    if(checked[gid]>=check_cnt)return;
    checked[gid] = check_cnt;

    int f_x[8] = {-1, 0, 0, 1, 1, -1, 1, -1};
    int f_y[8] = {0, -1, 1, 0, 1, 1, -1, -1};
    for(int i=0;i<connect;i++){
        int x2 = x+f_x[i];
        int y2 = y+f_y[i];
        if((x2<0)||(x2>=square_len)||(y2<0)||(y2>=square_len))continue;
        int gid2 = x2*square_len+y2;
        checkConnectDFS(grid_asses, cluster_labels, checked, check_cnt, x2, y2, startgid, lb, num, square_len, connect);
    }
}

//check connectivity of "gid" grid with the cluster of label "lb"
double checkConnect(
const int grid_asses[],
const int cluster_labels[],
int checked[],
const int &gid, const int &lb,
const int &N, const int &num, const int &square_len, const int &maxLabel,
int check_cnt=1) {
    if(check_cnt==1)
        for(int i=0;i<N;i++)checked[i] = 0;
    int x = gid/square_len;
    int y = gid%square_len;
    int connect = 4;
    checkConnectDFS(grid_asses, cluster_labels, checked, check_cnt, x, y, gid, lb, num, square_len, connect);

    int labelCheckCnt = 0;
    int tot = 0;
    double cluster_x=0, cluster_y=0;
    for(int gid1=0;gid1<N;gid1++) {
        int id1 = grid_asses[gid1];
        if((id1<num)&&(cluster_labels[id1]==lb)){
            tot += 1;
            double x = 1.0*(gid1/square_len)/square_len;
            double y = 1.0*(gid1%square_len)/square_len;
            cluster_x += x;
            cluster_y += y;
        }
    }
    cluster_x /= tot;
    cluster_y /= tot;

    for(int gid1=0;gid1<N;gid1++){
        if(checked[gid1]>=check_cnt)labelCheckCnt += 1;
    }
    if(tot==0)return 0;
    double per = 1.0*labelCheckCnt/tot;

    if(per>0.5)return 0;
    else {
        double x = 1.0*(gid/square_len)/square_len;
        double y = 1.0*(gid%square_len)/square_len;
        double dist = getDist(x, y, cluster_x, cluster_y);
        return dist*dist;
    }
}

//check connectivity of every grid with the cluster of its label
double checkConnectForAll(
    const int grid_asses[],
    const int cluster_labels[],
    int checked[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    int connect=8,
    bool if_disconn[]=nullptr) {

    int *a = new int[N];
    int *b = new int[maxLabel];
    for(int i=0;i<N;i++)a[i] = 0;
    for(int i=0;i<maxLabel;i++)b[i] = 0;

    int check_cnt = 1;
    for(int i=0;i<N;i++)checked[i] = 0;
    if(if_disconn!=nullptr){
        for(int i=0;i<N;i++)if_disconn[i] = false;
    }

    for(int gid=0;gid<N;gid++)
    if(grid_asses[gid]<num){
        int id = grid_asses[gid];
        int lb = cluster_labels[id];
        if(checked[gid]>0)continue;
        int x = gid/square_len;
        int y = gid%square_len;
        checkConnectDFS(grid_asses, cluster_labels, checked, check_cnt++, x, y, gid, lb, num, square_len, connect);
    }

    for(int gid=0;gid<N;gid++)
    if(grid_asses[gid]<num){
        int id = grid_asses[gid];
        int lb = cluster_labels[id];
        b[lb] += 1;
        a[checked[gid]-1] += 1;
    }

    double *cluster_x = new double[maxLabel];
    double *cluster_y = new double[maxLabel];
    int *cluster_cnt = new int[maxLabel];
    for(int i=0;i<maxLabel;i++){
        cluster_x[i] = 0;
        cluster_y[i] = 0;
        cluster_cnt[i] = 0;
    }

    for(int gid=0;gid<N;gid++){    // every positon
        double x = 1.0*(gid/square_len)/square_len;
        double y = 1.0*(gid%square_len)/square_len;
        int id = grid_asses[gid];
        if(id >= num)continue;
        int label = cluster_labels[id];
        cluster_x[label] += x;
        cluster_y[label] += y;
        cluster_cnt[label] += 1;
    }
    for(int i=0;i<maxLabel;i++){
        if(cluster_cnt[i]==0)continue;
        cluster_x[i] /= cluster_cnt[i];
        cluster_y[i] /= cluster_cnt[i];
    }

    double ret = 0;
    for(int gid=0;gid<N;gid++)
    if(grid_asses[gid]<num){
        int id = grid_asses[gid];
        int lb = cluster_labels[id];
        double per = 1.0*a[checked[gid]-1]/b[lb];
        if(per<=0.5){
//            ret += 1;
            double x = 1.0*(gid/square_len)/square_len;
            double y = 1.0*(gid%square_len)/square_len;
            double dist = getDist(x, y, cluster_x[lb], cluster_y[lb]);
            ret += dist*dist;
            if(if_disconn!=nullptr){
                if_disconn[gid] = true;
            }
            // printf("cn wrong %d\n", gid);
        }
        if(b[lb]<=1){
            if(if_disconn!=nullptr){
                if_disconn[gid] = true;
            }
        }
    }

    delete[] cluster_x;
    delete[] cluster_y;
    delete[] cluster_cnt;
    delete[] a;
    delete[] b;
    return ret;
}

//get connectivity cost matrix（constraint）
void getConnectCostMatrixArrayToArray(
    int grid_asses[],
    int cluster_labels[],
    double ret_a[],
    const int &N, const int &num, const int &square_len, const int &maxLabel) {

    double *cost_matrix = new double[N*(maxLabel+1)];  //elements of same label matching a same grid, have same cost
    for(int i=0;i<N*(maxLabel+1);i++)cost_matrix[i] = -1;

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid2=0;gid2<N;gid2++){    // grids
        int x2 = gid2/square_len;
        int y2 = gid2%square_len;

        int *checked = new int[N];

        for(int gid1=0;gid1<N;gid1++) {    // position of elements in layout now
            int id1 = grid_asses[gid1];
            int x1 = gid1/square_len;
            int y1 = gid1%square_len;

            int lb = maxLabel;
            if(id1<num)lb = cluster_labels[id1];    // label(cluster) of the element
            if(cost_matrix[gid2*(maxLabel+1)+lb]<0){
                // int tmp = grid_asses[gid2];
                // grid_asses[gid2] = id1;
                double cost = 0;
                if(id1<num){
                    cost = checkConnect(grid_asses, cluster_labels, checked, gid2, lb, N, num, square_len, maxLabel);
                }
                cost_matrix[gid2*(maxLabel+1)+lb] = cost;
            }
            ret_a[gid2*N+id1] = cost_matrix[gid2*(maxLabel+1)+lb];
        }

        delete[] checked;
    }

    delete[] cost_matrix;
    return;
}

//connectivity cost matrix, python interface
std::vector<std::vector<double>> getConnectCostMatrix(
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

    std::vector<std::vector<double>> ret(N, std::vector<double>(N, 0));
    double *ret_a = new double [N*N];
    getConnectCostMatrixArrayToArray(grid_asses, cluster_labels, ret_a, N, num, square_len, maxLabel);

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid=0;gid<N;gid++){
        int bias = gid*N;
        for(int id=0;id<N;id++){
            ret[gid][id] = ret_a[bias+id];
        }
    }

    delete[] grid_asses;
    delete[] cluster_labels;
    delete[] ret_a;
    return ret;
}

//origin cost matrix(distance to the origin layout)
void getOriginCostMatrixArrayToArray(
//int grid_asses[],
double ori_embedded[][2],
int cluster_labels[],
double ret_a[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid1=0;gid1<N;gid1++) {
        double x1 = 1.0*(gid1/square_len)/square_len;
        double y1 = 1.0*(gid1%square_len)/square_len;
        int bias = gid1*N;
//        for(int gid2=0;gid2<N;gid2++){
//            int id = grid_asses[gid2];
//            double x2 = 1.0*(gid2/square_len)/square_len;
//            double y2 = 1.0*(gid2%square_len)/square_len;
        for(int id=0;id<N;id++){
            double x2 = ori_embedded[id][0];
            double y2 = ori_embedded[id][1];
            int lb = cluster_labels[id];
            ret_a[bias+id] = getDist(x1, y1, x2, y2);
            ret_a[bias+id] = ret_a[bias+id]*ret_a[bias+id];
        }
    }

    return;
}

//compact cost matrix(distance from element to center of cluster)
void getCompactCostMatrixArrayToArray(
int grid_asses[],
int cluster_labels[],
double ret_a[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    double *cluster_x = new double[maxLabel];
    double *cluster_y = new double[maxLabel];
    int *cluster_cnt = new int[maxLabel];
    for(int i=0;i<maxLabel;i++){
        cluster_x[i] = 0;
        cluster_y[i] = 0;
        cluster_cnt[i] = 0;
    }

    for(int gid=0;gid<N;gid++){    // every positon
        double x = 1.0*(gid/square_len)/square_len;
        double y = 1.0*(gid%square_len)/square_len;
        int id = grid_asses[gid];
        if(id >= num)continue;
        int label = cluster_labels[id];
        cluster_x[label] += x;
        cluster_y[label] += y;
        cluster_cnt[label] += 1;
    }
    for(int i=0;i<maxLabel;i++){
        if(cluster_cnt[i]==0)continue;
        cluster_x[i] /= cluster_cnt[i];
        cluster_y[i] /= cluster_cnt[i];
    }

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid=0;gid<N;gid++) {    // grids
        int bias = gid*N;
        for(int id=0;id<num;id++){    //elements
            double x = 1.0*(gid/square_len)/square_len;
            double y = 1.0*(gid%square_len)/square_len;
            int lb = cluster_labels[id];
            ret_a[bias+id] = getDist(x, y, cluster_x[lb], cluster_y[lb]);
            ret_a[bias+id] = ret_a[bias+id]*ret_a[bias+id];
        }
        for(int id=num;id<N;id++)ret_a[bias+id] = 0;
    }

    delete[] cluster_x;
    delete[] cluster_y;
    delete[] cluster_cnt;
    return;
}

//compact cost matrix, python interface
std::vector<std::vector<double>> getCompactCostMatrix(
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

    std::vector<std::vector<double>> ret(N, std::vector<double>(N, 0));
    double *ret_a = new double [N*N];
    getCompactCostMatrixArrayToArray(grid_asses, cluster_labels, ret_a, N, num, square_len, maxLabel);

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid=0;gid<N;gid++) {
        int bias = gid*N;
        for(int id=0;id<num;id++){
            ret[gid][id] = ret_a[bias+id];
        }
    }

    delete[] grid_asses;
    delete[] cluster_labels;
    delete[] ret_a;
    return ret;
}

//count the blank-nonblank edges of grid "gid"
double checkBlankForGrid(
const int &gid,
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    int i = gid/square_len;
    int j = gid%square_len;

    int ans=0;

    int flag=0;
    if(grid_asses[gid]<num)flag=1;

    if(i-1>=0){    // up
        int gid1 = gid-square_len;
        int id1 = grid_asses[gid1];
        if(id1<num){
            ans += flag^1;
        }else ans += flag;
    }else ans += flag;

    if(i+1<square_len){    //down
        int gid1 = gid+square_len;
        int id1 = grid_asses[gid1];
        if(id1<num){
            ans += flag^1;
        }else ans += flag;
    }else ans += flag;

    if(j-1>=0){    //left
        int gid1 = gid-1;
        int id1 = grid_asses[gid1];
        if(id1<num){
            ans += flag^1;
        }else ans += flag;
    }else ans += flag;

    if(j+1<square_len){    //right
        int gid1 = gid+1;
        int id1 = grid_asses[gid1];
        if(id1<num){
            ans += flag^1;
        }else ans += flag;
    }else ans += flag;

    return ans;
}

//count the blank-nonblank edges of all grids
double checkBlankForAll(
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    int ans=0;

    for(int i=0;i<square_len;i++) {
        int bias = i*square_len;
        for(int j=0;j<square_len;j++) {
            int gid = bias+j;
            int flag=0;
            if(grid_asses[gid]<num)flag=1;

            if(i-1<0)ans += flag;    //up

            if(i+1<square_len){    //down
                int gid1 = gid+square_len;
                int id1 = grid_asses[gid1];
                if(id1<num){
                    ans += flag^1;
                }else ans += flag;
            }else ans += flag;

            if(j-1<0)ans += flag;    //left

            if(j+1<square_len){    //right
                int gid1 = gid+1;
                int id1 = grid_asses[gid1];
                if(id1<num){
                    ans += flag^1;
                }else ans += flag;
            }else ans += flag;
        }
    }

    return ans;
}

//count the edges of grid "gid" with each clusters
void checkEdgeSingleForLabel(
int labels_cnt[], const int &gid,
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    int i = gid/square_len;
    int j = gid%square_len;

    if(i-1>=0){    // up
        int gid1 = gid-square_len;
        int id1 = grid_asses[gid1];
        if(id1<num){
            int lb = cluster_labels[id1];
            labels_cnt[lb] += 1;
        }else labels_cnt[maxLabel] += 1;
    }else labels_cnt[maxLabel] += 1;

    if(i+1<square_len){    //down
        int gid1 = gid+square_len;
        int id1 = grid_asses[gid1];
        if(id1<num){
            int lb = cluster_labels[id1];
            labels_cnt[lb] += 1;
        }else labels_cnt[maxLabel] += 1;
    }else labels_cnt[maxLabel] += 1;

    if(j-1>=0){    //left
        int gid1 = gid-1;
        int id1 = grid_asses[gid1];
        if(id1<num){
            int lb = cluster_labels[id1];
            labels_cnt[lb] += 1;
        }else labels_cnt[maxLabel] += 1;
    }else labels_cnt[maxLabel] += 1;

    if(j+1<square_len){    //right
        int gid1 = gid+1;
        int id1 = grid_asses[gid1];
        if(id1<num){
            int lb = cluster_labels[id1];
            labels_cnt[lb] += 1;
        }else labels_cnt[maxLabel] += 1;
    }else labels_cnt[maxLabel] += 1;

    return;
}

//count edges of every grid with the cluster of label "lb"
void checkEdgeArrayForSingleLabel(
int E_grid[], const int &lb,
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    for(int i=0;i<N;i++){
        E_grid[i] = 0;
    }

    // #pragma omp parallel for num_threads(THREADS_NUM)
    for(int i=0;i<square_len;i++){
        int bias0 = i*square_len;
        for(int j=0;j<square_len;j++){
            int gid = bias0+j;

            if(i-1>=0){
                int gid1 = gid-square_len;
                int id1 = grid_asses[gid1];
                if((id1<num)&&(cluster_labels[id1]==lb)){
                    E_grid[gid] += 1;
                }
            }

            if(i+1<square_len){
                int gid1 = gid+square_len;
                int id1 = grid_asses[gid1];
                if((id1<num)&&(cluster_labels[id1]==lb)){
                    E_grid[gid] += 1;
                }
            }

            if(j-1>=0){
                int gid1 = gid-1;
                int id1 = grid_asses[gid1];
                if((id1<num)&&(cluster_labels[id1]==lb)){
                    E_grid[gid] += 1;
                }
            }

            if(j+1<square_len){
                int gid1 = gid+1;
                int id1 = grid_asses[gid1];
                if((id1<num)&&(cluster_labels[id1]==lb)){
                    E_grid[gid] += 1;
                }
            }
        }
    }

    return;
}

//count edges of every grid with different cluster from it
void checkEdgeArray(
int E_grid[],
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    for(int i=0;i<N;i++){
        E_grid[i]= 0;
    }

    // #pragma omp parallel for num_threads(THREADS_NUM)
    for(int i=0;i<square_len;i++){
        int bias0 = i*square_len;
        for(int j=0;j<square_len;j++){
            int gid = bias0+j;
            int id = grid_asses[gid];
            int lb = maxLabel;
            if(id<num)lb = cluster_labels[id];

            if(i-1>=0){
                int gid1 = gid-square_len;
                int id1 = grid_asses[gid1];
                int lb1 = maxLabel;
                if(id1<num)lb1 = cluster_labels[id1];
                if(lb1!=lb){
                    E_grid[gid] += 1;
                }
            }else E_grid[gid] += 1;

            if(i+1<square_len){
                int gid1 = gid+square_len;
                int id1 = grid_asses[gid1];
                int lb1 = maxLabel;
                if(id1<num)lb1 = cluster_labels[id1];
                if(lb1!=lb){
                    E_grid[gid] += 1;
                }
            }else E_grid[gid] += 1;

            if(j-1>=0){
                int gid1 = gid-1;
                int id1 = grid_asses[gid1];
                int lb1 = maxLabel;
                if(id1<num)lb1 = cluster_labels[id1];
                if(lb1!=lb){
                    E_grid[gid] += 1;
                }
            }else E_grid[gid] += 1;

            if(j+1<square_len){
                int gid1 = gid+1;
                int id1 = grid_asses[gid1];
                int lb1 = maxLabel;
                if(id1<num)lb1 = cluster_labels[id1];
                if(lb1!=lb){
                    E_grid[gid] += 1;
                }
            }else E_grid[gid] += 1;
        }
    }

    return;
}

//count edges of every grid with different cluster from it
void checkEdgeArray4(
int E_grid[][4],
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    for(int i=0;i<N;i++){
        for(int j=0;j<4;j++)
            E_grid[i][j] = 0;
    }

    // #pragma omp parallel for num_threads(THREADS_NUM)
    for(int i=0;i<square_len;i++){
        int bias0 = i*square_len;
        for(int j=0;j<square_len;j++){
            int gid = bias0+j;
            int id = grid_asses[gid];
            int lb = maxLabel;
            if(id<num)lb = cluster_labels[id];

            if(i-1>=0){
                int gid1 = gid-square_len;
                int id1 = grid_asses[gid1];
                int lb1 = maxLabel;
                if(id1<num)lb1 = cluster_labels[id1];
                if(lb1!=lb){
                    E_grid[gid][0] += 1;
                }
            }else E_grid[gid][0] += 1;

            if(i+1<square_len){
                int gid1 = gid+square_len;
                int id1 = grid_asses[gid1];
                int lb1 = maxLabel;
                if(id1<num)lb1 = cluster_labels[id1];
                if(lb1!=lb){
                    E_grid[gid][1] += 1;
                }
            }else E_grid[gid][1] += 1;

            if(j-1>=0){
                int gid1 = gid-1;
                int id1 = grid_asses[gid1];
                int lb1 = maxLabel;
                if(id1<num)lb1 = cluster_labels[id1];
                if(lb1!=lb){
                    E_grid[gid][2] += 1;
                }
            }else E_grid[gid][2] += 1;

            if(j+1<square_len){
                int gid1 = gid+1;
                int id1 = grid_asses[gid1];
                int lb1 = maxLabel;
                if(id1<num)lb1 = cluster_labels[id1];
                if(lb1!=lb){
                    E_grid[gid][3] += 1;
                }
            }else E_grid[gid][3] += 1;
        }
    }

    return;
}

std::vector<double> find_alpha_search(double delta1, double delta2, double delta3) {
//    printf("delta %.6lf %.6lf %.6lf\n",delta1, delta2, delta3);
    int order1=1, order2=2, order3=3;
    if(delta2*delta3>=0.00001){
    }else if(delta1*delta2>=0.00001){
        order1=3; order3=1;
        double tmp=delta1; delta1=delta3; delta3=tmp;
    }else if(delta1*delta3>=0.00001){
        order1=2; order2=1;
        double tmp=delta1; delta1=delta2; delta2=tmp;
    }

    if(std::abs(delta1)+std::abs(delta2)+std::abs(delta3)<0.0001) {
        std::vector<double> ret(2,-1);
        return ret;
    }
    double a = delta2 - delta1;
    double b = delta3 - delta1;
    double c = delta1;
    double best=std::abs(c),best1=0,best2=0;
    for(double i=0;i<=1;i+=0.01)
    for(double j=0;j<=1;j+=0.01){
        if(i+j>1)continue;
        if(std::abs(i*a+j*b+c)<best) {
            best = std::abs(i*a+j*b+c);
            best1 = i;
            best2 = j;
        }
    }
    std::vector<double> ret(2,0);
    ret[0] = best1;
    ret[1] = best2;
    if(a*a+b*b>0.0001) {
        double sq = std::sqrt(a*a+b*b);
        double d1 = b/sq;
        double d2 = -a/sq;
        for(double k=-1;k<=1;k+=0.01){
            double i = best1+k*d1;
            double j = best2+k*d2;
            if((i<0)||(i>1)||(j<0)||(j>1))continue;
            if(i+j>1)continue;
            if(std::abs(std::abs(j*delta3)-std::abs(i*delta2))<std::abs(std::abs(ret[1]*delta3)-std::abs(ret[0]*delta2))){
                ret[0] = i;
                ret[1] = j;
            }
        }
    }
    if(order1==2)ret[0]=1-ret[1]-ret[0];
    if(order1==3)ret[1]=1-ret[1]-ret[0];
    return ret;
}


std::vector<double> find_alpha_force(double delta1, double delta2, double delta3) {
//    printf("delta %.6lf %.6lf %.6lf\n",delta1, delta2, delta3);
    int order1=1, order2=2, order3=3;
    if(delta2*delta3>=0){
    }else if(delta1*delta2>=0){
        order1=3; order3=1;
        double tmp=delta1; delta1=delta3; delta3=tmp;
    }else if(delta1*delta3>=0){
        order1=2; order2=1;
        double tmp=delta1; delta1=delta2; delta2=tmp;
    }

    if(std::abs(delta1)+std::abs(delta2)+std::abs(delta3)<0.0001) {
        std::vector<double> ret(2,-1);
        return ret;
    }
    double a = delta2 - delta1;
    double b = delta3 - delta1;
    double c = delta1;
    double best=std::abs(c),best1=0,best2=0;
    for(double i=0;i<=1;i+=0.01)
    for(double j=0;j<=1;j+=0.01){
        if(i+j>1)continue;
        if(std::abs(i*a+j*b+c)<best) {
            best = std::abs(i*a+j*b+c);
            best1 = i;
            best2 = j;
        }
    }
    std::vector<double> ret(2,0);
    ret[0] = best1;
    ret[1] = best2;
    if(a*a+b*b>0.0001) {
        double sq = std::sqrt(a*a+b*b);
        double d1 = b/sq;
        double d2 = -a/sq;
        for(double k=-1;k<=1;k+=0.01){
            double i = best1+k*d1;
            double j = best2+k*d2;
            if((i<0)||(i>1)||(j<0)||(j>1))continue;
            if(i+j>1)continue;
            if(std::abs(j-i)<std::abs(ret[1]-ret[0])){
                ret[0] = i;
                ret[1] = j;
            }
        }
    }
    if(order1==2)ret[0]=1-ret[1]-ret[0];
    if(order1==3)ret[1]=1-ret[1]-ret[0];
    return ret;
}

double find_gamma(double delta1, double delta2) {
    double gamma=0;
    if(delta1*delta2>=delta1*delta1)gamma = 1;
    else if(delta1*delta2>=delta2*delta2)gamma = 0;
    else gamma = (delta2-delta1)*delta2/(delta1-delta2)/(delta1-delta2);
    return gamma;
}

std::vector<double> find_alpha(double delta1, double delta2, double delta3) {
//    printf("delta %.6lf %.6lf %.6lf\n",delta1, delta2, delta3);
    delta1 /= 1;
    delta2 /= 1.5;
    delta3 /= 20;
    if(delta1*delta1+delta2*delta2+delta3*delta3==0){
        std::vector<double>(2, -1);
    }
    if((delta1*delta2>=0)&&(delta2*delta3>=0)) {
        std::vector<double> ret(2, 0);
        ret[0] = delta1*delta3/(delta1*delta2+delta1*delta3+delta2*delta3);
        ret[1] = delta1*delta2/(delta1*delta2+delta1*delta3+delta2*delta3);
        return ret;
    }
    double delta[3];
    double M[3][3];
    delta[0] = delta1; delta[1] = delta2; delta[2] = delta3;
    std::vector<double> alpha(3, 1.0/3);
    for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)M[i][j] = delta[i]*delta[j];
    for(int it=0;it<15;it++){
        int best_t=-1;
        double best_score=-1;
        for(int r=0;r<3;r++){
            double now_score=0;
            for(int t=0;t<3;t++)now_score += alpha[t]*M[r][t];
            if((best_t==-1)||(now_score<best_score)){
                best_t = r;
                best_score = now_score;
            }
        }
        double d1=delta[best_t], d2=0;
        for(int i=0;i<3;i++)d2 += alpha[i]*delta[i];
        double gamma = find_gamma(d1, d2);
        for(int i=0;i<3;i++)alpha[i] *= 1-gamma;
        alpha[best_t] += gamma;
        for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        if((delta[i]==0)&&(delta[j]==0)) {
            alpha[i] = (alpha[i]+alpha[j])/2;
            alpha[j] = alpha[i];
        }
    }
    std::vector<double> ret(2, 0);
    ret[0] = alpha[1];
    ret[1] = alpha[2];
    return ret;
}

// calculate full cost
std::vector<double> checkCostForGlobal(
    const double Similar_cost_matrix[],
    const double Compact_cost_matrix[],
    const int grid_asses[], const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    const double &alpha, const double &beta) {

    double Convex_cost = 0;
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

int soft_choose(double a[], int n) {
    double *c = new double[n];
    double tot = 0;
    int last = 0;
    for(int i=0;i<n;i++) {
        double tmp = 0;
        if(a[i]>0) {
            tmp = a[i];
            last = i;
        }
        c[i] = tmp;
        tot = tot + tmp;
    }

    double rd = rand()/RAND_MAX*tot;
    int ans = last;
    for(int i=0;i<n;i++) {
        if(c[i]>0) {
            rd -= c[i];
            if(rd<=0){
                ans = i;
                break;
            }
        }
    }
    delete[] c;
    return ans;
}

int best_k_choose(double a[], int n, int k) {
    double *c = new double[n];
    int *o = new int[n];
    int last = 0;
    int cnt = 0;
    for(int i=0;i<n;i++) {
        o[i] = i;
        double tmp = 0;
        if(a[i]>0) {
            tmp = a[i];
            last = i;
            cnt += 1;
        }
        c[i] = tmp;
    }
    k = std::min(k, cnt);
    for(int i=0;i<n-1;i++)
    for(int j=i+1;j<n;j++)
    if(c[o[i]]<c[o[j]]){
        int tmp = o[i];
        o[i] = o[j];
        o[j] = tmp;
    }

    int rd = rand()%k;

    delete[] c;
    return o[rd];
}

#endif