#include <iostream>
#include <vector>
#include <utility>
#include <ctime>
#include <algorithm>
#include <math.h>

int checkTriple(
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num,
const int &gid1, const int &gid0, const int &gid2){
    int id1 = grid_asses[gid1];
    int id0 = grid_asses[gid0];
    int id2 = grid_asses[gid2];
    if((id1>=num)||(id2>=num))return -1;
    int lb1 = cluster_labels[id1];
    int lb2 = cluster_labels[id2];
    int lb0 = -1;
    if(id0<num)lb0 = cluster_labels[id0];
    if(lb1!=lb2)return -1;
    if(lb1!=lb0)return 1;
    return 0;
}

//扫描边界三元组 part3
int checkBorderTriplesOfLine3(
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
const int &gid1, const int &gid2,
int oriX, int oriY, int targetX, int startX, int dx, int dy){
    if((dx==0)||(dy==0))return 1;
    int nx = oriX;
    int ny = oriY;
    int x = startX;
    while(nx<targetX){
        int new_nx = std::min((x+1)*2, targetX);
        int new_ny = (new_nx-nx)*dy+ny;
        int bt = std::min(ny, new_ny);
        int tp = std::max(ny, new_ny);
        int y = bt/(2*dx);
        int ty = y*2*dx;
        while(ty<tp){
            int gid_mid = x*square_len+y;
            if(checkTriple(grid_asses, cluster_labels, N, num, gid1, gid_mid, gid2)==1)
                return 0;
            y += 1;
            ty += 2*dx;
        }
        nx = new_nx;
        ny = new_ny;
        x += 1;
    }
    return 1;
}

//扫描边界三元组 part2
int checkBorderTriplesOfLine2(
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
int startX, int startY, int endX, int endY,
const int &gid1, const int &gid2){
    if(startX>endX) {
        std::swap(startX, endX);
        std::swap(startY, endY);
    }
    if(startX==endX) {
        int bt = std::min(startY, endY)/2;
        int tp = (std::max(startY, endY)+1)/2;
        int y = bt;
        while(y<tp){
            if(startX%2==1) {
                int gid_mid = (startX/2)*square_len+y;
                if(checkTriple(grid_asses, cluster_labels, N, num, gid1, gid_mid, gid2)==1)
                    return 0;
            }
            else{
                int gid_mid1 = (startX/2)*square_len+y;
                int gid_mid2 = (startX/2-1)*square_len+y;
                if(((startX/2>=square_len)||(checkTriple(grid_asses, cluster_labels, N, num, gid1, gid_mid1, gid2)==1))
                 &&((startX/2-1<0)||(checkTriple(grid_asses, cluster_labels, N, num, gid1, gid_mid2, gid2)==1)))
                    return 0;
            }
            y += 1;
        }
        return 1;
    }else if(startY==endY) {
        int bt = std::min(startX, endX)/2;
        int tp = (std::max(startX, endX)+1)/2;
        int x = bt;
        while(x<tp){
            if(startY%2==1) {
                int gid_mid = x*square_len+(startY/2);
                if(checkTriple(grid_asses, cluster_labels, N, num, gid1, gid_mid, gid2)==1)
                return 0;
            }
            else{
                int gid_mid1 = x*square_len+(startY/2);
                int gid_mid2 = x*square_len+(startY/2-1);
                if(((startY/2>=square_len)||(checkTriple(grid_asses, cluster_labels, N, num, gid1, gid_mid1, gid2)==1))
                 &&((startY/2-1<0)||(checkTriple(grid_asses, cluster_labels, N, num, gid1, gid_mid2, gid2)==1)))
                    return 0;
            }
            x += 1;
        }
        return 1;
    }else {
        int dx = endX - startX;
        int dy = endY - startY;
        int oriX = startX;
        int oriY = startY*dx;
        int targetX = endX;
        return checkBorderTriplesOfLine3(grid_asses, cluster_labels, N, num, square_len, maxLabel, gid1, gid2, oriX, oriY, targetX, startX/2, dx, dy);
    }
    return 1;
}

//扫描边界三元组 part1
std::vector<int> checkBorderTriplesOfLine(
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
const int &gid1, const int &gid2){
    int x1 = gid1/square_len;
    int y1 = gid1%square_len;
    int x2 = gid2/square_len;
    int y2 = gid2%square_len;
    int lb = cluster_labels[grid_asses[gid1]];
    int fx[4] = {0, 0, -1, 1};
    int fy[4] = {1, -1, 0, 0};
    int T0=0, T1=0;
    for(int i=0;i<4;i++) {
        int nx1=x1+fx[i];
        int ny1=y1+fy[i];
        if((nx1<square_len)&&(nx1>=0)&&(ny1<square_len)&&(ny1>=0)) {
            int ngid1 = nx1*square_len+ny1;
            if((grid_asses[ngid1]<num)&&(cluster_labels[grid_asses[ngid1]]==lb))continue;
        }
        for(int j=0;j<4;j++) {
            int nx2=x2+fx[j];
            int ny2=y2+fy[j];
            if((nx2<square_len)&&(nx2>=0)&&(ny2<square_len)&&(ny2>=0)) {
                int ngid2 = nx2*square_len+ny2;
                if((grid_asses[ngid2]<num)&&(cluster_labels[grid_asses[ngid2]]==lb))continue;
            }
            int startX = x1*2+1+fx[i];
            int startY = y1*2+1+fy[i];
            int endX = x2*2+1+fx[j];
            int endY = y2*2+1+fy[j];
            int t = checkBorderTriplesOfLine2(grid_asses, cluster_labels, N, num, square_len, maxLabel, startX, startY, endX, endY, gid1, gid2);
            T0 += t;
            T1 += 1;
        }
    }
    std::vector<int> ret(2, 0);
    ret[0] = T0;
    ret[1] = T1;
    return ret;
}

//轮廓三元组方法，检查凸性
std::vector<double> checkConvexForTBArray(
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
bool save=false, bool load=false, double label_pairs[][2]=nullptr, int mainLabel1=-1, int mainLabel2=-1) {

    int *head = new int[maxLabel];
    int *last = new int[num];
    int *element_asses = new int[num];
    for(int i=0;i<maxLabel;i++)head[i] = -1;

    for(int gid=0;gid<N;gid++)
    if(grid_asses[gid]<num) {
        int id = grid_asses[gid];
        element_asses[id] = gid;
        int lb = cluster_labels[id];
        last[id] = head[lb]; head[lb] = id;
    }

    double T0 = 0, T1 = 0;
    for(int li=0;li<maxLabel;li++){
        if((load)&&(li!=mainLabel1)&&(li!=mainLabel2)){
            T0 += label_pairs[li][0];
            T1 += label_pairs[li][1];
            continue;
        }
        int cnt=0;
        double tmp_T0=0, tmp_T1=0;
        for(int id1=head[li];id1>=0;id1=last[id1]){
            int gid1 = element_asses[id1];
            int x1 = gid1/square_len;
            int y1 = gid1%square_len;
            for(int id2=head[li];id2>=0;id2=last[id2]){
                int gid2 = element_asses[id2];
                int x2 = gid2/square_len;
                int y2 = gid2%square_len;
                std::vector<int> pair = checkBorderTriplesOfLine(grid_asses, cluster_labels, N, num, square_len, maxLabel, gid1, gid2);
                tmp_T0 += pair[0];
                tmp_T1 += pair[1];
            }
        }
        T1 += tmp_T1;
        T0 += tmp_T0;
        if(save){
            label_pairs[li][0] = tmp_T0;
            label_pairs[li][1] = tmp_T1;
        }
    }
    std::vector<double> T_pair(2, 0);
    T_pair[0] = (T1-T0);
    T_pair[1] = T1;

    delete[] head;
    delete[] last;
    delete[] element_asses;
    return T_pair;
}

//轮廓三元组方法，计算整体代价（与原layout相似性、紧凑性、凸性）
std::vector<double> checkCostForTB(
    const double Similar_cost_matrix[],
    const double Compact_cost_matrix[],
    const int grid_asses[], const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    const double &alpha, const double &beta,
    bool save=false, bool load=false, double label_pairs[][2]=nullptr, int mainLabel1=-1, int mainLabel2=-1) {

    std::vector<double> T_pair = checkConvexForTBArray(grid_asses, cluster_labels, N, num, square_len, maxLabel, save, load, label_pairs, mainLabel1, mainLabel2);

    double correct=0, full=0;

    full = T_pair[1];
    correct = T_pair[1]-T_pair[0];

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