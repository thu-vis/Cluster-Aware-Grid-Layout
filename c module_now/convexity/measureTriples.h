#ifndef _MEASURE_TRIPLES_H
#define _MEASURE_TRIPLES_H

#include "../utils/base.h"
#include "../utils/convexHull.h"
#include "../utils/util.h"

int (*global_triples)[4];   //saved triples
int *global_triples_head;   //heads of link list saving triples of every grid
int (*global_triples_list)[2];   //link list saving triples of every grid
int global_cnt = -1;    //count of triples
int global_N = -1;    //grid size of triples saved now

//check and store triples in an oblique line segment
void checkTriplesOfLine2(
int triples[][4],
int triples_head[],
int triples_list[][2],
int &cnt,
const int &square_len,
const int &gid1, const int &gid2,
int oriX, int oriY, int targetX, int startX, int dx, int dy){
    // printf("line %d %d %d %d\n", oriX, oriY, targetX, startX);
    if((dx==0)||(dy==0))return;
    int nx = oriX;
    int ny = oriY;
    int x = startX;
    while(nx<targetX){
        int new_nx = min((x+1)*2, targetX);
        int new_ny = (new_nx-nx)*dy+ny;
        int bt = min(ny, new_ny);
        int tp = max(ny, new_ny);
        int y = bt/(2*dx);
        int ty = y*2*dx;
        // printf("line x y %d %d %d %d\n", x, y, bt, tp);
        while(ty<tp){
            int gid=x*square_len+y;
            // printf("line y %d %d\n", y, ty);
            // printf("line gid %d %d %d\n", gid1, gid, gid2);
            if((gid!=gid1)&&(gid!=gid2)){
                triples[cnt][0] = gid1;
                triples[cnt][1] = gid;
                triples[cnt][2] = gid2;
                triples[cnt][3] = 1;

                triples_list[cnt*3][0] = cnt;
                triples_list[cnt*3][1] = triples_head[gid1];
                triples_head[gid1] = cnt*3;
                triples_list[cnt*3+1][0] = cnt;
                triples_list[cnt*3+1][1] = triples_head[gid];
                triples_head[gid] = cnt*3+1;
                triples_list[cnt*3+2][0] = cnt;
                triples_list[cnt*3+2][1] = triples_head[gid2];
                triples_head[gid2] = cnt*3+2;
                cnt += 1;
            }
            y += 1;
            ty += 2*dx;
        }
        nx = new_nx;
        ny = new_ny;
        x += 1;
    }
    return;
}

//check and store triples in an line segment
void checkTriplesOfLine(
int triples[][4],
int triples_head[],
int triples_list[][2],
int &cnt,
const int &square_len,
const int &gid1, const int &gid2){
    int startX = gid1/square_len;
    int startY = gid1%square_len;
    int endX = gid2/square_len;
    int endY = gid2%square_len;
    if(startX>endX) {
        std::swap(startX, endX);
        std::swap(startY, endY);
    }
    if(startX==endX) {
        int bt = min(startY, endY);
        int tp = max(startY, endY);
        int x = startX;
        int y = bt + 1;
        while(y<tp){
            int gid = x*square_len+y;
            triples[cnt][0] = gid1;
            triples[cnt][1] = gid;
            triples[cnt][2] = gid2;
            triples[cnt][3] = 1;

            triples_list[cnt*3][0] = cnt;
            triples_list[cnt*3][1] = triples_head[gid1];
            triples_head[gid1] = cnt*3;
            triples_list[cnt*3+1][0] = cnt;
            triples_list[cnt*3+1][1] = triples_head[gid];
            triples_head[gid] = cnt*3+1;
            triples_list[cnt*3+2][0] = cnt;
            triples_list[cnt*3+2][1] = triples_head[gid2];
            triples_head[gid2] = cnt*3+2;
            cnt += 1;
            y += 1;
        }
    }else if(startY==endY) {
        int bt = min(startX, endX);
        int tp = max(startX, endX);
        int y = startY;
        int x = bt + 1;
        while(x<tp){
            int gid = x*square_len+y;
            triples[cnt][0] = gid1;
            triples[cnt][1] = gid;
            triples[cnt][2] = gid2;
            triples[cnt][3] = 1;
            triples_list[cnt*3][0] = cnt;
            triples_list[cnt*3][1] = triples_head[gid1];
            triples_head[gid1] = cnt*3;
            triples_list[cnt*3+1][0] = cnt;
            triples_list[cnt*3+1][1] = triples_head[gid];
            triples_head[gid] = cnt*3+1;
            triples_list[cnt*3+2][0] = cnt;
            triples_list[cnt*3+2][1] = triples_head[gid2];
            triples_head[gid2] = cnt*3+2;
            cnt += 1;
            x += 1;
        }
    }else {
        int dx = endX - startX;
        int dy = endY - startY;
        int oriX = startX*2+1;
        int oriY = (startY*2+1)*dx;
        int targetX = endX*2+1;
        checkTriplesOfLine2(triples, triples_head, triples_list,
        cnt, square_len, gid1, gid2, oriX, oriY, targetX, startX, dx, dy);
    }
}

//check and save triples
int getTriples(
int triples[][4],
int triples_head[],
int triples_list[][2],
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel){
    int cnt=0;
    for(int gid=0;gid<N;gid++)triples_head[gid] = -1;

    for(int gid1=0;gid1<N-1;gid1++)
    for(int gid2=gid1+1;gid2<N;gid2++)checkTriplesOfLine(triples, triples_head, triples_list, cnt, square_len, gid1, gid2);

    return cnt;
}

//check triples of grid "gid" by the saved preprocess-information
std::vector<int> checkTriplesByDict(
int innerDict[], int outerDict[],
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
const int gid, const int lb) {
    int fcnt = 0;
    int ecnt = 0;
    for(int lb1=0;lb1<maxLabel;lb1++) {
        fcnt += innerDict[gid*maxLabel+lb1];
        if(lb!=lb1)ecnt += innerDict[gid*maxLabel+lb1];
        if(lb==lb1){
            for(int lb2=0;lb2<2;lb2++){
                fcnt += outerDict[gid*maxLabel*2+lb1*2+lb2];
                if(lb2==1)ecnt += outerDict[gid*maxLabel*2+lb1*2+lb2];
            }
        }
    }
    std::vector<int> pair(2, 0);
    pair[0] = ecnt;
    pair[1] = fcnt;
    return pair;
}

//preprocess for triples of every gid
void updateTripleDict(
int innerDict[], int outerDict[],
const int triples[][4],
const int triples_head[],
const int triples_list[][2],
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
const int gid) {
    for(int ii=triples_head[gid];ii>=0;ii=triples_list[ii][1])
    {
        int i = triples_list[ii][0];
        int gid1 = triples[i][0];
        int gid0 = triples[i][1];
        int gid2 = triples[i][2];
        int times = triples[i][3];
        if(gid==gid0){
            int id1 = grid_asses[gid1], id2 = grid_asses[gid2];
            if((id1>=num)||(id2>=num))continue;
            if(cluster_labels[id1]==cluster_labels[id2]){
                int lb = cluster_labels[id1];
                innerDict[gid*maxLabel+lb] += times;
            }
        }else if(gid1==gid){
            int id0 = grid_asses[gid0], id2 = grid_asses[gid2];
            if(id2>=num)continue;
            int lb0 = 0;
            if((id0>=num)||(cluster_labels[id0]!=cluster_labels[id2]))
                lb0 = 1;
            int lb2 = cluster_labels[id2];
            outerDict[gid*maxLabel*2+lb2*2+lb0] += times;
        }else if(gid2==gid){
            int id0 = grid_asses[gid0], id1 = grid_asses[gid1];
            if(id1>=num)continue;
            int lb1 = cluster_labels[id1];
            int lb0 = 0;
            if((id0>=num)||(cluster_labels[id0]!=cluster_labels[id1]))
                lb0 = 1;
            outerDict[gid*maxLabel*2+lb1*2+lb0] += times;
        }
    }
    return;
}

//preprocess for triples of every gid with saved information of old layout
void updateTripleDictwithOld(
double innerDict[], double outerDict[],
const int change_triples[][3],
const double change_triples_times[],
const int change_triples_head[],
const int change_triples_list[][2],
const int grid_asses[],
const int old_grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
const int gid) {
    // printf("start update %d\n", gid);
    // printf("change triples head %d\n", change_triples_head[gid]);
    for(int ii=change_triples_head[gid];ii>=0;ii=change_triples_list[ii][1])
    {
        // printf("ii %d\n", ii);
        int i = change_triples_list[ii][0];
        // printf("update %d %d\n", gid, i);
        int gid1 = change_triples[i][0];
        int gid0 = change_triples[i][1];
        int gid2 = change_triples[i][2];
        double times = change_triples_times[i];
        // printf("gid %d %d %d %.2lf\n", gid1, gid0, gid2, times);
        if(gid==gid0){
            int id1 = grid_asses[gid1], id2 = grid_asses[gid2];
            if((id1>=num)||(id2>=num))continue;
            if(cluster_labels[id1]==cluster_labels[id2]){
                int lb = cluster_labels[id1];
                innerDict[gid*maxLabel+lb] += times;
            }
        }else if(gid1==gid){
            int id0 = grid_asses[gid0], id2 = grid_asses[gid2];
            if(id2>=num)continue;
            int lb2 = cluster_labels[id2];
            int lb0 = 0;
            if((id0>=num)||(cluster_labels[id0]!=cluster_labels[id2]))
                lb0 = 1;
            outerDict[gid*maxLabel*2+lb2*2+lb0] += times;
        }else if(gid2==gid){
            int id0 = grid_asses[gid0], id1 = grid_asses[gid1];
            if(id1>=num)continue;
            int lb1 = cluster_labels[id1];
            int lb0 = 0;
            if((id0>=num)||(cluster_labels[id0]!=cluster_labels[id1]))
                lb0 = 1;
            outerDict[gid*maxLabel*2+lb1*2+lb0] += times;
        }

        if(gid==gid0){
            int id1 = old_grid_asses[gid1], id2 = old_grid_asses[gid2];
            if((id1>=num)||(id2>=num))continue;
            if(cluster_labels[id1]==cluster_labels[id2]){
                int lb = cluster_labels[id1];
                innerDict[gid*maxLabel+lb] -= times;
            }
        }else if(gid1==gid){
            int id0 = old_grid_asses[gid0], id2 = old_grid_asses[gid2];
            if(id2>=num)continue;
            int lb2 = cluster_labels[id2];
            int lb0 = 0;
            if((id0>=num)||(cluster_labels[id0]!=cluster_labels[id2]))
                lb0 = 1;
            outerDict[gid*maxLabel*2+lb2*2+lb0] -= times;
        }else if(gid2==gid){
            int id0 = old_grid_asses[gid0], id1 = old_grid_asses[gid1];
            if(id1>=num)continue;
            int lb1 = cluster_labels[id1];
            int lb0 = 0;
            if((id0>=num)||(cluster_labels[id0]!=cluster_labels[id1]))
                lb0 = 1;
            outerDict[gid*maxLabel*2+lb1*2+lb0] -= times;
        }
        // printf("done update %d %d\n", gid, i);
    }
    return;
}

//check convexity of triples
//return T[0..1]: convexity(in the form of loss/cost) pair in layout now
std::vector<double> checkConvexForTArray(
const int triples[][4],
const int &cnt,
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    double T0 = 0;
    double T1 = 0;

    #pragma omp parallel for reduction(+:T0, T1) num_threads(THREADS_NUM)
    for(int starti=0;starti<N;starti++){
        for(int i=starti;i<cnt;i+=N){
            int gid1 = triples[i][0];
            int gid0 = triples[i][1];
            int gid2 = triples[i][2];
            int times = triples[i][3];

            // printf("triple %d %d %d %d\n", gid1, gid0, gid2, times);

            int id1 = grid_asses[gid1];
            int id0 = grid_asses[gid0];
            int id2 = grid_asses[gid2];
            if((id1>=num)||(id2>=num))continue;
            int lb1 = cluster_labels[id1];
            int lb2 = cluster_labels[id2];
            int lb0 = -1;
            if(id0<num)lb0 = cluster_labels[id0];
            if(lb1==lb2){
                T1 += times;
                if(lb1!=lb0)T0 += times;
            }
        }
    }

    std::vector<double> T_pair(2, 0);
    T_pair[0] = T0;
    T_pair[1] = T1;

    return T_pair;
}

//check convexity of triples, python interface
std::vector<double> checkConvexForT(
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

    int (*triples)[4];
    int *triples_head;
    int (*triples_list)[2];
    int cnt;
    if(global_N==N){    //use saved triples link list
        triples = global_triples;
        triples_head = global_triples_head;
        triples_list = global_triples_list;
        cnt = global_cnt;
    }else{    //calculate and store new link list
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

    std::vector<double> ret = checkConvexForTArray(triples, cnt, grid_asses, cluster_labels, N, num, square_len, maxLabel);

    delete[] grid_asses;
    delete[] cluster_labels;
    return ret;
}

//check convexity of triples with saved information of old layout
//old_T_pair[0..1]: convexity(in the form of loss/cost) pair in old grid layout
//return T[0..1]: convexity(in the form of loss/cost) pair in layout now
std::vector<double> checkConvexForTwithOld(
int triples[][4],
int triples_head[],
int triples_list[][2],
double old_T_pair[],
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
    #pragma omp parallel for reduction(+:T0, T1, dec_T0, dec_T1) num_threads(THREADS_NUM)
    for(int gid=0;gid<N;gid++)
    if(changed[gid]==1) {
        for(int ii=triples_head[gid];ii>=0;ii=triples_list[ii][1]){
            int i = triples_list[ii][0];
            int gid1 = triples[i][0];
            int gid0 = triples[i][1];
            int gid2 = triples[i][2];
            double times = triples[i][3];
            times = (1.0*times)/(changed[gid1]+changed[gid0]+changed[gid2]);

            int id1 = grid_asses[gid1];
            int id0 = grid_asses[gid0];
            int id2 = grid_asses[gid2];
            if((id1<num)&&(id2<num)) {
                int lb1 = cluster_labels[id1];
                int lb2 = cluster_labels[id2];
                int lb0 = -1;
                if(id0<num)lb0 = cluster_labels[id0];
                if(lb1==lb2){
                    T1 += times;
                    if(lb1!=lb0)T0 += times;
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
                    dec_T1 += times;
                    if(lb1!=lb0)dec_T0 += times;
                }
            }
        }
    }
    std::vector<double> T_pair(2, 0);
    T_pair[0] = old_T_pair[0]+T0-dec_T0;
    T_pair[1] = old_T_pair[1]+T1-dec_T1;

    delete[] changed;
    return T_pair;
}

// count changed triples(different from old layout)
int getChangeTriplesCnt(
const int cnt,
int triples[][4],
int triples_head[],
int triples_list[][2],
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

    int change_cnt = 0;
    #pragma omp parallel for reduction(+:change_cnt) num_threads(THREADS_NUM)
    for(int gid=0;gid<N;gid++) {
        if(changed[gid]==1){
            for(int ii=triples_head[gid];ii>=0;ii=triples_list[ii][1])change_cnt += 1;
        }
    }

    delete[] changed;
    return change_cnt;
}

// check and store changed triples(different from old layout)
// change_triples_xxx: link list to save changed triples of every grid
// return int: number of changed triples
int getChangeTriples(
int change_triples[][3],
double change_triples_times[],
int change_triples_head[],
int change_triples_list[][2],
const int cnt,
int triples[][4],
int triples_head[],
int triples_list[][2],
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

    int change_cnt = 0;
    for(int gid=0;gid<N;gid++)change_triples_head[gid] = -1;

    for(int gid=0;gid<N;gid++)
    if(changed[gid]==1) {
        for(int ii=triples_head[gid];ii>=0;ii=triples_list[ii][1]){
            int i = triples_list[ii][0];
            int gid1 = triples[i][0];
            int gid0 = triples[i][1];
            int gid2 = triples[i][2];
            int times = triples[i][3];
            change_triples[change_cnt][0] = gid1;
            change_triples[change_cnt][1] = gid0;
            change_triples[change_cnt][2] = gid2;
            change_triples_times[change_cnt] = (1.0*times)/(changed[gid1]+changed[gid0]+changed[gid2]);
            change_triples_list[change_cnt*3][0] = change_cnt;
            change_triples_list[change_cnt*3][1] = change_triples_head[gid1];
            change_triples_head[gid1] = change_cnt*3;
            change_triples_list[change_cnt*3+1][0] = change_cnt;
            change_triples_list[change_cnt*3+1][1] = change_triples_head[gid0];
            change_triples_head[gid0] = change_cnt*3+1;
            change_triples_list[change_cnt*3+2][0] = change_cnt;
            change_triples_list[change_cnt*3+2][1] = change_triples_head[gid2];
            change_triples_head[gid2] = change_cnt*3+2;
            change_cnt += 1;
        }
    }

    delete[] changed;
    return change_cnt;
}


//get convexity cost matrix of triples
//cost_matrix_a[0..N*N-1]: cost matrix will be stored in it
//save: if save calculated intermediate information
//load: if use saved information to accelerate
//old_innerDict, old_outerDict: intermediate information
//old_grid_asses: old grid layout
//now_T_pair: convexity(in the form of loss/cost) pair to save and load
void getCostMatrixForTArrayToArray(
int grid_asses[],
int cluster_labels[],
double cost_matrix_a[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
bool save=false, bool load=false, int old_innerDict[]=nullptr, int old_outerDict[]=nullptr,
int old_grid_asses[]=nullptr, double now_T_pair[]=nullptr) {

    int (*triples)[4];
    int *triples_head;
    int (*triples_list)[2];
    int cnt;

    int (*change_triples)[3];
    double *change_triples_times;
    int *change_triples_head;
    int (*change_triples_list)[2];
    int change_cnt;

    if(global_N==N){    //use saved triples link list
        triples = global_triples;
        triples_head = global_triples_head;
        triples_list = global_triples_list;
        cnt = global_cnt;
    }else{    //calculate and store new link list
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

    if(load){
        // change_cnt = getChangeTriplesCnt(cnt, triples, triples_head, triples_list, grid_asses, old_grid_asses, cluster_labels, N, num, square_len, maxLabel);
        change_cnt = cnt;
        change_triples = new int[change_cnt][3];
        change_triples_times = new double[change_cnt];
        change_triples_head = new int[N];
        change_triples_list = new int[change_cnt*3][2];
        change_cnt = getChangeTriples(change_triples, change_triples_times, change_triples_head, change_triples_list,
        cnt, triples, triples_head, triples_list, grid_asses, old_grid_asses, cluster_labels, N, num, square_len, maxLabel);
//        printf("cnt %d %d\n", cnt, change_cnt);
    }

    std::vector<double> T_pair1(2, 0);
    if(!load) {
        T_pair1 = checkConvexForTArray(triples, cnt, grid_asses, cluster_labels, N, num, square_len, maxLabel);
    }else {
        T_pair1[0] = now_T_pair[0];
        T_pair1[1] = now_T_pair[1];
    }

    // printf("convex checked\n");

    double correct1=0, full1=0;
    full1 = T_pair1[1];
    correct1 = T_pair1[1]-T_pair1[0];

    int *innerDict = new int[N*maxLabel];
    int *outerDict = new int[N*maxLabel*2];
    for(int i=0;i<N*maxLabel;i++)innerDict[i] = 0;
    for(int i=0;i<N*maxLabel*2;i++)outerDict[i] = 0;

    if(!load) {
        #pragma omp parallel for num_threads(THREADS_NUM)
        for(int gid=0;gid<N;gid++)    //preprocess for triples of every gid
            updateTripleDict(innerDict, outerDict, triples, triples_head, triples_list, grid_asses, cluster_labels, N, num, square_len, maxLabel, gid);
    }else {
        double *innerDict_d = new double[N*maxLabel];
        double *outerDict_d = new double[N*maxLabel*2];
        for(int i=0;i<N*maxLabel;i++)innerDict_d[i] = old_innerDict[i];
        for(int i=0;i<N*maxLabel*2;i++)outerDict_d[i] = old_outerDict[i];
        #pragma omp parallel for num_threads(THREADS_NUM)
        for(int gid=0;gid<N;gid++)    //preprocess for triples of every gid
            updateTripleDictwithOld(innerDict_d, outerDict_d,
            change_triples, change_triples_times, change_triples_head, change_triples_list,
            grid_asses, old_grid_asses, cluster_labels, N, num, square_len, maxLabel, gid);
        for(int i=0;i<N*maxLabel;i++)innerDict[i] = round(innerDict_d[i]);
        for(int i=0;i<N*maxLabel*2;i++)outerDict[i] = round(outerDict_d[i]);
        delete[] innerDict_d;
        delete[] outerDict_d;
    }

    if(save) {
        for(int i=0;i<N*maxLabel;i++)old_innerDict[i] = innerDict[i];
        for(int i=0;i<N*maxLabel*2;i++)old_outerDict[i] = outerDict[i];
    }

    // printf("dict updated\n");

    double *correct_matrix1 = new double[N];    //triples of every grid in layout now
    double *full_matrix1 = new double[N];

    for(int gid=0;gid<N;gid++){
        int lb = maxLabel;
        int id = grid_asses[gid];
        if(id<num)lb = cluster_labels[id];
        std::vector<int> pair = checkTriplesByDict(innerDict, outerDict, grid_asses, cluster_labels, N, num, square_len, maxLabel, gid, lb);
        correct_matrix1[gid] = pair[1]-pair[0];
        full_matrix1[gid] = pair[1];
    }

    double *correct_matrix2 = new double[N*(maxLabel+1)];    //elements of same label matching a same grid, have same cost
    double *full_matrix2 = new double[N*(maxLabel+1)];

    for(int i=0;i<N*(maxLabel+1);i++){
        correct_matrix2[i] = full_matrix2[i] = -1;
    }

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid2=0;gid2<N;gid2++){
        int bias = gid2*N;
        int bias2 = gid2*(maxLabel+1);
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
            if(full_matrix2[bias2+lb]<0){
                // int tmp = grid_asses[gid2];
                // grid_asses[gid2] = id1;
                std::vector<int> pair = checkTriplesByDict(innerDict, outerDict, grid_asses, cluster_labels, N, num, square_len, maxLabel, gid2, lb);
                correct_matrix2[bias2+lb] = pair[1]-pair[0];
                full_matrix2[bias2+lb] = pair[1];
                // grid_asses[gid2] = tmp;
            }

            double correct2 = correct1-correct_matrix1[gid2]+correct_matrix2[bias2+lb];
            double full2 = full1-full_matrix1[gid2]+full_matrix2[bias2+lb];
            // printf("%d %d %lf %lf %lf %lf\n", gid2, id1, correct2, full2, correct1, full1);
            if(correct2*full1>0){
                // cost_matrix[gid2][id1] = (1.0*(correct1-correct2)/correct2+1.0*(full2-full1)/full1)*N;
                cost_matrix_a[bias+id1] = (correct1*full2-correct2*full1)/(correct2*full1)*N;
            }else cost_matrix_a[bias+id1] = 0;
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

    // printf("cost matrix done\n");
    if(load){
        delete[] change_triples;
        delete[] change_triples_times;
        delete[] change_triples_head;
        delete[] change_triples_list;
    }
    delete[] innerDict;
    delete[] outerDict;
    delete[] correct_matrix1;
    delete[] full_matrix1;
    delete[] correct_matrix2;
    delete[] full_matrix2;
    return;
}

//get convexity cost matrix of triples, python interface
std::vector<std::vector<double>> getCostMatrixForT(
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels) {
    // printf("cost matrix start\n");
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
    getCostMatrixForTArrayToArray(grid_asses, cluster_labels, cost_matrix_a, N, num, square_len, maxLabel);

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid=0;gid<N;gid++) {
        int bias = gid*N;
        for(int id=0;id<num;id++){
            cost_matrix[gid][id] = cost_matrix_a[bias+id];
        }
    }

    delete[] grid_asses;
    delete[] cluster_labels;
    delete[] cost_matrix_a;
    return cost_matrix;
}

//calculate full cost, using convexity of triples
//alpha, beta: cost arguments
//return cost[0..3], full/Similar/Compact/Convex cost
std::vector<double> checkCostForT(
const double Similar_cost_matrix[],
const double Compact_cost_matrix[],
const int grid_asses[], const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel,
const double &alpha, const double &beta,
bool save=false, bool load=false, int old_grid_asses[]=nullptr, double old_T_pair[]=nullptr) {

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
        T_pair = checkConvexForTArray(triples, cnt, grid_asses, cluster_labels, N, num, square_len, maxLabel);
        // printf("T_pair %.2lf %.2lf\n", T_pair[0], T_pair[1]);
    }else {
        T_pair = checkConvexForTwithOld(triples, triples_head, triples_list,
        old_T_pair, grid_asses, old_grid_asses, cluster_labels, N, num, square_len, maxLabel);
    }

    if(save){
        old_T_pair[0] = T_pair[0];
        old_T_pair[1] = T_pair[1];
    }

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

#endif