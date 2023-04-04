//example.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <utility>
#include <ctime>

#include "utils/base.h"
#include "utils/simpleCluster.h"
#include "utils/lap.h"
#include "utils/util.h"
#include "convexity/measureEdge.h"
#include "convexity/measureTriples.h"
#include "convexity/measureCHS.h"
#include "convexity/measureCHC.h"
#include "convexity/newMeasure2020.h"
#include "convexity/measureDoubles.h"
#include "convexity/newMeasureTB.h"
#include "convexity/newMeasureBoundary.h"
#include "convexity/newMeasureDeviation.h"
#include "convexity/newMeasureAlphaTriples.h"
#include "convexity/newMeasureMoS.h"

// simple cluster, python interface
std::vector<int> getClusters(
const std::vector<int> &_grid_asses,
const std::vector<int> &_labels) {
    int N = _grid_asses.size();
    int num = _labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = max(maxLabel, _labels[i]+1);

    int *grid_asses = new int[N];
    int *labels = new int[num];
    int *cluster_labels = new int[num];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
    for(int i=0;i<num;i++)labels[i] = _labels[i];

    getClustersArray(cluster_labels, maxLabel+5, int(0.02*N), grid_asses, labels, N, num, square_len);

    std::vector<int> ret(num, 0);
    for(int i=0;i<num;i++)ret[i] = cluster_labels[i];

    delete[] grid_asses;
    delete[] labels;
    delete[] cluster_labels;
    return ret;
}

// lap jv solve bi-graph-match
std::vector<int> solveLapArray(
    const double dist[],
    const int N,
    int k=50) {
    float *cost_matrix = new float[N*N];
    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int i1=0;i1<N;i1++){
        int bias = i1*N;
        for(int i2=0;i2<N;i2++){
            int i = bias+i2;
            cost_matrix[i] = dist[i]+0.5;   // prevent 0 distance
        }
    }
    float *u = new float[N];
    float *v = new float[N];
    int *grid_asses = new int[N];
    int *element_asses = new int[N];
    float cost = lap(N, cost_matrix, false, grid_asses, element_asses, u, v, min(k,N));


    std::vector<int> ret(N, 0);
    for(int i=0;i<N;i++)ret[i] = grid_asses[i];

    delete[] cost_matrix;
    delete[] u;
    delete[] v;
    delete[] grid_asses;
    delete[] element_asses;

    return ret;
}

//lap jv, python interface
std::vector<int> solveLap(
    const std::vector<std::vector<double>> &_dist,
    const bool use_knn=false, int k_value=50) {

    int n_row = _dist.size();
    int n_col = _dist[0].size();

    double *dist = new double[n_row*n_col];
    for(int i=0;i<n_row;i++)
    for(int j=0;j<n_col;j++)dist[i*n_col+j]=_dist[i][j];
    if(!use_knn)k_value = n_col;
    std::vector<int> ret = solveLapArray(dist, n_row, k_value);

    delete[] dist;
    return ret;
}

// KM, deal the augmenting path
void workback(int x,int lb[][3], int p2[])
{
    while(x>1){ int y=lb[x][0]; x=lb[x][2]; p2[y]=lb[x][1]; }
}

// KM
std::vector<int> solveKMArray(
    const double dist[],
    const int n_row, const int n_col) {
    const double tiny = 0.000001;

    double amax=0;
    for(int i=0;i<n_row;i++){
        int bias = i*n_col;
        for(int j=0;j<n_col;j++){
            amax=max(dist[bias+j],amax);
        }
    }

    double *a = new double[n_row*n_col];    // dist matrix -> award matrix
    for(int i=0;i<n_row;i++) {
        int bias = i*n_col;
        for(int j=0;j<n_col;j++) {
            a[bias+j] = amax-dist[bias+j];
        }
    }

    double* d1=new double[n_row+5];    // vertex labels, left part
    double* d2=new double[n_col+5];    // vertex labels, right part
    for(int i=0;i<n_row;i++)
    {
        int bias = i*n_col;
        d1[i]=0;
        for(int j=0;j<n_col;j++)
            d1[i]=max(d1[i],a[bias+j]);
    }
    for(int i=0;i<n_col;i++)
        d2[i]=0;

    int* f1=new int[n_col+5];    // minarg x in d1[x]+d2[i]-a[x*n_col+i]
    int* f2=new int[n_col+5];    // position of x in queue
    for(int i=0;i<n_col+5;i++)f1[i]=-1;
    for(int i=0;i<n_col+5;i++)f2[i]=-1;

    int N=max(n_row,n_col)+5;
    int (*lb)[3]=new int[N+5][3];    // vertex queue in km algorithm
    int l,r;

    int* bo1=new int[n_row+5];    // if been accessed, left part
    int* bo2=new int[n_col+5];    // if been accessed, right part
    for(int i=0;i<n_row+5;i++)bo1[i]=-1;
    for(int i=0;i<n_col+5;i++)bo2[i]=-1;

    int* p1=new int[n_row+5];    // match object, left part
    int* p2=new int[n_col+5];    // match object, right part
    for(int i=0;i<n_row+5;i++)p1[i]=-1;
    for(int i=0;i<n_col+5;i++)p2[i]=-1;


    for(int ii=0;ii<n_row;ii++)    // km algorithm, find match object for left part vertex
    {
        for(int i=0;i<n_col;i++)f1[i]=f2[i]=-1;
        l=0; r=1; lb[r][1]=ii; lb[r][2]=0;
        while(r>0)    // loop until find the match object of ii
        {
            int flag=-1;
            while(l<r)    // queue iterate
            {
                int x;
                l++; x=lb[l][1]; bo1[x]=ii;
                for(int i=0;i<n_col;i++)if(bo2[i]<ii)    // search the right part
                {
                    if(abs(d1[x]+d2[i]-a[x*n_col+i])<tiny)
                    {
                        bo2[i]=ii;
                        if(p2[i]==-1){ p2[i]=x; flag=l; break; }
                        r++; lb[r][0]=i; lb[r][1]=p2[i]; lb[r][2]=l;
                    }else
                    if((f1[i]==-1)||(d1[x]+d2[i]-a[x*n_col+i]<d1[f1[i]]+d2[i]-a[f1[i]*n_col+i])){ f1[i]=x; f2[i]=l; }
                }
                if(flag!=-1)break;
            }
            if(flag!=-1){ workback(flag,lb,p2); break; }


            double dd=2147483647;
            int dd1, dd2;
            for(int i=0;i<n_col;i++)if((bo2[i]<ii)&&(f1[i]>=0))    // find minist decrease
            if(d1[f1[i]]+d2[i]-a[f1[i]*n_col+i]<dd){
                dd = d1[f1[i]]+d2[i]-a[f1[i]*n_col+i];
                dd1 = f1[i];
                dd2 = i;
            }
            for(int i=0;i<n_row;i++)if(bo1[i]==ii)d1[i]-=dd;    // update vertex label, left part
            for(int i=0;i<n_col;i++)if(bo2[i]==ii)d2[i]+=dd;    // update vertex label, right part
            for(int i=0;i<n_col;i++) {
                if((bo2[i]<ii)&&(f1[i]>=0)&&(abs(d1[f1[i]]+d2[i]-a[f1[i]*n_col+i])<tiny))
                {
                    bo2[i]=ii;
                    if(p2[i]==-1){
                        p2[i]=f1[i]; flag=f2[i]; break;
                    }
                    r++; lb[r][0]=i; lb[r][1]=p2[i]; lb[r][2]=f2[i];
                }
            }
            if(flag!=-1){ workback(flag,lb,p2); break; }    // deal the augmenting path
        }
    }

    double ans=0,ans1=0,ans2=0;
    for(int i=0;i<n_row;i++)ans+=d1[i];
    for(int i=0;i<n_col;i++)ans+=d2[i];
    for(int i=0;i<n_col;i++)ans1+=a[p2[i]*n_col+i];
    for(int i=0;i<n_col;i++)p1[p2[i]]=i;
    for(int i=0;i<n_row;i++)ans2+=dist[i*n_col+p1[i]];

    std::vector<int> ret(n_row, 0);
    for(int i=0;i<n_row;i++)ret[i] = p1[i];

    delete[] a;
    delete[] d1;
    delete[] d2;
    delete[] f1;
    delete[] f2;
    delete[] bo1;
    delete[] bo2;
    delete[] p1;
    delete[] p2;
    delete[] lb;

    return ret;
}

//KM, python interface
std::vector<int> solveKM(
    const std::vector<std::vector<double>> &_dist) {

    int n_row = _dist.size();
    int n_col = _dist[0].size();

    double *dist = new double[n_row*n_col];
    for(int i=0;i<n_row;i++)
    for(int j=0;j<n_col;j++)dist[i*n_col+j]=_dist[i][j];
    std::vector<int> ret = solveKMArray(dist, n_row, n_col);

    delete[] dist;
    return ret;
}

//bi-graph match, change partly
std::vector<int> solveBiMatchChange(
const double dist[],
const int N,
const bool change[],
const int grid_asses[],
const std::string &type="km") {
    int* changeList = new int[N];
    int N2 = 0;
    for(int i=0;i<N;i++)if(change[i]){
        changeList[N2] = i;
        N2 += 1;
    }
    if(N2<N){    // only solve bi-graph match for a part of vertex
        double * new_dist = new double[N2*N2];
        #pragma omp parallel for num_threads(THREADS_NUM)
        for(int i=0;i<N2;i++){
            int gid = changeList[i];
            int bias = i*N2;
            int bias0 = gid*N;
            for(int j=0;j<N2;j++){
                int id = grid_asses[changeList[j]];
                new_dist[bias+j] = dist[bias0+id];
            }
        }

        std::vector<int> ret0(N2, 0);
        if(type=="km")
            ret0 = solveKMArray(new_dist, N2, N2);
        else
            ret0 = solveLapArray(new_dist, N2, std::max(50, int(0.15*N)));

        std::vector<int> ret(N, 0);
        for(int gid=0;gid<N;gid++)ret[gid] = grid_asses[gid];
        for(int i=0;i<N2;i++){
            int gid = changeList[i];
            int j = ret0[i];
            ret[gid] = grid_asses[changeList[j]];
        }

        delete[] new_dist;
        delete[] changeList;
        return ret;
    }else {    // solve all vertex
        std::vector<int> ret(N, 0);
        if(type=="km")
            ret = solveKMArray(dist, N, N);
        else
            ret = solveLapArray(dist, N, std::max(50, int(0.15*N)));
        delete[] changeList;
        return ret;
    }
}

// optimize the layout by iterate of bi-graph match
std::vector<double> optimizeBA(
//const std::vector<int> &_ori_grid_asses,
const std::vector<std::vector<double>> &_ori_embedded,
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels,
const std::vector<bool> &_change,
const std::string &type,
double alpha, double beta,
bool alter=false, const std::vector<double> alter_best=std::vector<double>(3,1),
int maxit=10) {

    double start = clock();

//    struct timespec tot_time1 = {0, 0};
//    struct timespec tot_time2 = {0, 0};
//    clock_gettime(CLOCK_REALTIME, &tot_time1);

    // ----------------------------------preprocess step start----------------------------------------
    printf("preprocess step start\n");

    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = max(maxLabel, _cluster_labels[i]+1);
    int *grid_asses = new int[N];
//    int *ori_grid_asses = new int[N];
    double (*ori_embedded)[2] = new double[N][2];
    int *cluster_labels = new int[num];
    bool *change = new bool[N];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
//    for(int i=0;i<N;i++)ori_grid_asses[i] = _ori_grid_asses[i];
    for(int i=0;i<N;i++) {
        ori_embedded[i][0] = _ori_embedded[i][0];
        ori_embedded[i][1] = _ori_embedded[i][1];
    }
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];
    for(int i=0;i<N;i++)change[i] = _change[i];

    double *Similar_cost_matrix = new double[N*N];

    getOriginCostMatrixArrayToArray(ori_embedded, cluster_labels, Similar_cost_matrix, N, num, square_len, maxLabel);

    if((type=="Triple")&&(global_cnt>=0)){    // clear the memory of measure by triples
        delete[] global_triples;
        delete[] global_triples_head;
        delete[] global_triples_list;
        global_cnt = -1;
        global_N = -1;
    }

    double *old_cost_matrix = new double[N*N];
    double *Compact_cost_matrix = new double[N*N];
    double *old_Convex_cost_matrix = new double[N*N];
    double *Convex_cost_matrix = new double[N*N];
    double *Cn_cost_matrix = new double[N*N];
    double *new_cost_matrix = new double[N*N];

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int i1=0;i1<N;i1++){
        int bias = i1*N;
        for(int i2=0;i2<N;i2++){
            int i = bias+i2;
            old_cost_matrix[i] = 0;
            Compact_cost_matrix[i] = old_Convex_cost_matrix[i] = Convex_cost_matrix[i] = Cn_cost_matrix[i] = new_cost_matrix[i] = 0;
        }
    }

    int *ans = new int[N];
    for(int i=0;i<N;i++)ans[i] = grid_asses[i];

    double best = 2147483647;    // cost(loss) of the best ans
    double c_best = 0;    // connectivity cost(constraint) of the best ans

    std::vector<double> best_cost(4, 2147483647);    // full cost(loss) of the best ans
    double last_c_cost;    // connectivity cost(constraint) of the last ans
    std::vector<double> last_cost(4, 2147483647);    // full cost(loss) of the last ans
    std::vector<double> last_cost2(4, 2147483647);    // full cost(loss) of the last ans

    if((alter)||(beta>0))
        getCompactCostMatrixArrayToArray(grid_asses, cluster_labels, Compact_cost_matrix, N, num, square_len, maxLabel);

    int *checked = new int[N];

    double cv_time = 0;
//    struct timespec cv_start_time = {0, 0};
//    clock_gettime(CLOCK_REALTIME, &cv_start_time);

    if(type=="CutRatio")
        last_cost = checkCostForE(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="Triple")
        last_cost = checkCostForT(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="2020")
        last_cost = checkCostFor2020(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="Global")
        last_cost = checkCostForGlobal(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta);

    if(!alter)best = last_cost[0];  //if fix alpha and beta, look for the ans with minist cost
    else best = 0;
    best_cost = last_cost;
    last_cost2 = last_cost;
    // printf("cost %.6lf %.6lf %.6lf\n",last_cost[1], last_cost[2], last_cost[3]);

    c_best = N*checkConnectForAll(grid_asses, cluster_labels, checked, N, num, square_len, maxLabel, 4);
    if(type=="Global")c_best = 0;
    last_c_cost = c_best;

    double pre_time = (clock()-start)/CLOCKS_PER_SEC;

    // ----------------------------------preprocess step done----------------------------------------
    printf("preprocess step done\n");

//    clock_gettime(CLOCK_REALTIME, &tot_time2);
//    std::cout << "pre time passed is: " << (tot_time2.tv_sec - cv_start_time.tv_sec)*1000 + (tot_time2.tv_nsec - cv_start_time.tv_nsec)/1000000 << "ms" << std::endl;
//    cv_time += ((tot_time2.tv_sec - cv_start_time.tv_sec)*1000 + (tot_time2.tv_nsec - cv_start_time.tv_nsec)/1000000);

    // ----------------------------------iterater step start----------------------------------------
    printf("iterater step start\n");

    // printf("pre time %.2lf\n", pre_time);
    double a = 1;
    int downMax = 2;
    int downCnt = 0;

    double km_time = 0;
    double cm_time = 0;
    double m_time = 0;
    double a_time = 0;

    bool use_knn = false;
    int k_value = 20;

    if((type!="CutRatio")&&(type!="Triple")&&(type!="2020")&&(type!="Global"))maxit = 0;

    // information of triples measure to save and load
    int *save_innerDict = new int[N*maxLabel];
    int *save_outerDict = new int[N*maxLabel*2];
    int *old_grid_asses = new int[N];
    double *T_pair = new double[2];

    int change_num = N;
    int ori_maxit = maxit;
    double avg_alpha=0, avg_beta=0;
    double ori_alpha=alpha, ori_beta=beta;

    int tot_it = 0;

    for(int it=0;it<maxit;it++) {    //iterate optimization
        tot_it = it+1;

        double start, tmp, start0;
        start = clock();
        start0 = clock();

//        struct timespec time1 = {0, 0};
//        clock_gettime(CLOCK_REALTIME, &time1);

        if(type!="Global"){    // only adjust border
            int rand_tmp = rand()%2;
            for(int x1=0;x1<square_len;x1++) {    // grids[x1][y1]
                int bias = x1*square_len;
                for(int y1=0;y1<square_len;y1++) {
                    int gid = bias+y1;
                    change[gid] = _change[gid];
                    if(change[gid]) {
                        change[gid] = false;
                        int lb = -1;
                        if(grid_asses[gid]<num)lb = cluster_labels[grid_asses[gid]];
                        for(int xx=-1;xx<=1;xx++) {    // grids[x1+xx][y1+yy]
                            int x2 = x1+xx;
                            if((x2<0)||(x2>=square_len))continue;
                            int bias2 = x2*square_len;
                            for(int yy=-1;yy<=1;yy++) {
                                int y2 = y1+yy;
                                if((y2<0)||(y2>=square_len))continue;
                                int gid2 = bias2+y2;
                                int lb2 = -1;
                                if(grid_asses[gid2]<num)lb2 = cluster_labels[grid_asses[gid2]];
                                if(lb2!=lb) {    // different cluters, means grids[x1][y1] is in border
                                    change[gid] = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        // get convexity cost matrix
        if(type=="CutRatio")
            getCostMatrixForEArrayToArray(grid_asses, cluster_labels, Convex_cost_matrix, N, num, square_len, maxLabel);
        else if(type=="Triple") {
//            if(it<=1) {
            // printf("change num %d\n", change_num);
            if((it<=0)||(change_num>=N/30)) {    // too many changed grids
                getCostMatrixForTArrayToArray(grid_asses, cluster_labels, Convex_cost_matrix, N, num, square_len, maxLabel,
                true, false, save_innerDict, save_outerDict);
            }else {    // runtime accelerate, only re-calculate changed grids from old layout
                getCostMatrixForTArrayToArray(grid_asses, cluster_labels, Convex_cost_matrix, N, num, square_len, maxLabel,
                true, true, save_innerDict, save_outerDict, old_grid_asses, T_pair);
            }
        }
        else if(type=="2020")
            getCostMatrixFor2020ArrayToArray(grid_asses, cluster_labels, Convex_cost_matrix, N, num, square_len, maxLabel);

//        struct timespec time0 = {0, 0};
//        clock_gettime(CLOCK_REALTIME, &time0);
//        std::cout << "cost time passed is: " << (time0.tv_sec - time1.tv_sec)*1000 + (time0.tv_nsec - time1.tv_nsec)/1000000 << "ms" << std::endl;
//        cv_time += ((time0.tv_sec - time1.tv_sec)*1000 + (time0.tv_nsec - time1.tv_nsec)/1000000);

        tmp = (clock()-start)/CLOCKS_PER_SEC;
        cm_time += tmp;
        // printf("convex maxtrix time %.2lf\n", tmp);

        getConnectCostMatrixArrayToArray(grid_asses, cluster_labels, Cn_cost_matrix, N, num, square_len, maxLabel);

        if(alter&&(type=="Global")){    // auto adjust alpha and beta
            double dec_Similar = std::max(0.0, last_cost[1]-alter_best[0])/alter_best[2]+0.0001;
            double dec_Compact = std::max(0.0, last_cost[2]-alter_best[1])/alter_best[3]+0.0001;
            alpha = 0;
            beta = dec_Compact/(dec_Similar+dec_Compact);
//            printf("new alpha: %.2lf %.2lf\n", alpha, beta);
            if(it==0)beta = ori_beta;
        }

        // if((!alter)&&(type=="Global")&&(it>=maxit-1)&&(beta==1)) {
        if((!alter)&&(type=="Global")&&(beta==1)) {
            beta = 0.99;
        }

        printf("it: %d new alpha: %.2lf %.2lf\n", it, alpha, beta);

        // calculate full cost matrix
        #pragma omp parallel for num_threads(THREADS_NUM)
        for(int i1=0;i1<N;i1++){
            int bias = i1*N;
            for(int i2=0;i2<N;i2++){
                int i = bias+i2;
                //new_cost_matrix[i] = (1-beta-alpha)*Similar_cost_matrix[i]+beta*Compact_cost_matrix[i]+alpha*(old_Convex_cost_matrix[i]+Convex_cost_matrix[i])+Cn_cost_matrix[i]*N*it/maxit;
                a = 1/(it+1.0);
                new_cost_matrix[i] = (1-beta-alpha)*Similar_cost_matrix[i]+beta*Compact_cost_matrix[i];
                old_cost_matrix[i] = new_cost_matrix[i];
                if((type!="Global")||alter) {
//                if(type!="Global") {
                    old_Convex_cost_matrix[i] = old_Convex_cost_matrix[i]*(1-a) + alpha*Convex_cost_matrix[i]*a;
                    new_cost_matrix[i] += old_Convex_cost_matrix[i];

//                    double tmp = old_Convex_cost_matrix[i];
//                    old_Convex_cost_matrix[i] = alpha*Convex_cost_matrix[i];
//                    new_cost_matrix[i] += (tmp+old_Convex_cost_matrix[i])/2;

                    new_cost_matrix[i] += Cn_cost_matrix[i]*N*it/maxit;
                }
            }
        }

        tmp = (clock()-start)/CLOCKS_PER_SEC;
        m_time += tmp;
        // printf("cost maxtrix time %.2lf\n", tmp);

        start = clock();

        std::string b_type = "km";
        // if((type=="Global")&&(alter))b_type = "lap";
        // std::cout << "type: " << b_type << std::endl;

        std::vector<int> new_asses = solveBiMatchChange(new_cost_matrix, N, change, grid_asses, b_type);   // bi-graph match

        tmp = (clock()-start)/CLOCKS_PER_SEC;
        km_time += tmp;
        // printf("km time %.2lf\n", tmp);

        start = clock();

        for(int i=0;i<N;i++)old_grid_asses[i] = grid_asses[i];
        for(int i=0;i<N;i++)grid_asses[i] = new_asses[i];

        if(((alter)||(beta>0))&&((!alter)||(it<maxit-1))) {
        // if((alter)||(beta>0)) {
            printf("re calculate comp\n");
            getCompactCostMatrixArrayToArray(grid_asses, cluster_labels, Compact_cost_matrix, N, num, square_len, maxLabel);    // re calculate
        }

        double cost = 2147483647;    // calculate the cost after bi-graph match
        std::vector<double> new_cost(4, 2147483647);
        if(type=="CutRatio")
            new_cost = checkCostForE(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
        else if(type=="Triple") {
//            if(it<=1) {
            if((it<=0)||(change_num>=N/30)) {
                new_cost = checkCostForT(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta,
                true, false, nullptr, T_pair);
            }else {    // runtime accelerate, only re-calculate changed grids from old layout
                new_cost = checkCostForT(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta,
                true, true, old_grid_asses, T_pair);
            }
        }
        else if(type=="2020")
            new_cost = checkCostFor2020(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
        else if(type=="Global")
            new_cost = checkCostForGlobal(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta);

        if(!alter)cost = new_cost[0];
        else cost = 0;

        double c_cost = N*checkConnectForAll(grid_asses, cluster_labels, checked, N, num, square_len, maxLabel, 4);
        if(type=="Global")c_cost = 0;

        if(alter&&(type=="Global")){
            printf("cost1 %.2lf %.2lf\n", new_cost[1], last_cost[1]);
            printf("cost2 %.2lf %.2lf\n", new_cost[2], last_cost[2]);
            if((std::abs(new_cost[1]-last_cost2[1])+std::abs(new_cost[2]-last_cost2[2]))/N<0.00015) {
                printf("converge 1\n");
//                break;
                maxit = it+1;
            }
            if((std::abs(new_cost[1]-last_cost[1])+std::abs(new_cost[2]-last_cost[2]))/N<0.0015) {
                printf("converge 2\n");
//                break;
                maxit = it+1;
            }
        }

        last_cost2 = last_cost;
        last_cost = new_cost;
        last_c_cost = c_cost;

        if(((alter)&&(it==0))||(cost+c_cost<=best+c_best)) {    // update ans
//        if(true) {
            best = cost;
            c_best = c_cost;
            best_cost = new_cost;
            for(int i=0;i<N;i++)ans[i] = grid_asses[i];
            downCnt = 0;
        }else downCnt += 1;

        change_num = 0;
        for(int i=0;i<N;i++)if(cluster_labels[grid_asses[i]]!=cluster_labels[old_grid_asses[i]])change_num += 1;

        printf("cost %.2lf %.2lf %d\n", cost, c_cost, downCnt);
        printf("cost prox %.2lf comp %.2lf\n", new_cost[1], new_cost[2]);

        tmp = (clock()-start)/CLOCKS_PER_SEC;
        a_time += tmp;
        // printf("update time %.2lf\n", tmp);

//       struct timespec time2 = {0, 0};
//       clock_gettime(CLOCK_REALTIME, &time2);
//       std::cout << "it time passed is: " << (time2.tv_sec - time1.tv_sec)*1000 + (time2.tv_nsec - time1.tv_nsec)/1000000 << "ms" << std::endl;

        if((type=="Triple")||(type=="2020")) {   // case that need to ensure the connectivity in this function
//            if(((!alter)&&(c_best>0))||((alter)&&(last_c_cost>0))){
            if(c_best>0){
                if(it >= maxit-1)maxit += 1;
                if(maxit>2*ori_maxit)break;
                continue;
            }
        }
        if(downCnt>=downMax)break;
    }

    // printf("tot convex matrix time %.2lf\n", cm_time);
    // printf("tot cost matrix time %.2lf\n", m_time);
    // printf("tot km time %.2lf\n", km_time);
    // printf("tot update time %.2lf\n", a_time);

    // ----------------------------------iterater step done----------------------------------------

//    clock_gettime(CLOCK_REALTIME, &tot_time2);
//    std::cout << "it done time passed is: " << (tot_time2.tv_sec - tot_time1.tv_sec)*1000 + (tot_time2.tv_nsec - tot_time1.tv_nsec)/1000000 << "ms" << std::endl;

    std::vector<double> ret(N+6, 0);

    for(int i=0;i<N;i++)ret[i] = ans[i];
    for(int i=0;i<3;i++)ret[N+i] = best_cost[i+1];
    ret[N+3] = tot_it;
    ret[N+4] = beta;
    ret[N+5] = cv_time/1000;
    printf("cv time %.3lf\n", cv_time/1000);

//    clock_gettime(CLOCK_REALTIME, &tot_time2);
//    std::cout << "ret time passed is: " << (tot_time2.tv_sec - tot_time1.tv_sec)*1000 + (tot_time2.tv_nsec - tot_time1.tv_nsec)/1000000 << "ms" << std::endl;

    if((type=="Triple")&&(global_cnt>=0)){
        delete[] global_triples;
        delete[] global_triples_head;
        delete[] global_triples_list;
        global_cnt = -1;
        global_N = -1;
    }

//    clock_gettime(CLOCK_REALTIME, &tot_time2);
//    std::cout << "del time passed is: " << (tot_time2.tv_sec - tot_time1.tv_sec)*1000 + (tot_time2.tv_nsec - tot_time1.tv_nsec)/1000000 << "ms" << std::endl;

    delete[] save_innerDict;
    delete[] save_outerDict;
    delete[] old_grid_asses;
    delete[] T_pair;

//    clock_gettime(CLOCK_REALTIME, &tot_time2);
//    std::cout << "del2 time passed is: " << (tot_time2.tv_sec - tot_time1.tv_sec)*1000 + (tot_time2.tv_nsec - tot_time1.tv_nsec)/1000000 << "ms" << std::endl;

    delete[] checked;
    delete[] grid_asses;
//    delete[] ori_grid_asses;
    delete[] ori_embedded;
    delete[] cluster_labels;
    delete[] change;

//    clock_gettime(CLOCK_REALTIME, &tot_time2);
//    std::cout << "del3 time passed is: " << (tot_time2.tv_sec - tot_time1.tv_sec)*1000 + (tot_time2.tv_nsec - tot_time1.tv_nsec)/1000000 << "ms" << std::endl;

    delete[] Similar_cost_matrix;
    delete[] Compact_cost_matrix;
    delete[] Convex_cost_matrix;
    delete[] old_Convex_cost_matrix;
    delete[] Cn_cost_matrix;
    delete[] new_cost_matrix;
    delete[] ans;

//    clock_gettime(CLOCK_REALTIME, &tot_time2);
//    std::cout << "tot time passed is: " << (tot_time2.tv_sec - tot_time1.tv_sec)*1000 + (tot_time2.tv_nsec - tot_time1.tv_nsec)/1000000 << "ms" << std::endl;

    return ret;
}

// optimize by search bars to swap
std::vector<double> optimizeSwap(
//const std::vector<int> &_ori_grid_asses,
const std::vector<std::vector<double>> &_ori_embedded,
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels,
const std::vector<bool> &_change,
const std::string &type,
double alpha, double beta,
int maxit=10, int choose_k=1, int seed=10, bool innerBiMatch=true, int swap_cnt=214748347) {

    // ----------------------------------preprocess step start----------------------------------------

    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = max(maxLabel, _cluster_labels[i]+1);
    int *grid_asses = new int[N];
//    int *ori_grid_asses = new int[N];
    double (*ori_embedded)[2] = new double[N][2];
    int *cluster_labels = new int[num];
    bool *change = new bool[N];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
//    for(int i=0;i<N;i++)ori_grid_asses[i] = _ori_grid_asses[i];
    for(int i=0;i<N;i++) {
        ori_embedded[i][0] = _ori_embedded[i][0];
        ori_embedded[i][1] = _ori_embedded[i][1];
    }
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];
    for(int i=0;i<N;i++)change[i] = _change[i];

    double *Similar_cost_matrix = new double[N*N];
    getOriginCostMatrixArrayToArray(ori_embedded, cluster_labels, Similar_cost_matrix, N, num, square_len, maxLabel);

    double *Compact_cost_matrix = new double[N*N];
    getCompactCostMatrixArrayToArray(grid_asses, cluster_labels, Compact_cost_matrix, N, num, square_len, maxLabel);

    int *ans = new int[N];
    for(int i=0;i<N;i++)ans[i] = grid_asses[i];
    double best = 2147483647;    // cost of ans
    double c_best = 0;    // connectivity cost(constraint) of ans

    double *old_T_pair = new double[2];    // convexity of triples
    double *old_D_pair = new double[N*N*2];
    int *old_grid_asses = new int[N];
    for(int i=0;i<N;i++)old_grid_asses[i] = grid_asses[i];

    bool *if_disconn = new bool[N];   // disconnect of now grids
    for(int i=0;i<N;i++)if_disconn[i] = false;
    int *checked = new int[N];

    if(type=="CutRatio")
        best = checkCostForE(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
    else if(type=="AreaRatio")
        best = checkCostForS(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
    else if(type=="PerimeterRatio")
        best = checkCostForC(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
    else if(type=="2020")
        best = checkCostFor2020(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
    else if(type=="TB")
        best = checkCostForTB(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
    else if(type=="Triple")
        best = checkCostForT(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, false, old_grid_asses, old_T_pair)[0];
    else if(type=="Global")
        best = checkCostForGlobal(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
    else if(type=="T2")
        best = checkCostForT2(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, false, old_grid_asses, old_T_pair, old_D_pair)[0];

    c_best = N*checkConnectForAll(grid_asses, cluster_labels, checked, N, num, square_len, maxLabel, 8, if_disconn);

    // ----------------------------------preprocess step done----------------------------------------

    // ----------------------------------swap step start----------------------------------------

    int downMax = 3;
    int downCnt = 0;
    srand(seed);

    int check_turns1 = 0;
    int check_turns2 = 0;
    int once_flag = 1;
    double s_time = 0;
    double sp_time = 0;
    double a_time = 0;

    bool use_knn = false;
    int k_value = 50;

    int max_num = 2;    // max bar length to search
    for(int it=0;it<maxit;it++){    //枚举轮次
        if(type=="CutRatio"||type=="AreaRatio"||type=="PerimeterRatio"||type=="2020"||type=="Triple"||type=="T2"||type=="TB"){

            int update_flag = 0;    // if layout been updated
            int improve_flag = 0;    // if find a better layout
            // if(type=="AreaRatio")max_num = 1;
            // if(type=="Triple")max_num = 1;
            // if(type=="2020")max_num = 1;
            // if(type=="PerimeterRatio")max_num = 6;

            double (*label_pairs)[2] = new double[maxLabel][2];    // convexity of every cluster
            int *worst_gid = new int[2*N*max_num];    // bar2
            int worst_cnt = 0;    // number of bar2
            int *order1 = new int[2*N];    // enumerate order of bar1
            int *order2 = new int[2*N];    // enumerate order of bar2
            int *E_grid = new int[N];    // count of edges with a different clusters of every grid
            int *E_grid2 = new int[N];    // count of edges with the bar1 clusters information of every grid
            int *labels_cnt = new int[maxLabel+1];    // count of edges with every clusters, bar1
            int *labels_cnt2 = new int[maxLabel+1];    // count of edges with every clusters, bar2
            int check_cnt = 1;
            for(int i=0;i<N;i++)checked[i]=0;

            double *all_cost = new double[2*N];    // cost of swapping
            double *all_c_cost = new double[2*N];    // connectivity cost(constraint) of swapping
            double *all_b_cost = new double[2*N];    // edges with blank grids
            double *all_dec = new double[2*N];    // edges with blank grids

            double start = clock();

            for(int now_num=max_num;now_num>=1;now_num--){    // long bar first

                for(int i=0;i<N;i++){
                    order1[i] = i;
                }

                std::random_shuffle(order1, order1+N);    //shuffle the order

                for(int gid1_o=0;gid1_o<N;gid1_o++){    //enumerate bar1
                    int gid1 = order1[gid1_o];
                    int id1 = grid_asses[gid1];
                    if(id1>=num)continue;
                    int mainLabel1 = cluster_labels[id1];
                    int x1 = gid1/square_len;
                    int y1 = gid1%square_len;
                    for(int ori=0;ori<2;ori++){
                        if((ori>0)&&(now_num==1))continue;
                        int ori_gid1[100];
                        if((ori==0)&&(x1>=now_num-1)){    //V bar
                            for(int i=0;i<now_num;i++)ori_gid1[i] = gid1 - i*square_len;
                        }else
                        if((ori==1)&&(y1>=now_num-1)){    //H bar
                            for(int i=0;i<now_num;i++)ori_gid1[i] = gid1 - i;
                        }else continue;

                        for(int i=0;i<maxLabel+1;i++)labels_cnt[i] = 0;
                        for(int i=0;i<maxLabel+1;i++)labels_cnt2[i] = 0;

                        int flag=0;
                        int dc_flag=0;

                        for(int i=0;i<now_num;i++){
                            int gid = ori_gid1[i];
                            int id = grid_asses[gid];
                            if((id>=num)||(cluster_labels[id]!=cluster_labels[id1])){
                                flag = 1;
                                break;
                            }
                            if(!change[gid]){
                                flag = 1;
                                break;
                            }
                            checkEdgeSingleForLabel(labels_cnt, gid, grid_asses, cluster_labels, N, num, square_len, maxLabel);
                        }
                        if(flag>0)continue;    //illegal bar

                        for(int i=0;i<now_num;i++){
                            int gid = ori_gid1[i];
                            dc_flag += if_disconn[gid];
                        }

                        int dcnt = 0, mainLabel2 = maxLabel;
                        for(int i=0;i<maxLabel+1;i++){
                            if(i==cluster_labels[id1])continue;
                            dcnt += labels_cnt[i];
                            if(labels_cnt[i]>labels_cnt[mainLabel2])mainLabel2 = i;
                        }
                        if(dcnt<=now_num)continue;    //not enough edges with a different clusters

                        // if((type=="Triple")&&(dcnt<=now_num+1))continue;

                        check_turns1 += 1;

                        checkEdgeArrayForSingleLabel(E_grid, mainLabel1, grid_asses, cluster_labels, N, num, square_len, maxLabel);
                        checkEdgeArray(E_grid2, grid_asses, cluster_labels, N, num, square_len, maxLabel);
                        worst_cnt = 0;

                        for(int gid2=0;gid2<N;gid2++){    // enumerate bar2
                            int id2 = grid_asses[gid2];
                            int lb2 = maxLabel;
                            if(id2<num)lb2 = cluster_labels[id2];
                            if((id2<num)&&(labels_cnt[lb2]==0))continue;   // cluster that have no edge with bar1
                            int x2 = gid2/square_len;
                            int y2 = gid2%square_len;
                            for(int ori=0;ori<2;ori++){
                                if((ori>0)&&(now_num==1))continue;
                                int now_gid2[100];
                                if((ori==0)&&(x2>=now_num-1)){    //V bar
                                    for(int i=0;i<now_num;i++)now_gid2[i] = gid2 - i*square_len;
                                }else
                                if((ori==1)&&(y2>=now_num-1)){    //H bar
                                    for(int i=0;i<now_num;i++)now_gid2[i] = gid2 - i;
                                }else continue;

                                int flag=0;
                                int dcnt=0;
                                int dcnt2=0;
                                for(int i=0;i<now_num;i++){
                                    int gid = now_gid2[i];
                                    int id = grid_asses[gid];
                                    int lb = maxLabel;
                                    if(id<num)lb = cluster_labels[id];
                                    if(lb!=lb2){
                                        flag = 1;
                                        break;
                                    }
                                    if(!change[gid]){
                                        flag = 1;
                                        break;
                                    }
                                    dcnt += E_grid[gid];
                                    dcnt2 += E_grid2[gid];
                                }

                                if(flag>0)continue;
                                if(dc_flag==0) {    // bar1 is connective with cluster
                                    if(dcnt<=max(now_num-2,0))continue;    // not enough edges with the bar1 cluster
                                    if(dcnt2<=now_num)continue;    // not enough edges with a different cluster
                                }else {
                                    if(dcnt2<=max(now_num-1,0))continue;    // not enough edges with a different cluster
                                }

                                // if((type=="Triple")&&(x1!=x2)&&(y1!=y2))continue;

                                for(int i=0;i<now_num;i++)worst_gid[worst_cnt*now_num+i] = now_gid2[i];    // be candidate
                                worst_cnt += 1;
                            }
                        }

                        if(worst_cnt==0)continue;

                        for(int i=0;i<worst_cnt;i++){
                            order2[i] = i;
                        }
                        // std::random_shuffle(order2, order2+worst_cnt);

                        double now_best = 0;    // cost of swapping with best bar2
                        double now_c_best = 0;    // connectivity cost(constraint) of swapping with best bar2
                        double now_b_best = 0;    // edges with blank grids
                        double ori_cost = 0;    // cost
                        double ori_c_cost = 0;    // connectivity cost(constraint)

                        if(type=="CutRatio")
                            now_best =checkCostForE(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
                        else if(type=="AreaRatio")
                            now_best =checkCostForS(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, false, label_pairs)[0];
                        else if(type=="PerimeterRatio")
                            now_best =checkCostForC(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, false, label_pairs)[0];
                        else if(type=="2020")
                            now_best =checkCostFor2020(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, false, label_pairs)[0];
                        else if(type=="TB")
                            now_best =checkCostForTB(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, false, label_pairs)[0];
                        else if(type=="Triple") {
                            now_best =checkCostForT(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, true, old_grid_asses, old_T_pair)[0];
                            for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = grid_asses[tmp_gid];
                        }
                        else if(type=="T2") {
                            now_best =checkCostForT2(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, true, old_grid_asses, old_T_pair, old_D_pair)[0];
                            for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = grid_asses[tmp_gid];
                        }

                        now_c_best = N*checkConnectForAll(grid_asses, cluster_labels, checked, N, num, square_len, maxLabel);

                        ori_cost = now_best;
                        ori_c_cost = now_c_best;

                        int best_gid = -1;    // best bar2

//                        printf("worst cnt %d\n", worst_cnt);

                        double tmp_start = clock();

                        int tmp_th = std::min(THREADS_NUM, worst_cnt);
                        #pragma omp parallel for num_threads(8)
                        for(int jj=0;jj<worst_cnt;jj++){    // enumerate candidate bar2
                            check_turns2 += 1;
                            int j = order2[jj];

                            all_dec[j] = -1;

                            int flag = 0;
                            for(int k1=0;k1<now_num;k1++)
                            for(int k2=0;k2<now_num;k2++){
                                if(ori_gid1[k1]==worst_gid[j*now_num+k2]){
                                    flag += 1;
                                }
                            }
                            if(flag>0)continue;    // have share grid

                            int id1 = grid_asses[ori_gid1[0]];
                            int id2 = grid_asses[worst_gid[j*now_num]];
                            int lb1 = cluster_labels[id1];
                            int lb2 = maxLabel;
                            if(id2<num)lb2 = cluster_labels[id2];
                            if(lb1==lb2)continue;    // same cluster

                            int *tmp_grid_asses = new int[N];
                            int *tmp_checked = new int[N];
                            memcpy(tmp_grid_asses, grid_asses, N*sizeof(int));

                            for(int k=0;k<now_num;k++){    //swap
                                std::swap(tmp_grid_asses[ori_gid1[k]], tmp_grid_asses[worst_gid[j*now_num+k]]);
                            }

                            double cost = 0;    // new cost
                            double c_cost = 0;    // new connectivity cost
                            double b_cost = 0;

                            if(type=="CutRatio")
                                cost = checkCostForE(Similar_cost_matrix, Compact_cost_matrix, tmp_grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
                            else if(type=="AreaRatio")
                                cost = checkCostForS(Similar_cost_matrix, Compact_cost_matrix, tmp_grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, false, true, label_pairs, lb1, lb2)[0];
                            else if(type=="PerimeterRatio")
                                cost = checkCostForC(Similar_cost_matrix, Compact_cost_matrix, tmp_grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, false, true, label_pairs, lb1, lb2)[0];
                            else if(type=="2020")
                                cost = checkCostFor2020(Similar_cost_matrix, Compact_cost_matrix, tmp_grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, false, true, label_pairs, lb1, lb2)[0];
                            else if(type=="TB")
                                cost = checkCostForTB(Similar_cost_matrix, Compact_cost_matrix, tmp_grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, false, true, label_pairs, lb1, lb2)[0];
                            else if(type=="Triple")
                                cost = checkCostForT(Similar_cost_matrix, Compact_cost_matrix, tmp_grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, false, true, old_grid_asses, old_T_pair)[0];
                            else if(type=="T2")
                                cost = checkCostForT2(Similar_cost_matrix, Compact_cost_matrix, tmp_grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, false, true, old_grid_asses, old_T_pair, old_D_pair)[0];

                            c_cost = N*checkConnectForAll(tmp_grid_asses, cluster_labels, tmp_checked, N, num, square_len, maxLabel);

                            for(int k=0;k<now_num;k++){
                                b_cost -= 0.000001*checkBlankForGrid(ori_gid1[k], tmp_grid_asses, cluster_labels, N, num, square_len, maxLabel);
                            }
                            for(int k=0;k<now_num;k++){
                                b_cost -= 0.000001*checkBlankForGrid(worst_gid[j*now_num+k], tmp_grid_asses, cluster_labels, N, num, square_len, maxLabel);
                            }

                            for(int k=0;k<now_num;k++){    //swap back
                                std::swap(tmp_grid_asses[ori_gid1[k]], tmp_grid_asses[worst_gid[j*now_num+k]]);
                            }

                            for(int k=0;k<now_num;k++){
                                b_cost += 0.000001*checkBlankForGrid(ori_gid1[k], tmp_grid_asses, cluster_labels, N, num, square_len, maxLabel);
                            }
                            for(int k=0;k<now_num;k++){
                                b_cost += 0.000001*checkBlankForGrid(worst_gid[j*now_num+k], tmp_grid_asses, cluster_labels, N, num, square_len, maxLabel);
                            }

                            all_cost[j] = cost;
                            all_c_cost[j] = c_cost;
                            all_b_cost[j] = b_cost;
                            all_dec[j] = ori_cost + ori_c_cost - cost - c_cost - b_cost;

//                            if(cost+c_cost+b_cost<now_best+now_c_best+now_b_best){
//                                now_best = cost;
//                                now_c_best = c_cost;
//                                now_b_best = b_cost;
//                                best_gid = j;
//                            }

                            // printf("done swap\n");

                            delete[] tmp_grid_asses;
                            delete[] tmp_checked;
                        }

                        for(int jj=0;jj<worst_cnt;jj++){    // enumerate candidate bar2

                            check_turns2 += 1;
                            int j = order2[jj];

                            double cost = all_cost[j];
                            double c_cost = all_c_cost[j];
                            double b_cost = all_b_cost[j];

                            if(ori_cost + ori_c_cost - all_dec[j]<now_best+now_c_best+now_b_best){
                                now_best = cost;
                                now_c_best = c_cost;
                                now_b_best = b_cost;
                                best_gid = j;
                            }
                        }

                        sp_time += (clock()-tmp_start)/CLOCKS_PER_SEC;

                        if(best_gid!=-1){    // choose the best bar2 to swap

//                            best_gid = soft_choose(all_dec, worst_cnt);
                            best_gid = best_k_choose(all_dec, worst_cnt, choose_k);

                            if(swap_cnt<=0) {
                                break;
                            }

                            update_flag = 1;

                            double cost = now_best;
                            double c_cost = now_c_best;

                            for(int k=0;k<now_num;k++){
                                std::swap(grid_asses[ori_gid1[k]], grid_asses[worst_gid[best_gid*now_num+k]]);
                            }

                            if(type=="CutRatio")
                                cost = checkCostForE(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
                            else if(type=="AreaRatio")
                                cost = checkCostForS(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
                            else if(type=="PerimeterRatio")
                                cost = checkCostForC(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
                            else if(type=="2020")
                                cost = checkCostFor2020(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
                            else if(type=="TB")
                                cost = checkCostForTB(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta)[0];
                            else if(type=="Triple") {
                                cost = checkCostForT(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, true, old_grid_asses, old_T_pair)[0];
                                for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = grid_asses[tmp_gid];
                            }
                            else if(type=="T2") {
                                cost = checkCostForT2(Similar_cost_matrix, Compact_cost_matrix, grid_asses, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, true, old_grid_asses, old_T_pair, old_D_pair)[0];
                                for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = grid_asses[tmp_gid];
                            }

                            c_cost = N*checkConnectForAll(grid_asses, cluster_labels, checked, N, num, square_len, maxLabel, 8, if_disconn);  // update disconnect

                            if(cost+c_cost+now_b_best<best+c_best){    // find a better layout
                                // printf("swap %d %d\n", grid_asses[ori_gid1[0]], grid_asses[worst_gid[best_gid*now_num]]);

                                swap_cnt -= 1;

                                best = cost;
                                c_best = c_cost;

                                for(int k=0;k<N;k++)ans[k] = grid_asses[k];
                                improve_flag = 1;
                            }
                        }
                    }
                }
            }

            delete[] label_pairs;
            delete[] E_grid;
            delete[] E_grid2;
            delete[] worst_gid;
            delete[] order1;
            delete[] order2;
            delete[] labels_cnt;
            delete[] labels_cnt2;

            delete[] all_cost;
            delete[] all_c_cost;
            delete[] all_b_cost;
            delete[] all_dec;

            double tmp = (clock()-start)/CLOCKS_PER_SEC;
            s_time += tmp;
            // printf("swap time %.2lf\n", tmp);

            start = clock();

            max_num += 2;

            if(improve_flag==1){
                downCnt=0;
            }else downCnt += 1;

            // printf("downCnt %d %d\n", it, downCnt);
            // printf("------------------------------\n");
            // printf("------------------------------\n");
            // printf("------------------------------\n");
            // printf("------------------------------\n");

            tmp = (clock()-start)/CLOCKS_PER_SEC;
            a_time += tmp;
            // printf("addition time %.2lf\n", tmp);

            printf("downCnt %d\n", downCnt);
            if(downCnt>=downMax) {
//                if(type=="PerimeterRatio")continue;
                break;
            }
        }
    }

    // ----------------------------------swap step done----------------------------------------

    if(innerBiMatch) {
        double *cost_matrix = new double[N*N];    //bi-graph match in each clusters to ensure the minist Similar cost

        #pragma omp parallel for num_threads(THREADS_NUM)
        for(int gid2=0;gid2<N;gid2++){
            int lb2 = -1;
            int id2 = ans[gid2];
            if(id2<num)lb2 = cluster_labels[id2];

            for(int gid1=0;gid1<N;gid1++){
                int lb1 = -1;
                int id1 = ans[gid1];
                if(id1<num)lb1 = cluster_labels[id1];
                if(lb1==lb2){
                    cost_matrix[gid2*N+id1] = Similar_cost_matrix[gid2*N+id1];
                }else {
                    cost_matrix[gid2*N+id1] = N;
                }
            }
        }

        std::vector<int> new_grid_asses = solveBiMatchChange(cost_matrix, N, change, ans);

        delete[] cost_matrix;
        for(int i=0;i<N;i++)ans[i] = new_grid_asses[i];
    }
    // printf("check turns %d %d\n", check_turns1, check_turns2);
//     printf("tot swap time %.2lf\n", s_time);
//     printf("tot swap2 time %.2lf\n", sp_time);
//     printf("tot addition time %.2lf\n", a_time);

    // printf("final cost\n");
    std::vector<double> cost(4, -1);
    if(type=="CutRatio")
        cost = checkCostForE(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="AreaRatio")
        cost = checkCostForS(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="PerimeterRatio")
        cost = checkCostForC(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="2020")
        cost = checkCostFor2020(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="TB")
        cost = checkCostForTB(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="Triple") {
        cost = checkCostForT(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, true, old_grid_asses, old_T_pair);
        for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = ans[tmp_gid];
    }
    else if(type=="T2") {
        cost = checkCostForT2(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, true, old_grid_asses, old_T_pair, old_D_pair);
        for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = ans[tmp_gid];
    }
    else if(type=="Global")
        cost = checkCostForGlobal(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    // printf("cost %.6lf %.6lf %.6lf\n", cost[1], cost[2], cost[3]);

    std::vector<double> ret(N+3, 0);
    for(int i=0;i<N;i++)ret[i] = ans[i];
    for(int i=0;i<3;i++)ret[N+i] = cost[i+1];

    delete[] checked;
    delete[] if_disconn;

    delete[] old_grid_asses;
    delete[] old_T_pair;
    delete[] old_D_pair;

    delete[] grid_asses;
    delete[] ori_embedded;
    delete[] cluster_labels;
    delete[] change;
    delete[] Similar_cost_matrix;
    delete[] Compact_cost_matrix;
    delete[] ans;

    return ret;
}

// check cost
std::vector<double> checkCostForOne(
const std::vector<std::vector<double>> &_ori_embedded,
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels,
const std::string &type,
double alpha, double beta,
int min_grids=0) {
    min_grids = std::max(1, min_grids);

    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = max(maxLabel, _cluster_labels[i]+1);
    int *grid_asses = new int[N];
//    int *ori_grid_asses = new int[N];
    double (*ori_embedded)[2] = new double[N][2];
    int *cluster_labels = new int[num];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
//    for(int i=0;i<N;i++)ori_grid_asses[i] = _ori_grid_asses[i];
    for(int i=0;i<N;i++) {
        ori_embedded[i][0] = _ori_embedded[i][0];
        ori_embedded[i][1] = _ori_embedded[i][1];
    }
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];

    double *Similar_cost_matrix = new double[N*N];
    getOriginCostMatrixArrayToArray(ori_embedded, cluster_labels, Similar_cost_matrix, N, num, square_len, maxLabel);

    double *Compact_cost_matrix = new double[N*N];
    getCompactCostMatrixArrayToArray(grid_asses, cluster_labels, Compact_cost_matrix, N, num, square_len, maxLabel);

    double *old_T_pair = new double[2];    // convexity of triples
    double *old_D_pair = new double[N*N*2];
    int *old_grid_asses = new int[N];
    for(int i=0;i<N;i++)old_grid_asses[i] = grid_asses[i];

    bool *if_disconn = new bool[N];   // disconnect of now grids
    for(int i=0;i<N;i++)if_disconn[i] = false;
    int *checked = new int[N];
    double c_cost = checkConnectForAll(grid_asses, cluster_labels, checked, N, num, square_len, maxLabel, 8, if_disconn);

    std::vector<double> ret(0, 0);
    for(int lb=0;lb<maxLabel;lb++) {
//        printf("start lb %d\n", lb);
        int *ans = new int[N];
        for(int gid=0;gid<N;gid++)ans[gid] = -1;
        int cnt = 0;
        for(int gid=0;gid<N;gid++)
        if((if_disconn[gid]==false)&&((grid_asses[gid]<num)&&(cluster_labels[grid_asses[gid]]==lb))) {
            ans[gid] = cnt;
            cnt += 1;
        }
        if(cnt<min_grids) {
            delete[] ans;
            continue;
        }

        int cnt2 = cnt;
        int *ans_labels = new int[cnt];

        for(int id=0;id<cnt;id++)ans_labels[id] = 0;
        for(int gid=0;gid<N;gid++)
        if(ans[gid]==-1) {
            ans[gid] = cnt2;
            cnt2 += 1;
        }

        std::vector<double> cost(4, -1);
        if(type=="CutRatio")
            cost = checkCostForE(Similar_cost_matrix, Compact_cost_matrix, ans, ans_labels, N, cnt, square_len, 1, alpha, beta);
        else if(type=="AreaRatio")
            cost = checkCostForS(Similar_cost_matrix, Compact_cost_matrix, ans, ans_labels, N, cnt, square_len, 1, alpha, beta);
        else if(type=="PerimeterRatio")
            cost = checkCostForC(Similar_cost_matrix, Compact_cost_matrix, ans, ans_labels, N, cnt, square_len, 1, alpha, beta);
        else if(type=="2020")
            cost = checkCostFor2020(Similar_cost_matrix, Compact_cost_matrix, ans, ans_labels, N, cnt, square_len, 1, alpha, beta);
        else if(type=="TB")
            cost = checkCostForTB(Similar_cost_matrix, Compact_cost_matrix, ans, ans_labels, N, cnt, square_len, 1, alpha, beta);
        else if(type=="Triple") {
            cost = checkCostForT(Similar_cost_matrix, Compact_cost_matrix, ans, ans_labels, N, cnt, square_len, 1, alpha, beta, true, false, old_grid_asses, old_T_pair);
            for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = ans[tmp_gid];
        }
        else if(type=="T2") {
            cost = checkCostForT2(Similar_cost_matrix, Compact_cost_matrix, ans, ans_labels, N, cnt, square_len, 1, alpha, beta, true, false, old_grid_asses, old_T_pair, old_D_pair);
            for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = ans[tmp_gid];
        }
        else if(type=="B")
            cost = checkCostForB(Similar_cost_matrix, Compact_cost_matrix, ans, ans_labels, N, cnt, square_len, 1, alpha, beta);
        else if(type=="Dev")
            cost = checkCostForDev(Similar_cost_matrix, Compact_cost_matrix, ans, ans_labels, N, cnt, square_len, 1, alpha, beta);
        else if(type=="AlphaT")
            cost = checkCostForAlphaT(Similar_cost_matrix, Compact_cost_matrix, ans, ans_labels, N, cnt, square_len, 1, alpha, beta);
        else if(type=="MoS")
            cost = checkCostForMoS(Similar_cost_matrix, Compact_cost_matrix, ans, ans_labels, N, cnt, square_len, 1, alpha, beta);
        // printf("cost %.6lf %.6lf %.6lf\n", cost[1], cost[2], cost[3]);


        ret.push_back(cost[3]);

        delete[] ans;
        delete[] ans_labels;
    }

    delete[] if_disconn;
    delete[] checked;

    delete[] old_grid_asses;
    delete[] old_T_pair;
    delete[] old_D_pair;

    delete[] grid_asses;
    delete[] ori_embedded;
    delete[] cluster_labels;
    delete[] Similar_cost_matrix;
    delete[] Compact_cost_matrix;

    return ret;
}

// check cost
std::vector<double> checkCostForAllShapes(
const std::vector<std::vector<double>> &_ori_embedded,
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels,
const std::string &type,
double alpha, double beta,
int min_grids=0) {
    min_grids = std::max(1, min_grids);

    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = max(maxLabel, _cluster_labels[i]+1);
    int *grid_asses = new int[N];
//    int *ori_grid_asses = new int[N];
    double (*ori_embedded)[2] = new double[N][2];
    int *cluster_labels = new int[num];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
//    for(int i=0;i<N;i++)ori_grid_asses[i] = _ori_grid_asses[i];
    for(int i=0;i<N;i++) {
        ori_embedded[i][0] = _ori_embedded[i][0];
        ori_embedded[i][1] = _ori_embedded[i][1];
    }
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];

    double *Similar_cost_matrix = new double[N*N];
    getOriginCostMatrixArrayToArray(ori_embedded, cluster_labels, Similar_cost_matrix, N, num, square_len, maxLabel);

    double *Compact_cost_matrix = new double[N*N];
    getCompactCostMatrixArrayToArray(grid_asses, cluster_labels, Compact_cost_matrix, N, num, square_len, maxLabel);

    double *old_T_pair = new double[2];    // convexity of triples
    double *old_D_pair = new double[N*N*2];
    int *old_grid_asses = new int[N];
    for(int i=0;i<N;i++)old_grid_asses[i] = grid_asses[i];

    bool *if_disconn = new bool[N];   // disconnect of now grids
    for(int i=0;i<N;i++)if_disconn[i] = false;
    int *checked = new int[N];
    double c_cost = checkConnectForAll(grid_asses, cluster_labels, checked, N, num, square_len, maxLabel, 8, if_disconn);

    std::vector<double> ret(0, 0);
    int clusters_cnt = 0;
    int all_cnt = 0;
    int *all_ans = new int[N];
    int *all_labels = new int[N];
    for(int gid=0;gid<N;gid++)all_ans[gid] = -1;

    for(int lb=0;lb<maxLabel;lb++) {
//        printf("start lb %d\n", lb);
        int cnt = 0;
        for(int gid=0;gid<N;gid++)
        if((if_disconn[gid]==false)&&((grid_asses[gid]<num)&&(cluster_labels[grid_asses[gid]]==lb))) {
            cnt += 1;
        }
        if(cnt<min_grids) {
            continue;
        }

        for(int gid=0;gid<N;gid++)
        if((if_disconn[gid]==false)&&((grid_asses[gid]<num)&&(cluster_labels[grid_asses[gid]]==lb))) {
            all_ans[gid] = all_cnt;
            all_labels[all_cnt] = clusters_cnt;
            all_cnt += 1;
        }
        clusters_cnt += 1;
    }

    int all_cnt2 = all_cnt;

    for(int gid=0;gid<N;gid++)
    if(all_ans[gid]==-1) {
        all_ans[gid] = all_cnt2;
        all_cnt2 += 1;
    }

    if(all_cnt==0) {
        delete[] all_ans;
        delete[] all_labels;

        delete[] if_disconn;
        delete[] checked;

        delete[] old_grid_asses;
        delete[] old_T_pair;
        delete[] old_D_pair;

        delete[] grid_asses;
        delete[] ori_embedded;
        delete[] cluster_labels;
        delete[] Similar_cost_matrix;
        delete[] Compact_cost_matrix;
        return ret;
    }

    std::vector<double> cost(4, -1);

    if(type=="CutRatio")
        cost = checkCostForE(Similar_cost_matrix, Compact_cost_matrix, all_ans, all_labels, N, all_cnt, square_len, clusters_cnt, alpha, beta);
    else if(type=="AreaRatio")
        cost = checkCostForS(Similar_cost_matrix, Compact_cost_matrix, all_ans, all_labels, N, all_cnt, square_len, clusters_cnt, alpha, beta);
    else if(type=="PerimeterRatio")
        cost = checkCostForC(Similar_cost_matrix, Compact_cost_matrix, all_ans, all_labels, N, all_cnt, square_len, clusters_cnt, alpha, beta);
    else if(type=="2020")
        cost = checkCostFor2020(Similar_cost_matrix, Compact_cost_matrix, all_ans, all_labels, N, all_cnt, square_len, clusters_cnt, alpha, beta);
    else if(type=="TB")
        cost = checkCostForTB(Similar_cost_matrix, Compact_cost_matrix, all_ans, all_labels, N, all_cnt, square_len, clusters_cnt, alpha, beta);
    else if(type=="Triple") {
        cost = checkCostForT(Similar_cost_matrix, Compact_cost_matrix, all_ans, all_labels, N, all_cnt, square_len, clusters_cnt, alpha, beta, true, false, old_grid_asses, old_T_pair);
        for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = all_ans[tmp_gid];
    }
    else if(type=="T2") {
        cost = checkCostForT2(Similar_cost_matrix, Compact_cost_matrix, all_ans, all_labels, N, all_cnt, square_len, clusters_cnt, alpha, beta, true, false, old_grid_asses, old_T_pair, old_D_pair);
        for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = all_ans[tmp_gid];
    }
    else if(type=="B")
        cost = checkCostForB(Similar_cost_matrix, Compact_cost_matrix, all_ans, all_labels, N, all_cnt, square_len, clusters_cnt, alpha, beta);
    else if(type=="Dev")
        cost = checkCostForDev(Similar_cost_matrix, Compact_cost_matrix, all_ans, all_labels, N, all_cnt, square_len, clusters_cnt, alpha, beta);
    else if(type=="AlphaT")
        cost = checkCostForAlphaT(Similar_cost_matrix, Compact_cost_matrix, all_ans, all_labels, N, all_cnt, square_len, clusters_cnt, alpha, beta);
    else if(type=="MoS")
        cost = checkCostForMoS(Similar_cost_matrix, Compact_cost_matrix, all_ans, all_labels, N, all_cnt, square_len, clusters_cnt, alpha, beta);
        // printf("cost %.6lf %.6lf %.6lf\n", cost[1], cost[2], cost[3]);

    ret.push_back(cost[3]);

    delete[] all_ans;
    delete[] all_labels;

    delete[] if_disconn;
    delete[] checked;

    delete[] old_grid_asses;
    delete[] old_T_pair;
    delete[] old_D_pair;

    delete[] grid_asses;
    delete[] ori_embedded;
    delete[] cluster_labels;
    delete[] Similar_cost_matrix;
    delete[] Compact_cost_matrix;

    return ret;
}

// check convexity one cluster
std::vector<double> checkCostForAll(
const std::vector<std::vector<double>> &_ori_embedded,
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels,
const std::string &type,
double alpha, double beta) {

    // ----------------------------------preprocess step start----------------------------------------

    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = max(maxLabel, _cluster_labels[i]+1);
    int *grid_asses = new int[N];
//    int *ori_grid_asses = new int[N];
    double (*ori_embedded)[2] = new double[N][2];
    int *cluster_labels = new int[num];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
//    for(int i=0;i<N;i++)ori_grid_asses[i] = _ori_grid_asses[i];
    for(int i=0;i<N;i++) {
        ori_embedded[i][0] = _ori_embedded[i][0];
        ori_embedded[i][1] = _ori_embedded[i][1];
    }
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];

    double *Similar_cost_matrix = new double[N*N];
    getOriginCostMatrixArrayToArray(ori_embedded, cluster_labels, Similar_cost_matrix, N, num, square_len, maxLabel);

    double *Compact_cost_matrix = new double[N*N];
    getCompactCostMatrixArrayToArray(grid_asses, cluster_labels, Compact_cost_matrix, N, num, square_len, maxLabel);

    int *ans = new int[N];
    for(int i=0;i<N;i++)ans[i] = grid_asses[i];
    double best = 2147483647;    // cost of ans

    double *old_T_pair = new double[2];    // convexity of triples
    double *old_D_pair = new double[N*N*2];
    int *old_grid_asses = new int[N];
    for(int i=0;i<N;i++)old_grid_asses[i] = grid_asses[i];

    printf("final cost\n");
    std::vector<double> cost(4, -1);
    if(type=="CutRatio")
        cost = checkCostForE(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="AreaRatio")
        cost = checkCostForS(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="PerimeterRatio")
        cost = checkCostForC(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="2020")
        cost = checkCostFor2020(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="TB")
        cost = checkCostForTB(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="Triple") {
        cost = checkCostForT(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, false, old_grid_asses, old_T_pair);
        for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = ans[tmp_gid];
    }
    else if(type=="T2") {
        cost = checkCostForT2(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta, true, false, old_grid_asses, old_T_pair, old_D_pair);
        for(int tmp_gid=0;tmp_gid<N;tmp_gid++)old_grid_asses[tmp_gid] = ans[tmp_gid];
    }
    else if(type=="B")
        cost = checkCostForB(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="Dev")
        cost = checkCostForDev(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="Global")
        cost = checkCostForGlobal(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="AlphaT")
        cost = checkCostForAlphaT(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);
    else if(type=="MoS")
        cost = checkCostForMoS(Similar_cost_matrix, Compact_cost_matrix, ans, cluster_labels, N, num, square_len, maxLabel, alpha, beta);

    std::vector<double> ret(N+4, 0);
    for(int i=0;i<N;i++)ret[i] = ans[i];
    for(int i=0;i<3;i++)ret[N+i] = cost[i+1];
    int* checked = new int[N];
    double tmp = N*checkConnectForAll(ans, cluster_labels, checked, N, num, square_len, maxLabel);
    ret[N+3] = tmp;
    delete[] checked;

    delete[] old_grid_asses;
    delete[] old_T_pair;
    delete[] old_D_pair;

    delete[] grid_asses;
//    delete[] ori_grid_asses;
    delete[] ori_embedded;
    delete[] cluster_labels;
    delete[] Similar_cost_matrix;
    delete[] Compact_cost_matrix;
    delete[] ans;

    return ret;
}

// bi-graph match in each cluster to ensure the minist Similar cost
std::vector<int> optimizeInnerCluster(
//const std::vector<int> &_ori_grid_asses,
const std::vector<std::vector<double>> &_ori_embedded,
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels,
const std::vector<bool> &_change) {
    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = max(maxLabel, _cluster_labels[i]+1);
    int *grid_asses = new int[N];
//    int *ori_grid_asses = new int[N];
    double (*ori_embedded)[2] = new double[N][2];
    int *cluster_labels = new int[num];
    bool *change = new bool[N];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
//    for(int i=0;i<N;i++)ori_grid_asses[i] = _ori_grid_asses[i];
    for(int i=0;i<N;i++) {
        ori_embedded[i][0] = _ori_embedded[i][0];
        ori_embedded[i][1] = _ori_embedded[i][1];
    }
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];
    for(int i=0;i<N;i++)change[i] = _change[i];

    double *Similar_cost_matrix = new double[N*N];
    getOriginCostMatrixArrayToArray(ori_embedded, cluster_labels, Similar_cost_matrix, N, num, square_len, maxLabel);

    double *cost_matrix = new double[N*N];    //cluster内部各自进行二分图匹配，代价矩阵
    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid2=0;gid2<N;gid2++){
        int lb2 = -1;
        int id2 = grid_asses[gid2];
        if(id2<num)lb2 = cluster_labels[id2];
        for(int gid1=0;gid1<N;gid1++){
            int lb1 = -1;
            int id1 = grid_asses[gid1];
            if(id1<num)lb1 = cluster_labels[id1];
            if(lb1==lb2){
                cost_matrix[gid2*N+id1] = Similar_cost_matrix[gid2*N+id1];
            }else {
                cost_matrix[gid2*N+id1] = N;
            }
        }
    }

    std::vector<int> ret = solveBiMatchChange(cost_matrix, N, change, grid_asses);   //cluster内部各自进行二分图匹配

    delete[] cost_matrix;
    delete[] grid_asses;
//    delete[] ori_grid_asses;
    delete[] ori_embedded;
    delete[] cluster_labels;
    delete[] change;
    delete[] Similar_cost_matrix;

    return ret;
}

// bi-graph match in each cluster to ensure the minist Similar cost, with must-link
std::vector<int> optimizeInnerClusterWithMustLink(
//const std::vector<int> &_ori_grid_asses,
const std::vector<std::vector<double>> &_ori_embedded,
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels,
const std::vector<bool> &_change,
const std::vector<std::vector<int>> &_must_links,
const std::vector<std::vector<int>> &_must_links2,
const int maxit = 1) {
    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = max(maxLabel, _cluster_labels[i]+1);
    int *grid_asses = new int[N];
//    int *ori_grid_asses = new int[N];
    double (*ori_embedded)[2] = new double[N][2];
    int *cluster_labels = new int[num];
    bool *change = new bool[N];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
//    for(int i=0;i<N;i++)ori_grid_asses[i] = _ori_grid_asses[i];
    for(int i=0;i<N;i++) {
        ori_embedded[i][0] = _ori_embedded[i][0];
        ori_embedded[i][1] = _ori_embedded[i][1];
    }
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];
    for(int i=0;i<N;i++)change[i] = _change[i];

    double *Similar_cost_matrix = new double[N*N];
    getOriginCostMatrixArrayToArray(ori_embedded, cluster_labels, Similar_cost_matrix, N, num, square_len, maxLabel);

    double *cost_matrix = new double[N*N];    //cluster内部各自进行二分图匹配，代价矩阵
    int ml_N = _must_links.size();
    int ml_N2 = _must_links2.size();
    int *element_asses = new int[N];
    double *cluster_dist = new double[N*(maxLabel+1)];
    bool *is_border = new bool[N];
    for(int id=0;id<N;id++)is_border[id] = false;
    for(int i=0;i<ml_N2;i++)is_border[_must_links2[i][0]] = true;


    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid2=0;gid2<N;gid2++){
        int bias = gid2*(maxLabel+1);
        int lb2 = maxLabel;
        int id2 = grid_asses[gid2];
        if(id2<num)lb2 = cluster_labels[id2];
        int x2 = gid2/square_len;
        int y2 = gid2%square_len;
        for(int lb=0;lb<=maxLabel;lb++)cluster_dist[bias+lb] = 2147483647;
        for(int gid1=0;gid1<N;gid1++){
            int lb1 = maxLabel;
            int id1 = grid_asses[gid1];
            if(id1<num)lb1 = cluster_labels[id1];
            int x1 = gid1/square_len;
            int y1 = gid1%square_len;
            cluster_dist[bias+lb1] = min(cluster_dist[bias+lb1], getDist(x1, y1, x2, y2));
        }
    }

    for(int it=0;it<maxit;it++) {
        #pragma omp parallel for num_threads(THREADS_NUM)
        for(int gid2=0;gid2<N;gid2++){
            int bias = gid2*N;
            int lb2 = maxLabel;
            int id2 = grid_asses[gid2];
            element_asses[id2] = gid2;
            if(id2<num)lb2 = cluster_labels[id2];
            for(int gid1=0;gid1<N;gid1++){
                int lb1 = maxLabel;
                int id1 = grid_asses[gid1];
                if(id1<num)lb1 = cluster_labels[id1];
                if(lb1==lb2){
                    cost_matrix[bias+id1] = Similar_cost_matrix[bias+id1];
                }else {
                    cost_matrix[bias+id1] = N;
                }
            }
        }

        for(int i=0;i<ml_N;i++) {
            int ml_size = _must_links[i].size();
            double x = 0;
            double y = 0;
            double tot_weight = 0;
            for(int j=0;j<ml_size;j++) {
                int id = _must_links[i][j];
                double weight = 1;
                if(is_border[id])weight = 2*(ml_size-1);
                int gid = element_asses[id];
                x += weight*(gid/square_len);
                y += weight*(gid%square_len);
                tot_weight += weight;
            }
            x /= tot_weight;
            y /= tot_weight;
            #pragma omp parallel for num_threads(THREADS_NUM)
            for(int j=0;j<ml_size;j++) {
                int id = _must_links[i][j];
                for(int gid=0;gid<N;gid++) {
                    int x0 = gid/square_len;
                    int y0 = gid%square_len;
                    cost_matrix[gid*N+id] += 1 * ((x0-x)*(x0-x) + (y0-y)*(y0-y));
                }
            }
        }

        for(int i=0;i<ml_N2;i++) {
            int id = _must_links2[i][0];
            int lb = _must_links2[i][1];
            for(int gid=0;gid<N;gid++){
                int tmp = gid*(maxLabel+1)+lb;
                cost_matrix[gid*N+id] += 10*cluster_dist[tmp]*cluster_dist[tmp];
            }
        }

        std::vector<int> ret = solveBiMatchChange(cost_matrix, N, change, grid_asses);   //cluster内部各自进行二分图匹配
        for(int i=0;i<N;i++)grid_asses[i] = ret[i];
    }
    std::vector<int> ret(N, 0);
    for(int i=0;i<N;i++)ret[i] = grid_asses[i];

    delete[] cost_matrix;
    delete[] grid_asses;
//    delete[] ori_grid_asses;
    delete[] ori_embedded;
    delete[] cluster_labels;
    delete[] change;
    delete[] Similar_cost_matrix;
    delete[] element_asses;
    delete[] cluster_dist;
    delete[] is_border;

    return ret;
}

void testomp() {
    int sum=0;
    #pragma omp parallel for reduction(+:sum) num_threads(8)
    for(int i=0;i<1000;i++){
        sum += i;
        if(i%50==0)printf("%d %d\n", i, omp_get_thread_num());
    }
    printf("%d\n",sum);
}

int testrand() {
    srand(10);
    return rand();
}

void testmem() {
    srand(10);
    int* a = new int[5];
    memset(a, 0, 5*sizeof(int));
    for(int i=0;i<5;i++)printf("%d\n", a[i]);
    int* b = new int[6];
    memset(b, 0, 6*sizeof(int));
    memset(b, 2, 4*sizeof(int));
    memcpy(b, a, 3*sizeof(int));
    for(int i=0;i<6;i++)printf("%d\n", b[i]);

    delete[] a;
    delete[] b;
    return;
}

PYBIND11_MODULE(gridlayoutOpt, m) {
    m.doc() = "Gridlayout Optimizer"; // optional module docstring
    m.def("getClusters", &getClusters, "A function");
    m.def("solveKM", &solveKM, "A function");
    m.def("solveLap", &solveLap, "A function");
    m.def("getConnectCostMatrix", &getConnectCostMatrix, "A function to get cost matrix");
    m.def("getCompactCostMatrix", &getCompactCostMatrix, "A function to get cost matrix");
//    m.def("checkConvexForE", &checkConvexForE, "A function to check convexity");
//    m.def("checkConvexForT", &checkConvexForT, "A function to check convexity");
//    m.def("getCostMatrixForE", &getCostMatrixForE, "A function to get cost matrix");
//    m.def("getCostMatrixForT", &getCostMatrixForT, "A function to get cost matrix");
    m.def("optimizeBA", &optimizeBA, "A function to optimize");
    m.def("optimizeSwap", &optimizeSwap, "A function to optimize");
    m.def("checkCostForAll", &checkCostForAll, "A function to check cost");
    m.def("checkCostForOne", &checkCostForOne, "A function to check cost");
    m.def("checkCostForAllShapes", &checkCostForAllShapes, "A function to check cost");
    m.def("optimizeInnerCluster", &optimizeInnerCluster, "A function to optimize");
    m.def("optimizeInnerClusterWithMustLink", &optimizeInnerClusterWithMustLink, "A function to optimize");
    m.def("find_alpha", &find_alpha, "A function");
    m.def("testomp", &testomp, "A function");
    m.def("testrand", &testrand, "A function");
    m.def("testmem", &testmem, "A function");
}