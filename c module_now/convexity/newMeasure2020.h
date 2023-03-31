#ifndef _MEASURE_2020_H
#define _MEASURE_2020_H

#include "../utils/base.h"
#include "../utils/convexHull.h"
#include "../utils/util.h"
#include <functional>
#include <float.h>

const double lambda = 1;    // tuning parameter
const int num_per_grid = 15; // 计算积分时，每个格子采样数为 num_per_grid ** 2
// 将 num_per_grid 调得越大，计算结果越精确

// 用于 debug. 输出图形格点坐标
void outputPolygon(int n, double grids[][2]) {
    std::cout << "point number: " << n << '\n';
    for (int i = 0; i < n; ++i) {
        std::cout << grids[i][0] << ", " << grids[i][1] << '\n';
    }
}

std::function<double(double, double)> psi(double x1, double y1, double x2, double y2) {
    return [x1, y1, x2, y2](double x, double y) {
        return (y2*x - y1*x + x1*y - x2*y) / (x1*y2 - x2*y1);
    };
}

std::function<double(double, double)> PsiFunction(int n, double grids[][2]) {
    return [n, grids](double x, double y) {
        auto func = psi(grids[n-1][0], grids[n-1][1], grids[0][0], grids[0][1]);
        double res = func(x, y);
        for (int i = 0; i < n-1; ++i) {
            func = psi(grids[i][0], grids[i][1], grids[i+1][0], grids[i+1][1]);
            res = std::max(func(x, y), res);
        }
        return std::pow(res, lambda);
    };
}

std::vector<double> getAreaAndCentroid(int n, double grids[][2]) {
    // get area & centroid of a convex polygon
    // split the polygon into several triangles
    // calculate area & centroid for all the triangles
    // and combine them together as the result

    double area = 0;
    double cx = 0;  // centroid.x()
    double cy = 0;  // centroid.y()

    if (n >= 3) {
        double (*copy)[2] = new double[n-1][2];
        for (int i = 0; i < n-1; ++i) {
            // translate grids[0] to be coincide with (0, 0)
            // for the sake of further calculation
            copy[i][0] = grids[i+1][0] - grids[0][0];
            copy[i][1] = grids[i+1][1] - grids[0][1];
        }
        double* areas = new double[n-2];
        for (int i = 0; i < n-2; ++i) {
            areas[i] = std::abs(copy[i][0] * copy[i+1][1] - copy[i][1] * copy[i+1][0]) / 2;
            area += areas[i];
            copy[i][0] = (copy[i][0] + copy[i+1][0]) / 3;
            copy[i][1] = (copy[i][1] + copy[i+1][1]) / 3;
            cx += areas[i] * copy[i][0];
            cy += areas[i] * copy[i][1];
        }
        cx /= area;
        cy /= area;
        delete [] copy;
        delete [] areas;
    }

    if (n == 2) {
        cx = (grids[0][0] + grids[1][0]) / 2;
        cy = (grids[0][1] + grids[1][1]) / 2;
    }

    std::vector<double> result;
    result.resize(3);
    result[0] = area;
    result[1] = cx;
    result[2] = cy;
    return result;
}

double checkConvex2020(
int grids[][2],
int grids_n,
bool save=false, double save_convex_hull[][2]=nullptr,
int save_CH_cnt[]=nullptr, double translate[][5]=nullptr) {
// grids_n: 格子个数
// grids[0~(grids_n-1)][2]: 每个格子的坐标
// 输出一个double 0~1，代表该图形的凸性，越大凸性越好

// 提醒1：convexHull.h中含有一些凸包相关的函数可供调用。
//        例如getConvexHull(N, nodes[][2])，该方法用于求解散点的凸包。
//        详见代码。不要修改convexHull.h中的代码。
// 提醒2：坐标为(i,j)的格子，实际范围应为(i,j)--(i,j+1)--(i+1,j+1)--(i+1,j)这个正方形。
//        如果想求整个图形的凸包，应当把所有格子的4个顶点都传入getConvexHull函数参数中。
// 提醒3：对于凸性计算式中的积分部分，可以采用均匀间隔枚举求和的方式。
//        例如，以1/5的枚举间隔，对于坐标为(0,0)的格子，
//        可以枚举[1/10, 3/10, 5/10, 7/10, 9/10] X [1/10, 3/10, 5/10, 7/10, 9/10]这25个点用于计算。
//        最好不要枚举边角上的点（x或y坐标为整数），这些点并不只在一个格子的范围内。
//        只枚举格子内部的点，可以使得每个格子的计算相对独立，方便方法拓展。
//        另外，因为枚举的精度问题，原图形和凸图形的面积虽然相同，在其范围内枚举到的点个数可能不同，可以通过乘一个倍数或分别取平均来相应地消除影响。
// 提醒4：一些可能有用的图形学算法：
//            扫描线法：用于枚举不规则区域中的所有点。
//            射线法：用于判断一个点是否在一个区域内。
//        希望该函数的复杂度控制在O(grids_n * k * m)，其中k为每个格子中枚举的点的个数，m为凸包的边数。

    if(grids_n==0)return 0;
    // STEP1: get centroid & area of the input shape
    int n = grids_n;
    double area = n;    // assume that a single grid has area of 1
    double centroid[2] = {0, 0};
    for (int i = 0; i < n; ++i) {
        centroid[0] += grids[i][0];
        centroid[1] += grids[i][1];
    }
    centroid[0] = centroid[0] / n + 0.5;
    centroid[1] = centroid[1] / n + 0.5;

    // STEP2: get convex hull
    double (*convex_hull)[2] = new double[4*n][2];
    for (int i = 0; i < n; ++i) {
        convex_hull[4*i][0] = grids[i][0];
        convex_hull[4*i][1] = grids[i][1];
        convex_hull[4*i+1][0] = grids[i][0]+1;
        convex_hull[4*i+1][1] = grids[i][1];
        convex_hull[4*i+2][0] = grids[i][0];
        convex_hull[4*i+2][1] = grids[i][1]+1;
        convex_hull[4*i+3][0] = grids[i][0]+1;
        convex_hull[4*i+3][1] = grids[i][1]+1;
    }
    int vertex_num = getConvexHull(4*n, convex_hull);
    // outputPolygon(vertex_num, convex_hull);

    // STEP3: get centroid & area of the convex hull
    auto res = getAreaAndCentroid(vertex_num, convex_hull);
    double hull_area = res[0];
    double hull_centroid[2];
    hull_centroid[0] = res[1];
    hull_centroid[1] = res[2];

    // STEP4: move & scale
    double dx = hull_centroid[0] + convex_hull[0][0];
    double dy = hull_centroid[1] + convex_hull[0][1];
    double scale_rate = std::sqrt(hull_area);
    for (int i = 0; i < vertex_num; ++i) {
        convex_hull[i][0] = (convex_hull[i][0] - dx) / scale_rate;
        convex_hull[i][1] = (convex_hull[i][1] - dy) / scale_rate;
    }

    if(save) {
        for(int i=0;i<vertex_num;i++) {
            save_convex_hull[i][0] = convex_hull[i][0];
            save_convex_hull[i][1] = convex_hull[i][1];
        }
        save_CH_cnt[0] = vertex_num;
        translate[0][0] = dx;
        translate[0][1] = dy;
        translate[0][2] = centroid[0];
        translate[0][3] = centroid[1];
        translate[0][4] = scale_rate;
    }

    scale_rate = std::sqrt(area);
    double (*moved_grids)[2] = new double[n][2];
    for (int i = 0; i < n; ++i) {
        moved_grids[i][0] = (grids[i][0] - centroid[0]) / scale_rate;
        moved_grids[i][1] = (grids[i][1] - centroid[1]) / scale_rate;
    }

    // STEP5: calculate Psi function & measure convexity
    auto Psi = PsiFunction(vertex_num, convex_hull);
    double p = 0;
    double gap = (1.0 / scale_rate) / (num_per_grid);
    double cur_x, cur_y;
    for (int i = 0; i < n; ++i) {
        for (int x_step = 0; x_step < num_per_grid; ++x_step) {
            for (int y_step = 0; y_step < num_per_grid; ++y_step) {
                cur_x = moved_grids[i][0] + (gap/2) + x_step * gap;
                cur_y = moved_grids[i][1] + (gap/2) + y_step * gap;
                p += Psi(cur_x, cur_y);
            }
        }
    }
    p = p / (num_per_grid * num_per_grid * n);  // normalize

    delete [] convex_hull;
    delete [] moved_grids;

    return std::min(2.0 / (2.0+lambda) / p, 1.0);
}

void getCostMatrixFor2020ArrayToArray(
const int grid_asses[],
const int cluster_labels[],
double cost_matrix_a[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

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

    int (*nodes)[2] = new int[num*4][2];

    int *CH_bias = new int[maxLabel];
    double (*convex_hull)[2] = new double[N][2];
    int *CH_cnt = new int[maxLabel];
    double (*translate)[5] = new double[maxLabel][5];
    double (*S_pairs)[2] = new double[maxLabel][2];
    double (*calc)[2] = new double[maxLabel][2];

    double S0 = 0, S1 = 0;
    int bias_cnt = 0;
    for(int li=0;li<maxLabel;li++){
        CH_bias[li] = bias_cnt;
        int cnt=0;
        double tmp_S0=0, tmp_S1=0;
        for(int id=head[li];id>=0;id=last[id]){
            int gid = element_asses[id];
            int x = gid/square_len;
            int y = gid%square_len;
            nodes[cnt][0] = x; nodes[cnt][1] = y; cnt++;
            tmp_S1 += 1;
        }
        double tmp = checkConvex2020(nodes, cnt, true, convex_hull+bias_cnt, CH_cnt+li, translate+li);
        tmp_S0 = tmp*tmp_S1;
        S_pairs[li][0] = tmp_S0;
        S_pairs[li][1] = tmp_S1;
        calc[li][1] = cnt*num_per_grid*num_per_grid;
        calc[li][0] = 2.0/(2.0+lambda)/tmp*calc[li][1];
        S1 += tmp_S1;
        S0 += tmp_S0;
        bias_cnt += cnt;
    }

    double *S0_matrix = new double[N*(maxLabel+1)];
    double *S1_matrix = new double[N*(maxLabel+1)];
    for(int i=0;i<N*(maxLabel+1);i++){
        S0_matrix[i] = S1_matrix[i] = -1;
    }

    #pragma omp parallel for num_threads(THREADS_NUM)
    for(int gid2=0;gid2<N;gid2++){
        int bias = gid2*N;
        int bias1 = gid2*(maxLabel+1);
        for(int gid1=0;gid1<N;gid1++) {
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
            if(S1_matrix[bias1+lb]<0){
                int lb2 = maxLabel;
                int id2 = grid_asses[gid2];
                if(id2<num)lb2 = cluster_labels[id2];
                if(lb==lb2) {
                    S0_matrix[bias1+lb] = S0;
                    S1_matrix[bias1+lb] = S1;
                }else {
                    S0_matrix[bias1+lb] = S0;
                    S1_matrix[bias1+lb] = S1;
                    if(lb2<maxLabel) {
                        S0_matrix[bias1+lb] -= S_pairs[lb2][0];
                        S1_matrix[bias1+lb] -= S_pairs[lb2][1];
                        auto Psi = PsiFunction(CH_cnt[lb2], convex_hull+CH_bias[lb2]);
                        double p = 0;
                        double scale_rate = translate[lb2][4];
                        double cx = translate[lb2][2];
                        double cy = translate[lb2][3];
                        double gap = (1.0 / scale_rate) / (num_per_grid);
                        double cur_x, cur_y;
                        for (int x_step = 0; x_step < num_per_grid; ++x_step) {
                            for (int y_step = 0; y_step < num_per_grid; ++y_step) {
                                cur_x = (x2-cx)/scale_rate + (gap/2) + x_step * gap;
                                cur_y = (y2-cy)/scale_rate + (gap/2) + y_step * gap;
                                p += Psi(cur_x, cur_y);
                            }
                        }
                        double new_p = calc[lb2][0]-p;
                        double new_p_cnt = calc[lb2][1]-num_per_grid*num_per_grid;
                        S1_matrix[bias1+lb] += S_pairs[lb2][1]-1;
                        S0_matrix[bias1+lb] += (S_pairs[lb2][1]-1)*std::min(2.0/(2.0+lambda)/(new_p/new_p_cnt), 1.0);
                    }
                    if(lb<maxLabel) {
                        S0_matrix[bias1+lb] -= S_pairs[lb][0];
                        S1_matrix[bias1+lb] -= S_pairs[lb][1];
                        auto Psi = PsiFunction(CH_cnt[lb], convex_hull+CH_bias[lb]);
                        double p = 0;
                        double scale_rate = translate[lb][4];
                        double cx = translate[lb][2];
                        double cy = translate[lb][3];
                        double gap = (1.0 / scale_rate) / (num_per_grid);
                        double cur_x, cur_y;
                        for (int x_step = 0; x_step < num_per_grid; ++x_step) {
                            for (int y_step = 0; y_step < num_per_grid; ++y_step) {
                                cur_x = (x2-cx)/scale_rate + (gap/2) + x_step * gap;
                                cur_y = (y2-cy)/scale_rate + (gap/2) + y_step * gap;
                                p += Psi(cur_x, cur_y);
                            }
                        }
                        double new_p = calc[lb][0]+p;
                        double new_p_cnt = calc[lb][1]+num_per_grid*num_per_grid;
                        S1_matrix[bias1+lb] += S_pairs[lb][1]+1;
                        S0_matrix[bias1+lb] += (S_pairs[lb][1]+1)*std::min(2.0/(2.0+lambda)/(new_p/new_p_cnt), 1.0);
                    }
                }
            }
            cost_matrix_a[bias+id1] = ((S0/S1)-(S0_matrix[bias1+lb]/S1_matrix[bias1+lb]))*N *2;
        }
    }

    delete[] head;
    delete[] last;
    delete[] element_asses;
    delete[] nodes;
    delete[] CH_bias;
    delete[] convex_hull;
    delete[] CH_cnt;
    delete[] translate;
    delete[] S_pairs;
    delete[] calc;
    delete[] S0_matrix;
    delete[] S1_matrix;
    return;
}

//检查全局凸性
//grid_asses[N]: 0~N-1，为当前布局中每个格子上放置的元素的标号
//cluster_labels[num]: 0~maxLabel-1，为每个元素的所属的类簇标号
//N: 格子个数
//num: 元素个数，标号num~N-1的元素为虚拟元素（空格）
//square_len: sqrt(N)，网格布局的宽高
//maxLabel: 类簇个数
//save, load...: 保存或读取类簇凸性的计算结果
std::vector<double> checkConvexFor2020Array(
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

    int (*nodes)[2] = new int[num*4][2];
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
            nodes[cnt][0] = x; nodes[cnt][1] = y; cnt++;
            tmp_S1 += 1;
        }
        tmp_S0 = checkConvex2020(nodes, cnt)*tmp_S1;
        S1 += tmp_S1;
        S0 += tmp_S0;
//        printf("%d %.2lf %.2lf\n", li, tmp_S0, tmp_S1);
        if(save){
            label_pairs[li][0] = tmp_S0;
            label_pairs[li][1] = tmp_S1;
        }
    }
    std::vector<double> S_pair(2, 0);
    S_pair[0] = (S1-S0);
    S_pair[1] = S1;
//    printf("%d %.2lf %.2lf\n", maxLabel, S_pair[0], S_pair[1]);

    delete[] head;
    delete[] last;
    delete[] element_asses;
    delete[] nodes;
    return S_pair;
}

//计算整体代价（与原layout相似性、紧凑性、凸性）
std::vector<double> checkCostFor2020(
    const double Similar_cost_matrix[],
    const double Compact_cost_matrix[],
    const int grid_asses[], const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    const double &alpha, const double &beta,
    bool save=false, bool load=false, double label_pairs[][2]=nullptr, int mainLabel1=-1, int mainLabel2=-1) {

    std::vector<double> S_pair = checkConvexFor2020Array(grid_asses, cluster_labels, N, num, square_len, maxLabel, save, load, label_pairs, mainLabel1, mainLabel2);
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