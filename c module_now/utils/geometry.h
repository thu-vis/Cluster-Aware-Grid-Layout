#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

typedef std::vector<double> Point, Vector;
typedef std::vector<std::vector<double>> PointList, VectorList;

struct pair_hash
{
    template<class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const
    {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

// Geometric calculation of foundation
void getSurroundRect(
    double x1, double y1, double x2, double y2,
    double &xmin, double &xmax, double &ymin, double &ymax)
{
    if(x1 < x2) {
        xmin = x1; xmax = x2;
    }
    else {
        xmin = x2; xmax = x1;
    }
    if(y1 < y2) {
        ymin = y1; ymax = y2;
    }
    else {
        ymin = y2; ymax = y1;
    }
}

bool judgeSegmentsIntersect(
    double x1, double y1, double x2, double y2,
    double x3, double y3, double x4, double y4)
{
    double xmin1, xmax1, ymin1, ymax1;
    double xmin2, xmax2, ymin2, ymax2;
    getSurroundRect(x1, y1, x2, y2, xmin1, xmax1, ymin1, ymax1);
    getSurroundRect(x3, y3, x4, y4, xmin2, xmax2, ymin2, ymax2);
	if(xmax1 < xmin2 || xmax2 < xmin1 ||ymax1 < ymin2 || ymax2 < ymin1)   return false;

	int res1 = (x1 - x2) * (y3 - y2) - (y1 - y2) * (x3 - x2);
	int res2 = (x1 - x2) * (y4 - y2) - (y1 - y2) * (x4 - x2);

	int res3 = (x3 - x4) * (y1 - y4) - (y3 - y4) * (x1 - x4);
	int res4 = (x3 - x4) * (y2 - y4) - (y3 - y4) * (x2 - x4);
	if(res1 * res2 <= 0 && res3 * res4 <= 0) return true;
	else return false;
}

double getAngle(const std::vector<double> & vc1, const std::vector<double> vc2) {
    // calculate Angle in two vectors [-pi, pi]
    double dot = vc1[0] * vc2[0] + vc1[1] * vc2[1];
    double l1 = sqrt(vc1[0] * vc1[0] + vc1[1] * vc1[1]);
    double l2 = sqrt(vc2[0] * vc2[0] + vc2[1] * vc2[1]);
    double angle = acos(dot / (l1 * l2));
    // std::cout << angle << " "<< l1 << " "<<l2 <<"\n";
    double skew_product = vc1[0] * vc2[1] - vc1[1] * vc2[0];
    if(skew_product >= 0) return angle;
    else return -1 * angle;
}

double getPerimeter(const std::vector<std::vector<double>> & vcs, const int type = 0) {
    double res = 0;
    if(type == 0) {
        for(auto i = 1;i < vcs.size(); i++) {
            auto & vc1 = vcs[i];
            auto & vc2 = vcs[i-1];
            res += sqrt((vc1[0] - vc2[0]) * (vc1[0] - vc2[0]) + (vc1[1] - vc2[1]) * (vc1[1] - vc2[1]));
        }
        auto & vc1 = vcs[0];
        auto & vc2 = vcs[vcs.size()-1];
        res += sqrt((vc1[0] - vc2[0]) * (vc1[0] - vc2[0]) + (vc1[1] - vc2[1]) * (vc1[1] - vc2[1]));
    }
    else {
        for(auto & vc: vcs) {
            res += sqrt(vc[0] * vc[0] + vc[1] * vc[1]);
        }
    }
    return res;
}

double getArea(const std::vector<std::vector<double>> & polygon) {
    if(polygon.size() <= 2) return 0;
    double area = 0;
    for(int i=0;i < polygon.size() - 1; i++){
        double triArea = (polygon[i][0]*polygon[i+1][1] - polygon[i+1][0]*polygon[i][1]) / 2;
        area += triArea;
    }
    double fn = (polygon[polygon.size()-1][0]*polygon[0][1] - polygon[0][0]*polygon[polygon.size()-1][1])/2;
    return abs(area+fn);
}

double getVectorLength(const std::vector<double> vc){
    return sqrt(vc[0] * vc[0] + vc[1] * vc[1]);
}

// Judgment of simple polygon
bool judgePolygonSimple(
    const PointList & polygon,
    VectorList & rt_vectors,
    bool judge_selfinter = false)
{
    // judge non-coincide, non-reverse
    rt_vectors.clear();
    Vector before {0, 0};
    for (auto i = 1; i < polygon.size(); i++) {
        Vector temp {polygon[i][0] - polygon[i-1][0], polygon[i][1] - polygon[i-1][1]};
        if(temp[0] == 0 && temp[1] == 0) {
            std::cout << "JudgePolygonSimple - polygon coincides." << std::endl;
            return false;
        }
        if((temp[0] * before[1] == temp[1] * before[0]) && temp[0] * before[0] < 0) {
            std::cout << "JudgePolygonSimple - polygon reverses." << std::endl;
            return false;
        }
        rt_vectors.push_back(temp);
        before = temp;
    }
    rt_vectors.push_back(Vector {polygon[0][0] - polygon[polygon.size()-1][0], polygon[0][1] - polygon[polygon.size()-1][1]});

    // judge simple: has no self-intersections
    if(judge_selfinter) {
        for (auto i = 2; i < polygon.size(); i++) {
            for (auto j = 1; j < i - 1; j++) {
                if(judgeSegmentsIntersect(polygon[i][0], polygon[i][1], polygon[i-1][0], polygon[i-1][1],
                    polygon[j][0], polygon[j][1], polygon[j-1][0], polygon[j-1][1])) {
                    // std::cout<< polygon[i][0]<<" " <<polygon[i][1]<<" "<<polygon[i-1][0]<<" " <<polygon[i-1][0]<<std::endl;
                    // std::cout<< polygon[j][0]<<" " <<polygon[j][1]<<" "<<polygon[j-1][0]<<" " <<polygon[j-1][0]<<std::endl;
                    std::cout << "JudgePolygonSimple - polygon intersects." << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}


// polygon and grid transitions
static int grid_square = 1;
static int grid(int x, int y) {
    return x * grid_square + y;
}
static int point_square = 1;
static int point(int x, int y) {
    return x * point_square + y;
}
static std::vector<double> pos(int point){
    // std::cout << point <<" "<<point_square << " ";
    int x = point / point_square;
    int y = point % point_square;
    // std::cout << x <<" "<<y << " "<< 1.0 *x << " "<< 1.0 *y << std::endl;
    return std::vector<double> {1.0 * x, 1.0 * y};
}

void getClusterBoundary(
    const std::vector<int> & grid_asses,
    const std::vector<int> & cluster_labels,
    std::vector<PointList> & rt_boundary,
    std::vector<int> & rt_boundary_labels)
{
    // get cluster boundary segs
    int cur_id = 0;
    std::unordered_map<int, int> label2ids;
    std::vector<std::unordered_set<std::pair<int, int>, pair_hash>> edges;
    std::vector<std::unordered_map<int, std::vector<int>>> matrix;
    int square_len = ceil(sqrt(grid_asses.size()));
    grid_square = square_len;
    point_square = square_len + 1;

    for(auto i = 0;i < grid_asses.size(); i++) {
        int gid = grid_asses[i];
        int label = cluster_labels[gid];
        if(label2ids.find(label) == label2ids.end()){
            label2ids[label] = cur_id;
            edges.push_back(std::unordered_set<std::pair<int, int>, pair_hash>());
            matrix.push_back(std::unordered_map<int, std::vector<int>>());
            cur_id += 1;
        }
        int lid = label2ids[label];
        auto & ledges = edges[lid];
        auto & lmatrix = matrix[lid];

        int x = i / square_len;
        int y = i % square_len;
        int edge_start, edge_end;
        if(x == 0 || cluster_labels[grid_asses[grid(x-1, y)]] != label) {
            edge_start = point(x, y);
            edge_end = point(x, y+1);
            ledges.insert(std::pair<int, int>(edge_start, edge_end));
            lmatrix[edge_start].push_back(edge_end);
            lmatrix[edge_end].push_back(edge_start);
        }
        if(y == 0 || cluster_labels[grid_asses[grid(x, y-1)]] != label) {
            edge_start = point(x, y);
            edge_end = point(x+1, y);
            ledges.insert(std::pair<int, int>(edge_start, edge_end));
            lmatrix[edge_start].push_back(edge_end);
            lmatrix[edge_end].push_back(edge_start);
        }
        if(x == square_len-1 || cluster_labels[grid_asses[grid(x+1, y)]] != label) {
            edge_start = point(x+1, y);
            edge_end = point(x+1, y+1);
            ledges.insert(std::pair<int, int>(edge_start, edge_end));
            lmatrix[edge_start].push_back(edge_end);
            lmatrix[edge_end].push_back(edge_start);
        }
        if(y == square_len-1 || cluster_labels[grid_asses[grid(x, y+1)]] != label) {
            edge_start = point(x, y+1);
            edge_end = point(x+1, y+1);
            ledges.insert(std::pair<int, int>(edge_start, edge_end));
            lmatrix[edge_start].push_back(edge_end);
            lmatrix[edge_end].push_back(edge_start);
        }
    }

    // connect boundary
    rt_boundary.clear();
    rt_boundary_labels.clear();
    for(auto it = label2ids.begin(); it != label2ids.end(); it++) {
        int label = it->first;
        int lid = it->second;
        auto & ledges = edges[lid];
        auto & lmatrix = matrix[lid];
        while(!ledges.empty()) {
            auto edge_head = *ledges.begin();
            ledges.erase(ledges.begin());
            PointList points;
            points.push_back(pos(edge_head.first));
            int start = edge_head.first;
            int now = edge_head.second;
            int bf = start;
            while(now != start) {
                auto & link = lmatrix[now];
                int nxt = 0;
                if(link[0] == bf) {
                    nxt = link[1];
                }
                else nxt = link[0];
                if(now < nxt) ledges.erase(std::pair<int, int>(now, nxt));
                else ledges.erase(std::pair<int, int>(nxt, now));
                points.push_back(pos(now));
                bf = now;
                now = nxt;
            }
            rt_boundary.push_back(points);
            rt_boundary_labels.push_back(label);
        }
    }
}



#endif