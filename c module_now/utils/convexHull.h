#ifndef _CONVEXHULL_H
#define _CONVEXHULL_H

#include <iostream>
#include <vector>
#include <utility>
#include <ctime>
#include <algorithm>
#include <math.h>
#include "base.h"

struct CH_node {
    int id;
    double x;
    double y;
    double dist;
};

double Distance(const double nodes[][2], int a, int b)
{
    return std::sqrt((nodes[a][0]-nodes[b][0])*(nodes[a][0]-nodes[b][0])+(nodes[a][1]-nodes[b][1])*(nodes[a][1]-nodes[b][1]));
}

double Multiply0(const struct CH_node &b, const struct CH_node &c)
{
    return(b.x*c.y-b.y*c.x);
}

double Multiply(const struct CH_node &a, const struct CH_node &b, const struct CH_node &c)
{
    return((b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x));
}

bool cmp(const struct CH_node &x, const struct CH_node &y){
    double m;
    m = Multiply0(x, y);
    if(m>0)return true;
    else if(m==0&&(x.dist<y.dist))
        return true;
    else return false;
}


// N: num of vertexes
// nodes[0..N-1][2]: vertex
// return int M，it's the num of vertexes of convex hull
// nodes will be modified，nodes[0..M-1][2] is the vertexes of convex hull in clockwise order
int getConvexHull(int N, double nodes[][2]){
    if(N<=2)return N;
    
    double px, py;
    int p;
    py=-1;
    for(int i=0;i<N;i++)
    {
        if(py==-1||nodes[i][1]<py)
        {
            px=nodes[i][0];
            py=nodes[i][1];
            p=i;
        }
        else if(nodes[i][1]==py&&nodes[i][0]<px)
        {
            px=nodes[i][0];
            py=nodes[i][1];
            p=i;
        }
    }
    double tmp=nodes[0][0]; nodes[0][0]=nodes[p][0]; nodes[p][0]=tmp;
    tmp=nodes[0][1]; nodes[0][1]=nodes[p][1]; nodes[p][1]=tmp;

    double (*tmp_nodes)[2] = new double[N+1][2];
    for(int i=0;i<N;i++){
        tmp_nodes[i][0] = nodes[i][0];
        tmp_nodes[i][1] = nodes[i][1];
    }
    tmp_nodes[N][0] = tmp_nodes[0][0]; tmp_nodes[N][1] = tmp_nodes[0][1];

    struct CH_node *a = new struct CH_node[N+1];

    for(int i=0;i<N;i++){
        a[i].id = i;
        a[i].x = nodes[i][0]-nodes[0][0];
        a[i].y = nodes[i][1]-nodes[0][1];
        a[i].dist = getDist(a[i].x, a[i].y, 0, 0);
    }
    a[N].id = N; a[N].x = 0; a[N].y = 0; a[N].dist = 0;

    std::sort(a+1, a+N, cmp);

    int* order = new int[N+1];
    int* order2 = new int[N+1];
    order2[0] = 0;
    order2[1] = 1;
    order2[2] = 2;

    int top=2;
    for(int ii=3;ii<=N;ii++)
    {
        while(top>=1&&Multiply(a[order2[top-1]], a[order2[top]], a[ii])<=0)
            top--;
        order2[top+1] = ii;
        top++;
    }

    for(int i=0;i<top;i++){
        order[i] = a[order2[i]].id;
        nodes[i][0] = tmp_nodes[order[i]][0];
        nodes[i][1] = tmp_nodes[order[i]][1];
    }

    delete[] order;
    delete[] order2;
    delete[] tmp_nodes;
    delete[] a;

    return top;
}

// get the area of convex hull
double getSofPoly(int N, double nodes[][2]){
    if(N<=2)return 0;
    double area = 0;
    for(int i=0;i<N-1;i++){
        double triArea = (nodes[i][0]*nodes[i+1][1] - nodes[i+1][0]*nodes[i][1])/2;
        area += triArea;
    }
    double fn = (nodes[N-1][0]*nodes[0][1] - nodes[0][0]*nodes[N-1][1])/2;
    return abs(area+fn);
}

// get the perimeter of convex hull
double getCofPoly(int N, double nodes[][2]){
    if(N<=1)return 0;
    double C = 0;
    for(int i=0;i<N-1;i++){
        double edge = Distance(nodes, i, i+1);
        C += edge;
    }
    double edge = Distance(nodes, N-1, 0);
    return abs(C+edge);
}

#endif