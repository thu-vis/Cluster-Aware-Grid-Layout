#include <iostream>
#include <vector>
#include <utility>
#include <ctime>
#include <algorithm>
#include <math.h>
#include <string.h>

const int zkw_INF = 2147483647;
const double zkw_tiny = 0.000001;

void zkw_link(int i, int j, int k, double l, int &cnt, int next[], int head[], int p[], int d[], double w[])
{
    cnt++; next[cnt]=head[i]; head[i]=cnt; p[cnt]=j; d[cnt]=k; w[cnt]=l;
    cnt++; next[cnt]=head[j]; head[j]=cnt; p[cnt]=i; d[cnt]=0; w[cnt]=-l;
}

int zkw_dfs(int tt, double &ans, int x, int flow, double dis[], int o[], int v[], int arc[], int next[], int p[], int d[], double w[], double f[])
{
    if(x==tt){ ans+=dis[x]*flow; return flow; }
    int ans2=0; o[x]=1; v[x]=1;
    for(int i=arc[x];i;i=next[i])
    {
        int j = p[i];
        double l = dis[x]+w[i]-dis[j];
        int mi = std::min(flow,d[i]);

        if((mi>0)&&(l<f[j]))f[j] = l;
        if((mi>0)&&(l<zkw_tiny)&&(!o[j]))
        {
            l = zkw_dfs(tt, ans, j, mi, dis, o ,v, arc, next, p, d, w, f);
            ans2 += l; flow -= l; d[i] -= l; d[i^1] += l;
        }
        if(flow==0)break; arc[x] = next[arc[x]];
    }
    o[x]=0;
    return ans2;
}

double zkw(float edge[], int N, int ans_asses[])
{
    int *o = new int[N*2+5];
    int *v = new int[N*2+5];
    double *f = new double[N*2+5];
    int *head = new int[N*2+5];
    int *arc = new int[N*2+5];
    double *dis = new double[N*2+5];

    for(int i=0;i<N*2+5;i++){
        o[i] = 0; v[i] = 0; f[i] = 0;
        head[i] = 0; arc[i] = 0; dis[i] = 0;
    }

    int edge_N = 0;
    for(int i=0;i<N*N;i++)if(edge[i]>-0.5)edge_N += 1;
    edge_N += 2*N;

    int *next = new int[edge_N*2+5];
    int *p = new int[edge_N*2+5];
    int *d = new int[edge_N*2+5];
    double *w = new double[edge_N*2+5];

    for(int i=0;i<edge_N*2+5;i++){
        next[i] = 0; p[i] = 0; d[i] = 0; w[i] = 0;
    }

    int cnt = 1;

    for(int i=0;i<N;i++)zkw_link(2*N, i, 1, 0, cnt, next, head, p, d, w);
    for(int i=0;i<N;i++)zkw_link(N+i, 2*N+1, 1, 0, cnt, next, head, p, d, w);

    for(int i=0;i<N;i++){
        int bias = i*N;
        for(int j=0;j<N;j++){
            if(edge[bias+j]>-zkw_tiny){
                zkw_link(i, N+j, 1, edge[bias+j], cnt, next, head, p, d, w);
            }
        }
    }

    int flow = 0;
    double ans = 0;
    while(1)
    {
        for(int i=0;i<=2*N+1;i++)v[i] = 0;
        for(int i=0;i<=2*N+1;i++)arc[i] = head[i];
        for(int i=0;i<=2*N+1;i++)f[i] = zkw_INF;
        flow += zkw_dfs(2*N+1, ans, 2*N, zkw_INF, dis, o, v, arc, next, p, d, w, f);
        double imp = zkw_INF;
        for(int i=0;i<=2*N+1;i++)
            if((!v[i])&&(f[i]<imp))imp = f[i];
        if(imp>zkw_INF/2)break;
        for(int i=0;i<=2*N+1;i++)if(!v[i])dis[i] += imp;
    }

    for(int x=0;x<N;x++){
        for(int i=head[x];i;i=next[i]){
            int j = p[i];
            if((d[i]==0)&&(N<=j)&&(j<2*N)&&((i&1)==0)){
                ans_asses[x] = j-N;
                break;
            }
        }
    }
    return ans;
}