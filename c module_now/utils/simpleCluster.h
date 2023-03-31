#ifndef _SIMPLECLUSTER_H
#define _SIMPLECLUSTER_H

#include <iostream>
#include <vector>
#include <utility>
#include <ctime>

// using namespace std;

int find(int fa[], int x){
    if(fa[x]==x)return x;
    else{
        fa[x] = find(fa, fa[x]);
        return fa[x];
    }
}

int getClustersArray(
int cluster_labels[], int max_clusters, int min_size,
int grid_asses[],
int labels[],
const int &N, const int &num, const int &square_len){
    int *fa = new int[num];
    int *size = new int[num];
    int maxLabel = 0;
    for(int i=0;i<num;i++){
        fa[i] = i; size[i] = 1;
        cluster_labels[i] = -1;
    }
    for(int x1=0;x1<square_len;x1++){
        for(int y1=0;y1<square_len;y1++){
            int gid1 = x1*square_len+y1;
            int id1 = grid_asses[gid1];
            if(id1>=num)continue;
            int f_x[4] = {0, 1, 1, 1};
            int f_y[4] = {1, 0, 1, -1};
            for(int k=0;k<4;k++){
                int x2 = x1+f_x[k];
                int y2 = y1+f_y[k];
                if((x2<0)||(x2>=square_len)||(y2<0)||(y2>=square_len))continue;
                int gid2 = x2*square_len+y2;
                int id2 = grid_asses[gid2];
                if(id2>=num)continue;
                if(labels[id1]!=labels[id2])continue;
                int f1 = find(fa, id1);
                int f2 = find(fa, id2);
                if(f1==f2)continue;
                size[f1] += size[f2];
                fa[f2] = f1;
            }
        }
    }
    int *list = new int[N];
    int list_cnt = 0;
    for(int i=0;i<num;i++)
    if(fa[i]==i){
        list[list_cnt] = i;
        list_cnt++;
    }
    for(int i=0;i<list_cnt-1;i++)
    for(int j=i+1;j<list_cnt;j++)
    if(size[list[i]]<size[list[j]]){
        int t = list[i];
        list[i] = list[j];
        list[j] = t;
    }
    for(int i=0;i<std::min(max_clusters, list_cnt);i++){
        if(size[list[i]]<min_size)break;
        for(int id=0;id<num;id++)
        if(find(fa, id)==list[i])cluster_labels[id] = i;
        maxLabel = i+1;
    }
    int *votes = new int[maxLabel+1];
    int maxit = 5;
    for(int it=0;;it++) {
        int flag = 0;
        for(int x1=0;x1<square_len;x1++){
            for(int y1=0;y1<square_len;y1++){
                int gid1 = x1*square_len+y1;
                int id1 = grid_asses[gid1];
                if(id1>=num)continue;
                if(cluster_labels[id1]>=0)continue;

                for(int i=0;i<maxLabel+1;i++)votes[i] = 0;

                for(int i=-2;i<=2;i++){
                    int x2 = x1+i;
                    if((x2<0)||(x2>=square_len))continue;
                    for(int j=-2;j<=2;j++) {
                        if((i==0)&&(j==0))continue;
                        int y2 = y1+j;
                        if((y2<0)||(y2>=square_len))continue;
                        int gid2 = x2*square_len+y2;
                        int id2 = grid_asses[gid2];
                        if(id2>=num)continue;
                        if(cluster_labels[id2]<0)votes[maxLabel] += 1;
                        else votes[cluster_labels[id2]] += 1;
                    }
                }
                int maxlb = 0;
                for(int i=1;i<maxLabel;i++)if(votes[i]>votes[maxlb])maxlb=i;
                if(it>=maxit){
                    if(votes[maxlb]>0)cluster_labels[id1] = maxlb;
                    else flag = 1;
                }else {
                    if((votes[maxlb]>0)&&(votes[maxlb]>votes[maxLabel]*2))cluster_labels[id1] = maxlb;
                    else flag = 1;
                }
            }
        }
        if(flag==0)break;
    }

    delete[] fa;
    delete[] size;
    delete[] list;
    delete[] votes;
    return maxLabel;
}

#endif