#ifndef MESHSEGMENTATION_FLOW_H
#define MESHSEGMENTATION_FLOW_H

#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include "model/Point.h"

#define FLOW_EPSILON 1e-10
#define INF 1e12

struct Edge
{
    int to, next;
    int from;
    float cap;

    friend std::ostream &operator<<(std::ostream &os, const Edge &e)
    {
        return os << e.from << ", " << e.to << ", " << e.cap;
    }
};

/** Building FlowNet to solve Max Flow Problem */
class FlowNet
{
public:
    int num_v, edges_num;
    std::vector<int> head, last;
    std::vector<Edge> edges;
    std::vector<bool> visit;
    std::vector<float> flow;

    explicit FlowNet(int num_ver);

    void addEdge(int from, int to, float cap);

    bool bfs(int src, int dst);

    float EK(int src, int dst);

    void findCut(int src);
};

FlowNet::FlowNet(int num_ver) : num_v(num_ver), edges_num(0)
{
    head.resize(num_v);
    last.resize(num_v);
    flow.resize(num_v);
    visit.resize(num_v);
    for (int i = 0; i < num_v; i++)
    {
        head[i] = -1;
        last[i] = -1;
        visit[i] = false;
        flow[i] = 0.0f;
    }
    edges.resize(num_v * 2 + 10);
}

void FlowNet::addEdge(int from, int to, float cap)
{
    edges[edges_num].from = from;
    edges[edges_num].to = to;
    edges[edges_num].cap = cap;
    edges[edges_num].next = head[from];
    head[from] = edges_num++;

    edges[edges_num].from = to;
    edges[edges_num].to = from;
    edges[edges_num].cap = cap;
    edges[edges_num].next = head[to];
    head[to] = edges_num++;
    if (edges_num >= edges.size())
    {
        std::cout << "Edges exceeded, resize" << std::endl;
        edges.resize(edges.size() * 2);
    }
}

bool FlowNet::bfs(int s, int t)
{
    for (int i = 0; i < num_v; i++)
    {
        last[i] = -1;
    }
    std::queue<int> q;
    q.push(s);
    flow[s] = INF;
    while (!q.empty())
    {
        int p = q.front();
        q.pop();
        if (p == t)
            break;
        for (int eg = head[p]; eg != -1; eg = edges[eg].next)
        {
            int to = edges[eg].to;
            float vol = edges[eg].cap;
            if (vol > FLOW_EPSILON && last[to] == -1)
            {
                last[to] = eg;
                flow[to] = fmin(flow[p], vol);
                q.push(to);
            }
        }
    }
    return last[t] != -1;
}

float FlowNet::EK(int s, int t)
{
    float maxflow = 0;
    while (bfs(s, t))
    {
        maxflow += flow[t];
        for (int i = t; i != s; i = edges[last[i] ^ 1].to)
        {
            edges[last[i]].cap -= flow[t];
            edges[last[i] ^ 1].cap += flow[t];
        }
    }

    findCut(0);
    return maxflow;
}

void FlowNet::findCut(int src)
{
    visit[src] = true;
    for (int i = head[src]; i != -1; i = edges[i].next)
    {
        int to = edges[i].to;
        if (edges[i].cap > FLOW_EPSILON && !visit[to])
        {
            findCut(to);
        }
    }
}



#endif //MESHSEGMENTATION_FLOW_H
