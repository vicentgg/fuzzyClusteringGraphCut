#ifndef MESHSEGMENTATION_MESH_H
#define MESHSEGMENTATION_MESH_H

#include <vector>
#include <map>
#include <fstream>
#include <cmath>
#include <queue>
#include <set>
#include "Face.h"
#include "../Flow.h"
#include "../util/Util.h"
#include "Vertex.h"

#define EPSILON 1e-12
#define INF_FLOAT 1e12
#define FUZZY_REGION 0.25 //模糊区域概率范围
#define FUZZY -2
#define REPA -3
#define REPB -1

typedef std::pair<float, int> fipair;

class Mesh
{
public:
    std::vector<Face> faces;
    std::vector<Vertex> vertices;

    Mesh() = default;

    Mesh(const Mesh& an);


    //二维数组 distance[i][j]表示面片i到面片j的最短距离
    std::vector<std::vector<float>> distance;
    double avg_angle_dist;

    //专门用于记录面片隶属区域概率的数组
    std::vector<float> prob; 
    std::vector<int> partition;

    int num_faces;
    //-----------------------------------------------------------------------------------------------------------------
    /** Returns true if Faces a and b are neighbors */
    bool isNeighFace(const Face &a, const Face &b, std::vector<int> &common);
    //-----------------------------------------------------------------------------------------------------------------
    /* Computing distance and assign probabilities.
     * Reference: Hierarchical Mesh Decomposition using Fuzzy Clustering and Cuts (3.1 and 3.2)
     * The probability that a face belongs to a certain patch depends on its distance from other faces in this patch */
    //-----------------------------------------------------------------------------------------------------------------
    /** Compute angular distance (angle between face`s normals)*/
    double angleDist(const Face &a, const Face &b, float &angle);
    /** Compute geodesic distance (shortest path between face`s dual vertices on the dual graph)*/
    double geoDist(const Face &a, const Face &b, const std::vector<int> &common);
    /** Compute distances for each face and get neighbors */
    void getDual();
    /** Compute the shortest path for each face*/
    void compDist();
    /** Dijkstra algorithm */
    void dijkstra(int src, std::vector<float> &d);
    //-----------------------------------------------------------------------------------------------------------------
    /* Mesh Segmentation.
     * Reference: Hierarchical Mesh Decomposition using Fuzzy Clustering and Cuts (3.3 and 3.4) */
    /** Generating the final decomposition */
    void meshSeg(int depth, int id, std::vector<int> &vs, std::map<int, int> &part);
    /** Initializing before mesh segmentation */
    void solve();
    /** Generating a fuzzy decomposition */
    void fuzzyCluster(std::vector<int> &vs, std::map<int, int> &part, int &fuzzy,
                      int repA, int repB);
    /** Building Flow Net(Graph) to find a boundary between the components by applying a maximum flow algorithm */
    bool buildFlowGraph(std::vector<int> &vs, std::map<int, int> &part);
    //-----------------------------------------------------------------------------------------------------------------
    // Output
    void writeOff(const std::string& filename);
};

Mesh::Mesh(const Mesh &an)
{
    faces = an.faces;
    vertices = an.vertices;
}

/** Find |v| */
inline float vlen(const Point &v)
{
    float mag = v.x[0] * v.x[0] + v.x[1] * v.x[1] + v.x[2] * v.x[2];
    return sqrtf(mag);
}

//判断两个面片是否是相邻的
bool Mesh::isNeighFace(const Face &a, const Face &b, std::vector<int> &common)
{
    common.clear();
    Indices ia = a.indices, ib = b.indices;
    //统计三个顶点重复的个数
    int same = 0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (ia[i] == ib[j])
            {
                common.push_back(ia[i]);
                same++;
                break;
            }
        }
    }

    if (same == 2)
        return true;

    if (same == 3)
    {
        std::cout << "Same face !!!" << std::endl;
    }
    return false;
}

//计算两个面片的角度距离
double Mesh::angleDist(const Face &a, const Face &b, float &angle)
{
    float dot = util::dot(a.normal, b.normal);
    if (dot >= 1.0f)
        dot = 1.0f;
    else if (dot <= -1.0f)
        dot = -1.0f;
    angle = asin(dot);
    float angle_dist = 1 - cos(angle);
    //判断两个面片是否为凹边 当为凸边时 取0.2
    if (util::dot(b.normal - a.normal, a.center - b.center) < 1e-12)
    {
        angle_dist *= 0.2;
    }
    return angle_dist;
}

//计算测地距离 测地距离为相邻面片质心之间的距离
double Mesh::geoDist(const Face &a, const Face &b, const std::vector<int> &common)
{
    assert(common.size() == 2);
    Point p0 = vertices[common[0]].p, p1 = vertices[common[1]].p; 
    Point ce = p1 - p0, va = a.center - p0, vb = b.center - p0;
    float l = vlen(ce), la = vlen(va), lb = vlen(vb);

    double va_ce = util::dot(va, ce) / (la * l);
    va_ce = std::clamp(va_ce, -1., 1.); //防止出现越界
    double vb_ce = util::dot(vb, ce) / (lb * l);
    vb_ce = std::clamp(vb_ce, -1., 1.);

    double angle;
    angle = acos(va_ce) + acos(vb_ce);

    return la * la + lb * lb - 2 * la * lb * cos(angle);
}


/**
 * step1 计算原始网格的面片之间的距离（测地距离和角度距离），然后构建网格对偶图，计算并记录对偶图中所有顶点的最近距离。
 */
//1.1 计算每个面片与相邻面片之间的角度距离与测地距离 然后建立网格的对偶图 使用相邻面片的角度距离和测地距离计算对偶图中每个边缘的权重
void Mesh::getDual()
{
    int num_neigh = 0;
    num_faces = faces.size();
    double tot_angle_dist = 0.0f, tot_geo_dist = 0.0f;
    std::cout << num_faces << std::endl;
    for (int i = 0; i < num_faces; ++i)
    {
        for (int j = i + 1; j < faces.size(); ++j)
        {
            std::vector<int> common;
            if (isNeighFace(faces[i], faces[j], common))
            {
                num_neigh++;
                float angle;
                double angle_dist = angleDist(faces[i], faces[j], angle);
                double geo_dist = geoDist(faces[i], faces[j], common);
                tot_angle_dist += angle_dist;
                tot_geo_dist += geo_dist;
                assert(!isnan(tot_geo_dist));

                //将面片加入到对应的相邻面片数组中 角度距离 测地距离 二面角    
                faces[i].dedges.emplace_back(j, angle_dist, geo_dist, angle); 
                faces[j].dedges.emplace_back(i, angle_dist, geo_dist, angle);
            }
        }
    }
    avg_angle_dist = tot_angle_dist / num_neigh;
    double avg_geo_dist = tot_geo_dist / num_neigh;
    for (Face &f: faces)
    {
        for (DualEdge &de: f.dedges)
        {
            de.weight = 0.2 * de.angle_dist / avg_angle_dist +
                        0.8 * de.geo_dist / avg_geo_dist;
            std::cout << de << std::endl;

        }
    }
    std::cout << "Num neighbor faces " << num_neigh << std::endl;
    std::cout << "Avg Angle dist " << avg_angle_dist << std::endl;
    std::cout << "Avg Geo dist " << avg_geo_dist << std::endl;
}

//1.2 预处理阶段，计算对偶图中所有顶点之间的最近距离 计算每个面片到其他面片的最短路径
void Mesh::compDist()
{
    /**
     * 目前这块地方的性能损失较大
     * 还有优化的空间
     */
    for (int i = 0; i < num_faces; i++)
    {
        std::vector<float> tmp_dist(num_faces, INF_FLOAT); //建立最短路径数组 初始值设为最大 下标为对应的面片
        dijkstra(i, tmp_dist);//使用最短路径算法 计算当前面片到其他面片的最短距离
        distance.push_back(tmp_dist);
    }
    std::cout << "Finish computing shortest path" << std::endl;
}

//1.2.1 计算最短路径问题 这个地方可以进一步优化 迪斯捷特拉算法
void Mesh::dijkstra(int src, std::vector<float> &tmp_dist)
{
    //定义一个小顶堆
    std::priority_queue<fipair, std::vector<fipair>, std::greater<>> pq;
    pq.push(std::make_pair(0.0f, src));
    tmp_dist[src] = 0; //到本身面片的距离为0
    while (!pq.empty())
    {
        /**
         * 获取堆顶元素 并弹出
         * 堆顶为当前记录下src到面片u的最短路径
         * 遍历面片u的相邻面片 求出src到u相邻面片的最短路径
         * 并更新数据 加入到小顶堆中
         */
        int u = pq.top().second;
        pq.pop();

        for (DualEdge &de: faces[u].dedges)
        {
            int v = de.face;
            float weight = de.weight; //使用上一步计算两个相邻面片之间的角度距离和测地距离 得到的权重路径

            if (tmp_dist[v] > tmp_dist[u] + weight) //当已记录的距离 大于src-u的距离+u-v距离 则进行更新
            {
                tmp_dist[v] = tmp_dist[u] + weight; //更新为当前面到v面的路径长度 + 当前面的最短路径长度
                pq.push(std::make_pair(tmp_dist[v], v)); //更新小顶堆
            }
        }
    }
}


/**
 * step2 进行网格的分割
 * 选取对偶图中相距最远的两个面片作为种子
 * 使用模糊聚类更新种子，确定分割区域和模糊区域
 * 使用最小割进行模糊区域的划分
 * 对每个得到的分割区域重复上述步骤
 */
void Mesh::solve()
{
    partition.resize(num_faces); //记录每个面的最终分类结果
    std::map<int, int> part;  //记录每次分割时 每个面和对应的分类结果
    std::vector<int> vs;    
    vs.reserve(num_faces); //记录每次递归时 操作区域的面片下标
    for (int i = 0; i < num_faces; i++)
        vs.push_back(i);

    meshSeg(0, 0, vs, part);
}
// 2. 递归进行区域划分 基础划分为一分二 然后不断递归 
void Mesh::meshSeg(int depth, int id, std::vector<int> &vs,
                        std::map<int, int> &part)
{
    std::cout << "Depth and id: " << depth << " " << id << std::endl;//输出递归深度和标签值
    int fuzzy = 0; //记录模糊的区域面片数量
    int repA, repB; //种子面片的下标
    float max_dist = -1.0f;
    float max_angle = -100, min_angle = 100;

    std::set<int> cur_vs(vs.begin(), vs.end()); //记录面片的起始下标
    //遍历该区域内 找到最大的二面角和最小的二面角
    for (int v: vs)
    {
        for (DualEdge &de: faces[v].dedges)
        {
            if (cur_vs.find(de.face) == cur_vs.end()) continue;
            if (de.angle > max_angle) max_angle = de.angle;
            if (de.angle < min_angle) min_angle = de.angle;
        }
    }
    //2.1 找到该区域内距离最远的两个面片 当作两个分割区域的种子面片
    for (int i = 0; i < vs.size(); i++)
    {
        for (int j = 0; j < vs.size(); j++)
        {
            int m = vs[i], n = vs[j];
            if (distance[m][n] > max_dist)
            {
                max_dist = distance[m][n];
                if (m < n)
                {
                    repA = m;
                    repB = n;
                } else
                {
                    repA = n;
                    repB = m;
                }
            }
        }
    }
    //输出当前区域的阈值参数大小
    std::cout << "threshold " << distance[repA][repB] << " "
              << max_angle - min_angle << std::endl;
    //当最远距离和最大角度差 小于设定的阈值时 变停止分割 并给当前区域的所有面片打上分类标签          
    if (distance[repA][repB] < 30 || max_angle - min_angle < 1.1)
    {
        for (int v: vs)
        {
            partition[v] = id;
        }
        return;
    }

    //2.1 进行模糊聚类 将网格划分为三个区域 区域A 区域B 和模糊区域
    fuzzyCluster(vs, part, fuzzy, repA, repB);

    //2.3 模糊区域面片数量不为0时 使用最小割方法 寻找准确边界
    if (fuzzy != 0)
        //如果建立模糊区流图失败
        if (!buildFlowGraph(vs, part)) 
        {
            for (int v: vs)
            {
                //则直接对该网格区域进行标记 丢弃上一步的模糊分类结果
                partition[v] = id; 
            }
            return;
        };

    
    //当递归深度到达一定时 直接对区域的网格按照上述分割结果打标签 不再继续进行分割
    if (depth >= 5)
    {
        int pa = 0, pb = 0;
        for (int v: vs)
        {
            assert(part[v] != FUZZY);
            if (part[v] == REPA)
            {
                pa++;
                partition[v] = id * 2 + 1;
            } else
            {
                pb++;
                partition[v] = id * 2 + 2;
            }
        }
        std::cout << "Finish partition, pa and pb " << pa << " " << pb << std::endl;
        return;
    }

    // 聚合区域A和B
    std::vector<int> vs_a, vs_b;
    for (int v: vs)
    {
        assert(part[v] != FUZZY);//判断是否还存在模糊区域
        if (part[v] == REPA)
        {
            vs_a.push_back(v);
        } else
        {
            vs_b.push_back(v);
        }
        part.erase(v);
    }

    std::cout << "VS, VS A, VS B " << vs.size() << " " << vs_a.size() << " "
              << vs_b.size() << std::endl;

    meshSeg(depth + 1, id * 2 + 1, vs_a, part);

    meshSeg(depth + 1, id * 2 + 2, vs_b, part);
}

//2.1 在前面找到初始种子区域后 开始进行模糊聚类
void Mesh::fuzzyCluster(std::vector<int> &vs, std::map<int, int> &part, int &fuzzy, int repA, int repB)
{
    int nv = vs.size();
    //重置概率数组
    prob.reserve(nv);
    std::cout << "Representatives " << repA << " " << repB << std::endl;

    /**
     * 利用模糊聚类的方法迭代更新两个区域的种子面片
     * 1.计算每个面片分别属于两个区域的概率
     * 2.重新计算种子面片
     * 3.迭代上述步骤 直到种子面片不再发生变化
     * 迭代次数规定不超过十次 时间复杂度O(n^2)
     */
    for (int iter = 0; iter < 10; iter++)
    {
        //使用种子面片到每个面片的距离 计算每个面片属于属于区域A和B的概率
        for (int t = 0; t < nv; t++)
        {
            int v = vs[t];
            float a_fi = distance[v][repA], b_fi = distance[v][repB];
            prob[t] = a_fi / (a_fi + b_fi); //记录了当前面片属于区域B的概率
        }       

        //重新计算种子面片
        int lastA = repA, lastB = repB;
        float min_repB = INF_FLOAT, min_repA = INF_FLOAT; //记录统计的最小值
        /**
         * 统计当前面片i到每个面片之间的距离与两个区域概率的乘积之和
         * 当对应区域的乘积之和最小时 对应区域的种子面片更新为当前面片
         */
        for (int i = 0; i < nv; i++)
        {
            int v1 = vs[i];
            int tmp_repB = 0.0f, tmp_repA = 0.0f;
            for (int j = 0; j < nv; j++)
            {
                
                int v2 = vs[j];
                tmp_repB += prob[j] * distance[v1][v2];
                tmp_repA += (1.0 - prob[j]) * distance[v1][v2];
            }
            if (tmp_repB < min_repB)
            {
                min_repB = tmp_repB;
                repB = v1;
            }
            if (tmp_repA < min_repA)
            {
                min_repA = tmp_repA;
                repA = v1;
            }
        }
        std::cout << repA << " " << repB << std::endl;
        //当种子面片不再发生变化时 停止迭代
        if (lastA == repA && lastB == repB)
            break;
    }

    //当迭代收敛之后 将网格分为三个区域
    for (int i = 0; i < nv; i++)
    {
        int vertice = vs[i];

        if (prob[i] > 0.5 + FUZZY_REGION)
        {
            part[vertice] = REPB;
        } else if (prob[i] < 0.5 - FUZZY_REGION)
        {
            part[vertice] = REPA;
        } else
        {
            part[vertice] = FUZZY;
            fuzzy++; //记录模糊区域面片数量
        }
    }
    std::cout << "Finish fuzzy decomposition, found " << fuzzy
              << " fuzzy vertices" << std::endl;
}


// 2.2 使用最小割的方法寻找准确边界
// vs - 用于分解网格的所有顶点
// part - vs对应的分区
bool Mesh::buildFlowGraph(std::vector<int> &vs, std::map<int, int> &part)
{
    FlowNet flowNet(num_faces);
    // 原始顶点到流网络顶点
    std::map<int, int> ori2flow;
    std::set<int> v_ca, v_cb;
    std::set<int> cur_vs(vs.begin(), vs.end());

    int src = 0, dst;
    int vcnt = 1, fuzzy = 0;
    for (int v: vs)
    {
        if (part[v] == FUZZY)
        {
            fuzzy++;
            if (ori2flow.find(v) == ori2flow.end())
            {
                ori2flow[v] = vcnt++;
            }
            int m = 0;
            for (DualEdge &de: faces[v].dedges)
            {
                if (cur_vs.find(de.face) == cur_vs.end()) continue;
                int neibor = de.face, par = part[neibor];
                if (ori2flow.find(neibor) == ori2flow.end())
                {
                    ori2flow[neibor] = vcnt++;
                }
                bool add_edge = true;
                if (par == FUZZY)
                {
                    assert(v != neibor);
                    if (neibor < v)
                    {
                        add_edge = false;
                    }
                } else if (par == REPA)
                {
                    v_ca.insert(neibor);
                } else if (par == REPB)
                {
                    v_cb.insert(neibor);
                }
                if (add_edge)
                {
                    float cap = 1.0 / (1 + de.angle_dist / avg_angle_dist);
                    int from = ori2flow[v], to = ori2flow[neibor];
                    flowNet.addEdge(from, to, cap);
                }
            }
        }
    }
    dst = vcnt;
    std::cout << "Vertices all: " << vcnt + 1 << " fuzzy: " << fuzzy
              << " VCA: " << v_ca.size() << " VCB: " << v_cb.size() << std::endl;
    if (v_ca.size() == 0 || v_cb.size() == 0) return false;
    // Add edges for src and dst
    for (int vca: v_ca)
    {
        flowNet.addEdge(0, ori2flow[vca], INF_FLOAT);
    }
    for (int vcb: v_cb)
    {
        flowNet.addEdge(dst, ori2flow[vcb], INF_FLOAT);
    }
    flowNet.num_v = vcnt + 1;
    flowNet.EK(src, dst);
    // 分配流网络结果
    for (std::map<int, int>::const_iterator it = ori2flow.begin();
         it != ori2flow.end(); it++)
    {
        if (flowNet.visit[it->second])
        {
            partition[it->first] = REPA;
            part[it->first] = REPA;
        } else
        {
            partition[it->first] = REPB;
            part[it->first] = REPB;
        }
    }
    return true;
}

void Mesh::writeOff(const std::string &filename)
{
    std::ofstream file(filename);

    std::cout << "Writing in " << filename << std::endl;

    file << "OFF" << std::endl;
    file << vertices.size() << " " << faces.size() << " 0" << std::endl;

    for (auto &v : vertices)
        file << v.p.x[0] << " " << v.p.x[1] << " " << v.p.x[2] << std::endl;

    int i = 0;
    for (auto &f : faces)
    {
        file << "3 " << f.indices.x << " " << f.indices.y << " " << f.indices.z << " ";
        file << 60 * (partition[i] % 4 + 1) << " " << 80 * ((partition[i] + 1) % 3 + 1) << " " << 50 * ((partition[i] + 2) % 5 + 1) << std::endl;
        i++;
    }
}


#endif //MESHSEGMENTATION_GRAPHMESH_H

