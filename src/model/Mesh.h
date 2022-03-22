#ifndef MESHSEGMENTATION_MESH_H
#define MESHSEGMENTATION_MESH_H

#include <vector>
#include <time.h>
#include <cstdlib>
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
#define FUZZY_REGION 0 //模糊区域概率范围
#define FUZZY -2
#define REPA -3
#define REPB -1
#define PI acos(-1)

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
    std::vector<int> visited_par;

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
    void solve_feiliu(std::vector<int>& vs, int slog);

    /** Generating a fuzzy decomposition */
    void fuzzyCluster(std::vector<int> &vs, std::map<int, int> &part, int &fuzzy,
                      int repA, int repB);
    /** Building Flow Net(Graph) to find a boundary between the components by applying a maximum flow algorithm */
    bool buildFlowGraph(std::vector<int> &vs, std::map<int, int> &part);
    //-----------------------------------------------------------------------------------------------------------------
    // Output
    void writeOff(const std::string& filename);

    int reType(double rand_sig);

    //计算形状指数
    double ShapeIndex(); 
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

    double va_vb = util::dot(va, vb) / (la * lb);
    va_vb = std::clamp(va_vb, -1., 1.);


    // double va_ce = util::dot(va, ce) / (la * l); //求解两个向量之间的夹角
    // va_ce = std::clamp(va_ce, -1., 1.); //防止出现越界
    // double vb_ce = util::dot(vb, ce) / (lb * l);
    // vb_ce = std::clamp(vb_ce, -1., 1.);

    double angle;
    // angle = acos(va_ce) + acos(vb_ce);

    angle = acos(va_vb);

    return la * la + lb * lb - 2 * la * lb * cos(angle);
}


// 获取当前系统时间
inline int64_t getCurrentLocalTimeStamp()
{
    std::chrono::time_point<std::chrono::system_clock, std::chrono::milliseconds> tp = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
    auto tmp = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch());
    return tmp.count();
}


//计算顶点的形状指数 
//分别计算顶点的高斯曲率与平均曲率
double Mesh::ShapeIndex() {
    int64_t start_time = getCurrentLocalTimeStamp();

    int num_verts = vertices.size();

    for(int i = 0; i < num_verts; i ++) {
        double sum_angle = 0;
        double sum_s = 0;
        Point vert_n(0,0,0); //顶点法向量

        for(int j = 0, len = vertices[i].neighbor_face.size(); j < len; j ++) {
            // std::cout << vertices[i].neighbor_face[j] << ' ';
            //获取顶点的一阶邻域面片

            //1.计算每个面片的面积 与内角
            const Face &face = faces[vertices[i].neighbor_face[j]];
            Indices a = face.indices;

            int v;
            for(v = 0; v < 3; v ++) {
                if(a[v] == i) break;
            }
            int v1,v2;
            for(v1 = 0; v1 < 3; v1 ++) {
                if(v1 != v) break;
            }
            for(v2 = 0; v2 < 3; v2 ++) {
                if(v2 != v && v2 != v1) break;
            }
            
            Point p = vertices[a[v]].p, p1 = vertices[a[v1]].p, p2 = vertices[a[v2]].p;
            Point n1 = p1 - p, n2 = p2 - p;
            util::normalizeV(n1);
            util::normalizeV(n2);
            float l1 = vlen(n1), l2 = vlen(n2);


            double l1_l2 = util::dot(n1, n2) / (l1*l2);
            l1_l2 = std::clamp(l1_l2, -1., 1.);
            double angle = acos(l1_l2);
            vertices[i].neighbor_angle.emplace(vertices[i].neighbor_face[j], angle);

            vert_n = vert_n + face.normal;
            sum_angle += angle;

            if(angle > PI / 2) angle = PI - angle;
            double s = l1*l2*sin(angle)/2;

            sum_s += s;
        }

        std::cout << "sum_angle " << sum_angle << std::endl;
        vertices[i].gaus = 3 * (2 * PI - sum_angle) / sum_s; //高斯曲率
        std::cout << "gaus " << vertices[i].gaus << std::endl;
        vert_n = vert_n / vertices[i].neighbor_face.size();
        vertices[i].normal = vert_n;  //顶点法向量
        vertices[i].s = sum_s;


    }

    //计算平均曲率
    for(int i = 0; i < num_verts; i ++) {

        double avrage = 0;

        for(int j = 0, len = vertices[i].neighbor_face.size(); j < len; j ++) {

            const Face &face = faces[vertices[i].neighbor_face[j]];
            Indices a = face.indices;

            int v; //找到中心点 a[v]
            for(v = 0; v < 3; v ++) {
                if(a[v] == i) break;
            }


            for(int k = 0, len = vertices[i].neighbor_face.size(); k < len; k ++) {
                if(k != j) {
                    const Face &face = faces[vertices[i].neighbor_face[k]];
                    Indices b = face.indices;
                    std::vector<int> common;

                    if (isNeighFace(faces[j], faces[k], common)) {

                        int pair;
                        for(auto com : common) {
                            if(com != a[v]) pair = com;
                        }

                        int l,r;
                        for(int m = 0; m < 3; m ++) {
                            if(a[m] != a[v] && a[m] != pair) {
                                l = a[m];
                                break;
                            }
                        }
                        for(int m = 0; m < 3; m ++) {
                            if(b[m] != a[v] && b[m] != pair) {
                                r = b[m];
                                break;
                            }
                        }
                        // std::cout << "l " << l << " r " << r << std::endl;

                        Point n1 = vertices[pair].p - vertices[a[v]].p;
                        util::normalizeV(n1);
                        Point normal = vertices[a[v]].normal;
                        util::normalizeV(normal);
                        double vert_n_n1 = util::dot(normal, n1);
                        // std::cout << "vert_n_n1: " << vert_n_n1 << std::endl;
                        double l_angle = vertices[l].neighbor_angle[vertices[i].neighbor_face[j]];

                        double r_angle = vertices[r].neighbor_angle[vertices[i].neighbor_face[k]];

                        double cot_ab = (1/tan(l_angle) + 1/tan(r_angle)) / 2;
                        // std::cout << "cot_ab: " << cot_ab << std::endl;
                        avrage = avrage + (cot_ab*vert_n_n1);
                    }
                }
            }
        }
        // std::cout << "vertices[i].avrage: " << avrage << std::endl;
        // std::cout << "vertices[i].s: " << vertices[i].s << std::endl;
        vertices[i].avrage = avrage * 3 / (4 * vertices[i].s);
        
    }


//计算顶点的形状指数
    for(int i = 0; i < num_verts; i ++) {
        
        double t = 0.5 - atan(vertices[i].avrage / sqrt(abs(vertices[i].avrage * vertices[i].avrage - vertices[i].gaus))) / PI;
        // std::cout << " avrage:" << vertices[i].avrage << " gaus:" << vertices[i].gaus << " type:" << t << std::endl;
        vertices[i].index = t;
    }

    int64_t end_time = getCurrentLocalTimeStamp();
    // std::cout << "形状信息统计：" << end_time - start_time << std::endl; 

    return 0;
}

/**
 * step1 计算原始网格的面片之间的距离（测地距离和角度距离），然后构建网格对偶图，计算并记录对偶图中所有顶点的最近距离。
 */
//1.1 计算每个面片与相邻面片之间的角度距离与测地距离 然后建立网格的对偶图 使用相邻面片的角度距离和测地距离计算对偶图中每个边缘的权重
void Mesh::getDual()
{
    int64_t start_time = getCurrentLocalTimeStamp();

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

                // assert(geo_dist < INF_FLOAT);

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
            // std::cout << de << std::endl;

        }
    }
    std::cout << "Num neighbor faces " << num_neigh << std::endl;
    std::cout << "Avg Angle dist " << avg_angle_dist << std::endl;
    std::cout << "Avg Geo dist " << avg_geo_dist << std::endl;


    ShapeIndex();

    int64_t end_time = getCurrentLocalTimeStamp();

    std::cout << "包含曲率的特征提取耗时：" << end_time - start_time << "ms" << std::endl;
}

//1.2 预处理阶段，计算对偶图中所有顶点之间的最近距离 计算每个面片到其他面片的最短路径
void Mesh::compDist()
{
    /**
     * 目前这块地方的性能损失较大
     * 还有优化的空间
     */
    int64_t start_time = getCurrentLocalTimeStamp();
    for (int i = 0; i < num_faces; i++)
    {
        std::vector<float> tmp_dist(num_faces, INF_FLOAT); //建立最短路径数组 初始值设为最大 下标为对应的面片
        dijkstra(i, tmp_dist);//使用最短路径算法 计算当前面片到其他面片的最短距离
        distance.push_back(tmp_dist);
    }
    std::cout << "Finish computing shortest path" << std::endl;
    int64_t end_time = getCurrentLocalTimeStamp();

    std::cout << "对偶图构建耗时：" << end_time - start_time << "ms" << std::endl;

    float max_float = 0;
    for(int i = 0; i < num_faces; i ++) {
        for(int j = 0; j < num_faces; j ++) {
            if(distance[i][j] > 1e+11) {
                max_float = 1;
                break;
            }
        }
        if(max_float == 1) break;
    }

    if(max_float == 1) {
        std::cout << "非流形网格" << std::endl;
        std::vector<int> visited(num_faces, -1);

        for(int i = 0; i < num_faces; i ++) {
            if(visited[i] == -1) {
                for(int j = 0; j < num_faces; j ++) {
                    if(distance[i][j] <= 1e+11) {
                        visited[j] = i;
                    }   
                }
            } 
        }

        visited_par = visited;

        std::unordered_map<int, std::vector<int>> v_map;
        for(int i = 0; i < visited.size(); i ++) {
            if(v_map.find(visited[i]) == v_map.end()) {
                std::vector<int> num;
            
                for(int j = i; j < visited.size(); j ++) {
                    if(visited[j] == visited[i]) num.push_back(j);
                }

                v_map[visited[i]] = num;
            }
        }
        

        for(auto &map_visited : v_map) {
            std::vector<int> vs = map_visited.second;
            int slog = map_visited.first;
            solve_feiliu(vs,slog);
        }

        std::cout <<"map size: " << v_map.size() << std::endl;
    }
    else {
        std::cout << "流形网格" << std::endl;
        solve();
    }
    

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

//适用于非流形
void Mesh::solve_feiliu(std::vector<int>& vs, int slog)
{
    partition.resize(num_faces); //记录每个面的最终分类结果
    std::map<int, int> part;  //记录每次分割时 每个面和对应的分类结果
    meshSeg(0, slog, vs, part);
}

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
    std::cout << "分割深度和记号： " << depth << " " << id << std::endl;//输出递归深度和标签值
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
            if (distance[m][n] > max_dist && distance[m][n] < INF_FLOAT)
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
    std::cout << "面片最远距离 角度最大差值 " << distance[repA][repB] << " "
              << max_angle - min_angle << std::endl;
    //当最远距离和最大角度差 小于设定的阈值时 变停止分割 并给当前区域的所有面片打上分类标签          
    if (distance[repA][repB] < 30 || max_angle - min_angle < 1.1)
    {
        std::cout << "低于阈值停止分割" << std::endl;
        for (int v: vs)
        {
            partition[v] = id;
        }
        return;
    }
    
    int64_t start_time = getCurrentLocalTimeStamp();

    //2.1 进行模糊聚类 将网格划分为三个区域 区域A 区域B 和模糊区域
    fuzzyCluster(vs, part, fuzzy, repA, repB);

    int64_t end_time = getCurrentLocalTimeStamp();
    std::cout << "一次聚类耗时：" << end_time - start_time << "ms" << std::endl;

    // //2.3 模糊区域面片数量不为0时 使用最小割方法 寻找准确边界
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
    if (depth >= 3)
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
            // else {
            //     partition[v] = id * 2 + 3;
            // }
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

    std::cout << "全部面, A面, B面 " << vs.size() << " " << vs_a.size() << " "
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
    std::cout << "初始种子面片： " << repA << " " << repB << std::endl;

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
            // std::cout << "每个面片属于b的概率： " << prob[t] << std::endl;
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
                continue;
            }
            if (tmp_repA < min_repA)
            {
                min_repA = tmp_repA;
                repA = v1;
            }
        }

        std::cout << "新的种子点 " << repA << " " << repB << std::endl;
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
    std::cout << "完成模糊聚类, 发现 " << fuzzy
              << " 个模糊点" << std::endl;
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
    std::cout << "Vertices all: " << vcnt + 1 << " 模糊区域: " << fuzzy
              << " A区域: " << v_ca.size() << " B区域: " << v_cb.size() << std::endl;
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
        int t = (reType(vertices[f.indices[0]].index)+reType(vertices[f.indices[1]].index)+reType(vertices[f.indices[2]].index))/3;
        file << 10 * (t + 1) << " " << 15 * (t + 1) << " " << 20 * (t + 2) << std::endl;

        // file << 60 * (partition[i] % 4 + 1) << " " << 80 * ((partition[i] + 1) % 3 + 1) << " " << 50 * ((partition[i] + 2) % 5 + 1) << std::endl;
        i++;
    }
}

int Mesh::reType(double rand_sig) {
    int type = 0;

    if(rand_sig <= 0.0625) {
        type = 1;
    } else if(rand_sig <= 0.1875){
        type = 2;
    } else if(rand_sig <= 0.3125){
        type = 3;
    } else if(rand_sig <= 0.4375){
        type = 4;
    } else if(rand_sig <= 0.5625){
        type = 5;
    } else if(rand_sig <= 0.6875){
        type = 6;
    } else if(rand_sig <= 0.8125){
        type = 7;
    } else if(rand_sig <= 0.9375){
        type = 8;
    } else if(rand_sig <= 1){
        type = 9;
    } else type = 0;

    return type;
}


#endif //MESHSEGMENTATION_GRAPHMESH_H

