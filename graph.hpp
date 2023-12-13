#ifndef Graph_hpp
#define Graph_hpp

#include <string.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// #define _LINUX_

#ifdef _LINUX_
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
// #include <sys/resource.h>
#endif

using namespace std;

class Graph {
    unsigned int n_{};      
    unsigned int m_{};
    unsigned int m_max_{};
    unsigned int sup_max_{};
    long long idx_size_;

    FILE* log_f_;

    vector<long> t_new_to_old_;
    vector<int> edges_idx_;                          // 不同时间点对应的下标
    vector<pair<int, int>> edges_t_;                 // 边集按时间顺序存放
    vector<unordered_map<int, vector<int>>> nbr_t_;  // u < v 单向边，统计每条边及其出现时间

    vector<pair<int, int>> edges_;  // 单向边，每个边只存一次，u < v

    unordered_map<int, int>* nbr_;  // 单向边，value = eid,映射edges

    int* nbr_cnt_;                  // 单向边，计数，每个边只存一次，u < v

    vector<int>* cn_;                          // 共同邻居
    int* sup_;                                 // 支持度
    int* truss_;                               // e的最大truss (u<v,只存一条边)
    int* td_;                                  // truss_degrss
    vector<vector<pair<int, int>>>* truss_t_;  // 索引

    int* tt_cnt_;  // truss time degree

    bool* visit_;

    void init_nbr_cnt();
    void compute_common_neighbors();
    void truss_decomposition();
    void init_truss_deg();

    void compute_cn_time_bl(const int& t_s);
    void decremental_cn_bl(const int& t_s);

    void init_truss_time();
    void init_tt_cnt(const int& k);

    void print_graph_size();
    void print_idx_size();

   public:
    unsigned int t_{};
    int k_max_{};
    Graph();
    ~Graph();
    void load(const string& path);

    void index_baseline();
    void index();

    void write_idx_txt(const string& path);
    void write_idx(const string& path);
    void load_idx(const string& path);
    void init_log(const string& log_path);

    bool query(int u, int v, int t_s, int t_e, int k);
    int query_all(int t_s, int t_e, int k);
    int online_query(int t_s, int t_e, int k, int& snapshot_m);
};

bool cmp1(const pair<int, int>& a, const pair<int, int>& b);
bool cmp2(const pair<int, int>& a, const pair<int, int>& b);

#endif /* Graph_hpp */
