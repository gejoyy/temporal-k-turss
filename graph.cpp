#include "graph.hpp"

Graph::Graph() {
    idx_size_ = 0;
    log_f_ = nullptr;
    nbr_ = nullptr;
    nbr_cnt_ = nullptr;
    cn_ = nullptr;
    sup_ = nullptr;
    truss_ = nullptr;
    td_ = nullptr;
    truss_t_ = nullptr;
    tt_cnt_ = nullptr;
    visit_ = nullptr;
}

Graph::~Graph() {
    if (nbr_ != nullptr) {
        delete[] nbr_;
        nbr_ = nullptr;
    }
    if (nbr_cnt_ != nullptr) {
        delete[] nbr_cnt_;
        nbr_cnt_ = nullptr;
    }
    if (cn_ != nullptr) {
        delete[] cn_;
        cn_ = nullptr;
    }
    if (sup_ != nullptr) {
        delete[] sup_;
        sup_ = nullptr;
    }
    if (truss_ != nullptr) {
        delete[] truss_;
        truss_ = nullptr;
    }
    if (td_ != nullptr) {
        delete[] td_;
        td_ = nullptr;
    }
    if (truss_t_ != nullptr) {
        delete[] truss_t_;
        truss_t_ = nullptr;
    }
    // if (tt_cnt_ != nullptr) { // 复用td_ 地址，不重复释放空间
    //     delete[] tt_cnt_;
    //     tt_cnt_ = nullptr;
    // }
    if (visit_ != nullptr) {
        delete[] visit_;
        visit_ = nullptr;
    }
    if (log_f_ != nullptr) {
        fclose(log_f_);
        log_f_ = nullptr;
    }
}

void Graph::load(const string& path) {
    if (log_f_ != nullptr) fprintf(log_f_, "Graph path: %s\n", path.c_str());
    printf("Graph path: %s\n", path.c_str());

    ifstream ifs(path);
    if (!ifs.is_open()) {
        cerr << "open file failed!" << endl;
        exit(-1);
    }

    n_ = 0;
    m_max_ = 0;
    t_ = 0;

    int u, v;
    long ts, pre_ts = -1;
    while (!ifs.eof()) {
        ifs >> u >> v >> ts;
        if (ifs.fail()) break;  // 防止多读一个空行

        if (u == v) continue;

        if (u > v) swap(u, v);  // 默认 u < v
        edges_t_.emplace_back(make_pair(u, v));

        // adjust size of neighbor list if necessary.
        if (v + 1 > nbr_t_.size()) {
            nbr_t_.resize(v + 1);
        }

        if (ts != pre_ts) {
            // 一起添加的 两个下标一致
            t_new_to_old_.emplace_back(ts);
            edges_idx_.emplace_back(edges_t_.size() - 1);
            pre_ts = ts;
        }

        int format_t = t_new_to_old_.size() - 1;

        nbr_t_[u][v].emplace_back(format_t);  // u < v

        ++m_max_;  // 统计使用的时间边，不包括u=v
    }
    ifs.close();

    n_ = nbr_t_.size();
    t_ = t_new_to_old_.size();
    edges_idx_.emplace_back(edges_t_.size());

    edges_ = edges_t_;
    sort(edges_.begin(), edges_.end());
    edges_.erase(unique(edges_.begin(), edges_.end()), edges_.end());
    m_ = edges_.size();

    if (nbr_ == nullptr) nbr_ = new unordered_map<int, int>[n_];
    for (int i = 0; i < m_; ++i) {
        nbr_[edges_[i].first].insert(make_pair(edges_[i].second, i));
        nbr_[edges_[i].second].insert(make_pair(edges_[i].first, i));
    }

    init_nbr_cnt();
    // init_edges_time();
    printf("n = %d, m_max = %d, t = %d, effective_m = %d.\n", n_, m_max_, t_, m_);

    if (log_f_ != nullptr) {
        fprintf(log_f_, "n = %d, m_max = %d, t = %d, effective_m = %d.\n", n_, m_max_, t_, m_);
    }

    print_graph_size();
}

// 统计边出现的次数
void Graph::init_nbr_cnt() {
    if (nbr_cnt_ == nullptr) nbr_cnt_ = new int[m_]{};  // 初始化为 0
    for (int i = 0; i < m_; ++i) {
        nbr_cnt_[i] = nbr_t_[edges_[i].first][edges_[i].second].size();
    }
}

// 计算边的共同邻居
void Graph::compute_common_neighbors() {
    if (cn_ == nullptr) cn_ = new vector<int>[m_] {};
    if (sup_ == nullptr) sup_ = new int[m_]{};  // 初始化为 0
    for (int eid = 0; eid < m_; ++eid) {
        int u = edges_[eid].first;
        int v = edges_[eid].second;
        if (nbr_[u].size() > nbr_[v].size()) swap(u, v);
        for (const auto& p : nbr_[u]) {
            if (nbr_[v].find(p.first) != nbr_[v].end()) {
                cn_[eid].emplace_back(p.first);
            }
        }
        sort(cn_[eid].begin(), cn_[eid].end());
        sup_[eid] = cn_[eid].size();
        if (sup_[eid] > sup_max_) sup_max_ = sup_[eid];
    }
}

// truss分解
void Graph::truss_decomposition() {
    if (truss_ == nullptr) truss_ = new int[m_];
    vector<int>* bin = new vector<int>[sup_max_ + 1];
    memset(visit_, false, sizeof(bool) * m_);
    memset(truss_, 2, sizeof(int) * m_);

    for (int i = 0; i < m_; ++i) {
        bin[sup_[i]].emplace_back(i);
    }

    k_max_ = 0;
    for (int i = 0; i <= sup_max_; ++i) {
        for (int j = 0; j < bin[i].size(); ++j) {
            int eid = bin[i][j];
            if (visit_[eid]) continue;
            visit_[eid] = true;
            truss_[eid] = i + 2;
            int u = edges_[eid].first;
            int v = edges_[eid].second;
            k_max_ = max(k_max_, i + 2);
            for (const auto& w : cn_[eid]) {
                int e1 = nbr_[u][w];
                int e2 = nbr_[v][w];
                if (!visit_[e1] && !visit_[e2]) {
                    bin[--sup_[e1]].emplace_back(e1);
                    bin[--sup_[e2]].emplace_back(e2);
                }
            }
        }
    }

    // 恢复 support
    for (int i = 0; i < m_; ++i) {
        sup_[i] = cn_[i].size();
    }

    memset(visit_, false, sizeof(bool) * m_);
    delete[] bin;
}

// 计算 truss degree
void Graph::init_truss_deg() {
    if (td_ == nullptr) td_ = new int[m_]{};  // 初始化为 0

    for (int i = 0; i < m_; ++i) {
        int u = edges_[i].first;
        int v = edges_[i].second;
        int trussness = truss_[i];
        for (const auto& w : cn_[i]) {
            int e1 = nbr_[u][w];
            int e2 = nbr_[v][w];
            if (trussness <= truss_[e1] && trussness <= truss_[e2]) {
                ++td_[i];
            }
        }
        if (td_[i] < trussness - 2) {
            printf("!!!Wrong: u=%d, v=%d, truss=%d, truss_deg_=%d\n", u, v, trussness, td_[i]);
            exit(-1);
        }
    }
}

void Graph::compute_cn_time_bl(const int& t_s) {
    queue<int> q;
    int* cnt = new int[k_max_ + 1];
    for (int t_e = t_ - 1; t_e >= t_s; --t_e) {  // 从t-max 减到 t_s
        for (int i = edges_idx_[t_e]; i < edges_idx_[t_e + 1]; ++i) {
            int u = edges_t_[i].first;
            int v = edges_t_[i].second;
            int eid = nbr_[u][v];

            if (--nbr_cnt_[eid] != 0) continue;

            int trussness = truss_[eid];
            for (int j = 0, sz = sup_[eid]; j < sz; ++j) {
                int w = cn_[eid][j];
                int e1 = nbr_[u][w];
                int e2 = nbr_[v][w];

                // process (u,w)
                iter_swap(find(cn_[e1].begin(), cn_[e1].begin() + sup_[e1], v), cn_[e1].begin() + sup_[e1] - 1);
                --sup_[e1];
                if (truss_[e1] <= min(trussness, truss_[e2])) {
                    --td_[e1];
                    if (td_[e1] < truss_[e1] - 2 && !visit_[e1]) {
                        q.push(e1);
                        visit_[e1] = true;
                    }
                }

                // process (v,w)
                iter_swap(find(cn_[e2].begin(), cn_[e2].begin() + sup_[e2], u), cn_[e2].begin() + sup_[e2] - 1);
                --sup_[e2];
                if (truss_[e2] <= min(trussness, truss_[e1])) {
                    --td_[e2];
                    if (td_[e2] < truss_[e2] - 2 && !visit_[e2]) {
                        q.push(e2);
                        visit_[e2] = true;
                    }
                }
            }
            if (!visit_[eid]) {
                q.push(eid);
                visit_[eid] = true;
            }
            sup_[eid] = 0;
            // cn_[eid].clear();
            // td_[eid] = 0;
        }

        while (!q.empty()) {
            int eid = q.front();
            q.pop();
            visit_[eid] = false;
            int u = edges_[eid].first;
            int v = edges_[eid].second;

            int ot = truss_[eid];
            memset(cnt, 0, sizeof(int) * (ot + 1));

            // 重新计算 truss(u,v) & td_[u][v]
            for (int i = 0, sz = sup_[eid]; i < sz; ++i) {
                int w = cn_[eid][i];
                int e1 = nbr_[u][w];
                int e2 = nbr_[v][w];
                ++cnt[min({ot, truss_[e1], truss_[e2]})];
            }
            int td = 0;
            for (int k = ot; k >= 2; --k) {
                td += cnt[k];
                if (td >= k - 2) {
                    truss_[eid] = k;
                    break;
                }
            }
            td_[eid] = td;
            int nt = truss_[eid];

            // 添加受影响的边进入队列
            for (int i = 0, sz = sup_[eid]; i < sz; ++i) {
                int w = cn_[eid][i];
                int e1 = nbr_[u][w];
                int e2 = nbr_[v][w];

                // process (u,w)
                if (nt < truss_[e1] && truss_[e1] <= min(ot, truss_[e2])) {
                    --td_[e1];
                    if (td_[e1] < truss_[e1] - 2 && !visit_[e1]) {
                        q.push(e1);
                        visit_[e1] = true;
                    }
                }

                // process (v,w)
                if (nt < truss_[e2] && truss_[e2] <= min(ot, truss_[e1])) {
                    --td_[e2];
                    if (td_[e2] < truss_[e2] - 2 && !visit_[e2]) {
                        q.push(e2);
                        visit_[e2] = true;
                    }
                }
            }

            // 更新索引
            for (int k = ot; k > nt; --k) {  // if nt = 2, k_min = 3
                if (t_s == 0 || truss_t_[eid][k - 2].empty() || truss_t_[eid][k - 2].back().second < t_e) {
                    truss_t_[eid][k - 2].emplace_back(make_pair(t_s, t_e));
                }
            }
        }
    }
    delete[] cnt;
}

//  删除 t 时间的边
void Graph::decremental_cn_bl(const int& t_s) {
    queue<int> q;
    int* cnt = new int[k_max_ + 1];
    for (int i = edges_idx_[t_s]; i < edges_idx_[t_s + 1]; ++i) {
        int u = edges_t_[i].first;
        int v = edges_t_[i].second;
        int eid = nbr_[u][v];

        if (--nbr_cnt_[eid] != 0) continue;

        int trussness = truss_[eid];
        for (int w : cn_[eid]) {
            int e1 = nbr_[u][w];
            int e2 = nbr_[v][w];

            // process (u,w)
            cn_[e1].erase(find(cn_[e1].begin(), cn_[e1].end(), v));  // 在cn(u,w) 删除v
            --sup_[e1];
            if (truss_[e1] <= min(trussness, truss_[e2])) {
                --td_[e1];
                if (td_[e1] < truss_[e1] - 2 && !visit_[e1]) {
                    q.push(e1);
                    visit_[e1] = true;
                }
            }

            // process (v,w)
            cn_[e2].erase(find(cn_[e2].begin(), cn_[e2].end(), u));  // 在cn(v,w) 删除u
            --sup_[e2];
            if (truss_[e2] <= min(trussness, truss_[e1])) {
                --td_[e2];
                if (td_[e2] < truss_[e2] - 2 && !visit_[e2]) {
                    q.push(e2);
                    visit_[e2] = true;
                }
            }
        }
        if (!visit_[eid]) {
            q.push(eid);
            visit_[eid] = true;
        }
        cn_[eid].clear();
        sup_[eid] = 0;
        // td_[eid] = 0;
    }

    while (!q.empty()) {
        int eid = q.front();
        q.pop();
        visit_[eid] = false;
        int u = edges_[eid].first;
        int v = edges_[eid].second;

        int ot = truss_[eid];
        memset(cnt, 0, sizeof(int) * (ot + 1));

        // 重新计算 truss(u,v) & td_(u,v)
        for (const int w : cn_[eid]) {
            int e1 = nbr_[u][w];
            int e2 = nbr_[v][w];
            ++cnt[min({ot, truss_[e1], truss_[e2]})];
        }
        int td = 0;
        for (int k = ot; k >= 2; --k) {
            td += cnt[k];
            if (td >= k - 2) {
                truss_[eid] = k;
                break;
            }
        }
        td_[eid] = td;
        int nt = truss_[eid];

        // 添加受影响的边进入队列
        for (const int w : cn_[eid]) {
            int e1 = nbr_[u][w];
            int e2 = nbr_[v][w];

            // process (u,w)
            if (nt < truss_[e1] && truss_[e1] <= min(ot, truss_[e2])) {
                --td_[e1];
                if (td_[e1] < truss_[e1] - 2 && !visit_[e1]) {
                    q.push(e1);
                    visit_[e1] = true;
                }
            }
            // process (v,w)
            if (nt < truss_[e2] && truss_[e2] <= min(ot, truss_[e1])) {
                --td_[e2];
                if (td_[e2] < truss_[e2] - 2 && !visit_[e2]) {
                    q.push(e2);
                    visit_[e2] = true;
                }
            }
        }

        // 更新索引
        for (int k = ot; k > nt; --k) {
            truss_t_[eid][k - 2].emplace_back(make_pair(t_s + 1, t_));
        }
    }
    delete[] cnt;
}

void Graph::index_baseline() {
#ifdef _LINUX_
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif
    if (visit_ == nullptr) visit_ = new bool[m_]{};
    if (truss_t_ == nullptr) truss_t_ = new vector<vector<pair<int, int>>>[m_];

    compute_common_neighbors();
    truss_decomposition();

    // 提前申请空间
    for (int i = 0; i < m_; ++i) {
        truss_t_[i].resize(truss_[i] - 1);  // "-2+1" truss >= 2 索引从下标为1的位置开始
    }

    init_truss_deg();

    // 备份数据 sup_, nbr_cnt_, truss_, td_
    auto sup_copy = new int[m_];
    auto nbr_cnt_copy = new int[m_];
    auto truss_copy = new int[m_];
    auto td_copy = new int[m_];

    for (int t_s = 0; t_s < t_; ++t_s) {
        if (t_s % 100 == 0) printf("t = %d.\n", t_s);

        // 备份数据
        memcpy(sup_copy, sup_, sizeof(int) * m_);
        memcpy(nbr_cnt_copy, nbr_cnt_, sizeof(int) * m_);
        memcpy(truss_copy, truss_, sizeof(int) * m_);
        memcpy(td_copy, td_, sizeof(int) * m_);

        compute_cn_time_bl(t_s);
        if (t_s == t_ - 1) break;

        // 恢复数据
        memcpy(sup_, sup_copy, sizeof(int) * m_);
        memcpy(nbr_cnt_, nbr_cnt_copy, sizeof(int) * m_);
        memcpy(truss_, truss_copy, sizeof(int) * m_);
        memcpy(td_, td_copy, sizeof(int) * m_);

        decremental_cn_bl(t_s);
    }

#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec) / 1000;
    printf("Running time (baseline): %lld s, %lld mins\n", t_msec / 1000, t_msec / 1000 / 60);
    if (log_f_ != nullptr) fprintf(log_f_, "Indexing time (Baseline): %lld s, %lld mins\n", t_msec / 1000, t_msec / 1000 / 60);
#else
    clock_t end = clock();
    printf("Running time (baseline): %.2f s, %.2f min\n", (double)(end - start) / CLOCKS_PER_SEC,
           (double)(end - start) / CLOCKS_PER_SEC / 60);
    if (log_f_ != nullptr)
        fprintf(log_f_, "Running time (baseline): %.2f s, %.2f min\n", (double)(end - start) / CLOCKS_PER_SEC,
                (double)(end - start) / CLOCKS_PER_SEC / 60);
#endif

    print_idx_size();
    delete[] sup_copy;
    delete[] nbr_cnt_copy;
    delete[] truss_copy;
    delete[] td_copy;
}

void Graph::init_truss_time() {
    // 备份数据
    auto nbr_cnt_copy = new int[m_];
    auto truss_copy = new int[m_];
    memcpy(nbr_cnt_copy, nbr_cnt_, sizeof(int) * m_);
    memcpy(truss_copy, truss_, sizeof(int) * m_);
    auto cn_copy = new vector<int>[m_];
    for (int i = 0; i < m_; i++) cn_copy[i] = cn_[i];

    queue<int> q;
    int* cnt = new int[k_max_ + 1];
    for (int t_e = t_ - 1; t_e >= 0; --t_e) {
        for (int i = edges_idx_[t_e]; i < edges_idx_[t_e + 1]; ++i) {
            int u = edges_t_[i].first;
            int v = edges_t_[i].second;
            int eid = nbr_[u][v];

            if (--nbr_cnt_[eid] != 0) continue;

            if (!visit_[eid]) {
                q.push(eid);
                visit_[eid] = true;
            }

            int trussness = truss_[eid];
            for (int w : cn_[eid]) {
                int e1 = nbr_[u][w];
                int e2 = nbr_[v][w];

                // process (u,w)
                cn_[e1].erase(lower_bound(cn_[e1].begin(), cn_[e1].end(), v));  // 在cn(u,w) 删除v
                if (truss_[e1] <= trussness && truss_[e1] <= truss_[e2]) {
                    --td_[e1];
                    if (td_[e1] < truss_[e1] - 2 && !visit_[e1]) {
                        q.push(e1);
                        visit_[e1] = true;
                    }
                }

                // process (v,w)
                cn_[e2].erase(lower_bound(cn_[e2].begin(), cn_[e2].end(), u));  // 在cn(v,w) 删除u
                if (truss_[e2] <= trussness && truss_[e2] <= truss_[e1]) {
                    --td_[e2];
                    if (td_[e2] < truss_[e2] - 2 && !visit_[e2]) {
                        q.push(e2);
                        visit_[e2] = true;
                    }
                }
            }
            cn_[eid].clear();
            // td_[eid] = 0;
        }

        while (!q.empty()) {
            int eid = q.front();
            q.pop();
            visit_[eid] = false;
            int u = edges_[eid].first;
            int v = edges_[eid].second;

            int ot = truss_[eid];
            memset(cnt, 0, sizeof(int) * (ot + 1));

            // 重新计算 truss(u,v) & td_[u][v]
            for (const int w : cn_[eid]) {
                int e1 = nbr_[u][w];
                int e2 = nbr_[v][w];
                ++cnt[min({ot, truss_[e1], truss_[e2]})];
            }
            int td = 0;
            for (int k = ot; k >= 2; --k) {
                td += cnt[k];
                if (td >= k - 2) {
                    truss_[eid] = k;
                    break;
                }
            }
            td_[eid] = td;
            int nt = truss_[eid];

            // 添加受影响的边进入队列
            for (const int w : cn_[eid]) {
                int e1 = nbr_[u][w];
                int e2 = nbr_[v][w];

                // process (u,w)
                if (nt < truss_[e1] && truss_[e1] <= ot && truss_[e1] <= truss_[e2]) {
                    --td_[e1];
                    if (td_[e1] < truss_[e1] - 2 && !visit_[e1]) {
                        q.push(e1);
                        visit_[e1] = true;
                    }
                }

                // process (v,w)
                if (nt < truss_[e2] && truss_[e2] <= ot && truss_[e2] <= truss_[e1]) {
                    --td_[e2];
                    if (td_[e2] < truss_[e2] - 2 && !visit_[e2]) {
                        q.push(e2);
                        visit_[e2] = true;
                    }
                }
            }

            // 更新索引
            for (int k = ot; k > nt; --k) {  // nt_min = 2
                truss_t_[eid][k - 2].emplace_back(make_pair(0, t_e));
            }
        }
    }

    // 恢复数据
    for (int i = 0; i < m_; i++) cn_[i] = cn_copy[i];
    memcpy(nbr_cnt_, nbr_cnt_copy, sizeof(int) * m_);
    memcpy(truss_, truss_copy, sizeof(int) * m_);

    delete[] cnt;
    delete[] cn_copy;
    delete[] nbr_cnt_copy;
    delete[] truss_copy;
}

// 每次删除 k-1 的边，对于 k 边不影响
void Graph::init_tt_cnt(const int& k) {
    memset(tt_cnt_, 0, sizeof(int) * m_);  // 必须置零

    // 删除 truss < k 边，删除这些边不影响 k-truss子图
    for (int eid = 0; eid < m_; ++eid) {
        if (!truss_[eid]) continue;  // 先前 k 迭代已经被处理过了
        if (truss_[eid] < k) {
            int u = edges_[eid].first;
            int v = edges_[eid].second;
            for (const auto& w : cn_[eid]) {
                int e1 = nbr_[u][w];
                int e2 = nbr_[v][w];
                cn_[e1].erase(lower_bound(cn_[e1].begin(), cn_[e1].end(), v));
                cn_[e2].erase(lower_bound(cn_[e2].begin(), cn_[e2].end(), u));
            }
            truss_[eid] = 0;
        }
    }

    // 剩下边 满足 turss >= k
    for (int eid = 0; eid < m_; ++eid) {
        if (truss_[eid] < k) continue;
        int truss_time = truss_t_[eid][k - 2].back().second;  // 即初始化索引的第一个时间窗口
        int u = edges_[eid].first;
        int v = edges_[eid].second;
        for (const int& w : cn_[eid]) {
            int tt1 = truss_t_[nbr_[u][w]][k - 2].back().second;
            int tt2 = truss_t_[nbr_[v][w]][k - 2].back().second;
            if (tt1 <= truss_time && tt2 <= truss_time) ++tt_cnt_[eid];
        }
        if (tt_cnt_[eid] < k - 2) {
            printf("!!!Wrong: u=%d, v=%d, tt_cnt_=%d, k=%d\n", u, v, tt_cnt_[eid], k);
            exit(-1);
        }
    }
}

void Graph::index() {
#ifdef _LINUX_
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif

    if (visit_ == nullptr) visit_ = new bool[m_]{};
    if (truss_t_ == nullptr) truss_t_ = new vector<vector<pair<int, int>>>[m_];

    printf("starting truss decomposition...\n");
    compute_common_neighbors();
    truss_decomposition();

    // 提前申请空间
    for (int i = 0; i < m_; ++i) {
        truss_t_[i].resize(truss_[i] - 1);  // "-2+1" truss >= 2 索引从下标为1的位置开始
    }

    printf("Initialize truss time...\n");
    init_truss_deg();
    init_truss_time();

    if (tt_cnt_ == nullptr) tt_cnt_ = td_;  // 地址复用

    auto nbr_cnt_copy = new int[m_];
    auto truss_copy = new int[m_];
    auto cn_copy = new vector<int>[m_];
    
    queue<int> q;
    for (int k = 3; k <= k_max_; ++k) {
        printf("---- Iteration k = %d/%d ----\n", k, k_max_);

        init_tt_cnt(k);
        // 备份数据
        for (int i = 0; i < m_; i++) cn_copy[i] = cn_[i];
        memcpy(nbr_cnt_copy, nbr_cnt_, sizeof(int) * m_);
        memcpy(truss_copy, truss_, sizeof(int) * m_);

        for (int t_s = 1; t_s < t_; ++t_s) {
            if (t_s % 10000 == 0) printf("t = %d.\n", t_s);
            for (int i = edges_idx_[t_s - 1]; i < edges_idx_[t_s]; ++i) {
                int u = edges_t_[i].first;
                int v = edges_t_[i].second;
                int eid = nbr_[u][v];

                --nbr_cnt_[eid];  // 无论如何 边数减一

                if (truss_[eid] < k || (!truss_t_[eid][k - 2].empty() && truss_t_[eid][k - 2].back().second == t_)) continue;

                int truss_time = truss_t_[eid][k - 2].back().second;
                if (!nbr_cnt_[eid]) {  // 完全删除
                    if (!visit_[eid]) {
                        q.push(eid);
                        visit_[eid] = true;
                    }
                    for (const auto& w : cn_[eid]) {
                        int e1 = nbr_[u][w];
                        int e2 = nbr_[v][w];
                        cn_[e1].erase(lower_bound(cn_[e1].begin(), cn_[e1].end(), v));
                        cn_[e2].erase(lower_bound(cn_[e2].begin(), cn_[e2].end(), u));

                        if (truss_[e1] >= k && truss_[e2] >= k) {
                            int tt1 = truss_t_[e1][k - 2].back().second;
                            int tt2 = truss_t_[e2][k - 2].back().second;
                            if (truss_time <= tt1 && tt2 <= tt1) {
                                --tt_cnt_[e1];
                                if (tt_cnt_[e1] < k - 2 && !visit_[e1]) {
                                    q.push(e1);
                                    visit_[e1] = true;
                                }
                            }
                            if (truss_time <= tt2 && tt1 << tt2) {
                                --tt_cnt_[e2];
                                if (tt_cnt_[e2] < k - 2 && !visit_[e2]) {
                                    q.push(e2);
                                    visit_[e2] = true;
                                }
                            }
                        }
                    }
                    cn_[eid].clear();
                    // truss_[eid] = 0;
                } else {
                    // 部分删除
                    int aaa = nbr_cnt_[eid];
                    int szz = nbr_t_[u][v].size();
                    int new_e_t = nbr_t_[u][v][nbr_t_[u][v].size() - aaa];
                    if (new_e_t >= truss_time) { // ???
                        if (!visit_[eid]) {
                            q.push(eid);
                            visit_[eid] = true;
                        }
                        for (const auto& w : cn_[eid]) {
                            int e1 = nbr_[u][w];
                            int e2 = nbr_[v][w];
                            if (truss_[e1] >= k && truss_[e2] >= k) {
                                int tt1 = truss_t_[e1][k - 2].back().second;
                                int tt2 = truss_t_[e2][k - 2].back().second;

                                if (truss_time <= tt1 && tt2 <= tt1 && new_e_t > tt1) {
                                    --tt_cnt_[e1];
                                    if (tt_cnt_[e1] < k - 2 && !visit_[e1]) {
                                        q.push(e1);
                                        visit_[e1] = true;
                                    }
                                }
                                if (truss_time <= tt2 && tt1 <= tt2 && new_e_t > tt2) {
                                    --tt_cnt_[e2];
                                    if (tt_cnt_[e2] < k - 2 && !visit_[e2]) {
                                        q.push(e2);
                                        visit_[e2] = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            while (!q.empty()) {
                int eid = q.front();
                q.pop();
                visit_[eid] = false;

                int u = edges_[eid].first;
                int v = edges_[eid].second;
                int old_t = truss_t_[eid][k - 2].back().second;

                // 更新 truss time
                vector<int> tt_arr;
                int e_t = t_;
                if (nbr_cnt_[eid]) e_t = nbr_t_[u][v][nbr_t_[u][v].size() - nbr_cnt_[eid]];
                for (const auto& w : cn_[eid]) {
                    int e1 = nbr_[u][w];
                    int e2 = nbr_[v][w];
                    if (truss_[e1] >= k && truss_[e2] >= k) {
                        int tt1 = truss_t_[e1][k - 2].back().second;
                        int tt2 = truss_t_[e2][k - 2].back().second;
                        tt_arr.emplace_back(max({e_t, tt1, tt2}));
                    }
                }
                int new_t = t_;
                if (tt_arr.size() >= k - 2) {
                    nth_element(tt_arr.begin(), tt_arr.begin() + k - 2 - 1, tt_arr.end());
                    new_t = tt_arr[k - 2 - 1];
                } else {
                    truss_[eid] = 0;  // truss < k
                }

                // 更新 tt_cnt_
                tt_cnt_[eid] = 0;
                if (tt_arr.size() >= k - 2) {
                    for (int& t : tt_arr) {
                        if (t <= new_t) ++tt_cnt_[eid];
                    }
                }

                // 更新 index
                if (truss_t_[eid][k - 2].back().first == t_s) {
                    truss_t_[eid][k - 2].back().second = new_t;
                } else if (truss_t_[eid][k - 2].back().second < new_t) {
                    truss_t_[eid][k - 2].emplace_back(make_pair(t_s, new_t));
                }

                // 添加受影响的边入队
                for (const auto& w : cn_[eid]) {
                    int e1 = nbr_[u][w];
                    int e2 = nbr_[v][w];

                    if (truss_[e1] >= k && truss_[e2] >= k) {
                        int tt1 = truss_t_[e1][k - 2].back().second;
                        int tt2 = truss_t_[e2][k - 2].back().second;
                        if (old_t <= tt1 && tt2 <= tt1 && new_t > tt1) {
                            --tt_cnt_[e1];
                            if (tt_cnt_[e1] < k - 2 && !visit_[e1]) {
                                q.push(e1);
                                visit_[e1] = true;
                            }
                        }
                        if (old_t <= tt2 && tt1 <= tt2 && new_t > tt2) {
                            --tt_cnt_[e2];
                            if (tt_cnt_[e2] < k - 2 && !visit_[e2]) {
                                q.push(e2);
                                visit_[e2] = true;
                            }
                        }
                    }
                }
            }
        }

        // 恢复数据
        for (int i = 0; i < m_; i++) cn_[i] = cn_copy[i];
        memcpy(nbr_cnt_, nbr_cnt_copy, sizeof(int) * m_);
        memcpy(truss_, truss_copy, sizeof(int) * m_);
    }
    delete[] cn_copy;
    delete[] nbr_cnt_copy;
    delete[] truss_copy;

#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec) / 1000;
    printf("Running time: %lld s, %lld mins\n", t_msec / 1000, t_msec / 1000 / 60);
    if (log_f_ != nullptr) fprintf(log_f_, "Indexing time: %lld s, %lld mins\n", t_msec / 1000, t_msec / 1000 / 60);

#else
    clock_t end = clock();
    printf("Running time: %.2f s, %.2f min\n", (double)(end - start) / CLOCKS_PER_SEC,
           (double)(end - start) / CLOCKS_PER_SEC / 60);
    if (log_f_ != nullptr)
        fprintf(log_f_, "Running time: %.2f s, %.2f min\n", (double)(end - start) / CLOCKS_PER_SEC,
                (double)(end - start) / CLOCKS_PER_SEC / 60);
#endif
    print_idx_size();
}

void Graph::write_idx_txt(const string& path) {
    ofstream fout(path, ios::out);  // out 文件清空重写
    for (int i = 0; i < m_; ++i) {
        int u = edges_[i].first;
        int v = edges_[i].second;
        fout << "\n --------(" << u << "," << v << ")-------" << endl;
        for (int k = 1; k < truss_t_[i].size(); ++k) {
            fout << "truss = " << k + 2 << endl;
            for (const auto& p : truss_t_[i][k]) {
                fout << "[" << p.first << "," << p.second << "] ";
            }
            fout << endl;
        }
    }
}

void Graph::write_idx(const string& path) {
    auto fp = fopen(path.c_str(), "wb");
    fwrite(&m_, sizeof(unsigned int), 1, fp);  // edge

    for (int i = 0; i < m_; ++i) {
        int sz = truss_t_[i].size();
        fwrite(&sz, sizeof(int), 1, fp);
        for (int k = 1; k < sz; ++k) {  // 从小标为 1 开始写入
            int cnt = truss_t_[i][k].size();
            fwrite(&cnt, sizeof(int), 1, fp);
            fwrite(&truss_t_[i][k][0], sizeof(pair<int, int>) * cnt, 1, fp);
        }
    }

    fclose(fp);
    printf("Write index!\n");
}

void Graph::load_idx(const string& path) {
    auto fp = fopen(path.c_str(), "rb");
    if (fp == NULL) {
        printf("cannot open file!\n");
        exit(1);
    }

    fread(&m_, sizeof(unsigned int), 1, fp);
    k_max_ = 0;  // ??? 修改

    if (truss_t_ == nullptr) truss_t_ = new vector<vector<pair<int, int>>>[m_];

    for (int i = 0; i < m_; i++) {
        int sz;
        fread(&sz, sizeof(int), 1, fp);
        k_max_ = max(k_max_, sz - 1 + 2);
        truss_t_[i].resize(sz);
        for (int k = 1; k < sz; ++k) {
            int cnt;
            fread(&cnt, sizeof(int), 1, fp);
            truss_t_[i][k].resize(cnt);
            fread(&truss_t_[i][k][0], cnt * sizeof(pair<int, int>), 1, fp);
        }
    }

    fclose(fp);
    printf("load index.\n");
}

void Graph::print_graph_size() {
    printf("Graph size: %.2f MB.\n", (float)edges_t_.size() * sizeof(int) * 3 / 1024 / 1024);
    if (log_f_ != nullptr)
        fprintf(log_f_, "Graph size: %.2f MB.\n", (float)edges_t_.size() * sizeof(int) * 3 / 1024 / 1024);
}

void Graph::print_idx_size() {
    idx_size_ = 0;
    idx_size_ += sizeof(int);
    idx_size_ += sizeof(int) * m_;  // m_ 个指针

    double average_t = 0;  // 平均时间窗口数量
    long average_truss = 0;
    int max_t = 0;  // 最多时间窗口数量

    for (int i = 0; i < m_; ++i) {
        for (int k = 1; k < truss_t_[i].size(); ++k) {
            int sz = truss_t_[i][k].size();
            idx_size_ += sz * sizeof(pair<int, int>);
            average_t += sz;
            max_t = max(max_t, sz);
        }
        average_truss += truss_t_[i].size() + 1;  // +2-1
    }

    printf("Index size: %.2f MB.\n", (float)idx_size_ / 1024 / 1024);
    printf("Average truss = %.2f, max k-truss = %d, average T = %.2f, max T = %d.\n", double(average_truss) / m_, k_max_, average_t / m_, max_t);
    if (log_f_ != nullptr) {
        fprintf(log_f_, "Index size: %.2f MB\n", (float)idx_size_ / 1024 / 1024);
        fprintf(log_f_, "Average truss = %.2f, max k-truss = %d, average T = %.2f, max T = %d.\n", double(average_truss) / m_, k_max_, average_t / m_, max_t);
    }
}

void Graph::init_log(const string& log_path) {
    log_f_ = fopen(log_path.c_str(), "a");  // 追加
    fprintf(log_f_, "\n\n=====================================\n");
    time_t now = time(0);
    fprintf(log_f_, "%s\n", ctime(&now));
}

bool Graph::query(int u, int v, int t_s, int t_e, int k) {
    int eid = nbr_[u][v];
    if (truss_t_[eid].size() - 1 < k - 2) return false;  // 下标0位置不存数据
    auto it = upper_bound(truss_t_[eid][k - 2].begin(), truss_t_[eid][k - 2].end(), make_pair(t_s, t_e), cmp1);
    --it;
    return it->second <= t_e;
}

int Graph::query_all(int t_s, int t_e, int k) {
    int res = 0;
    unordered_set<int>* nbr = new unordered_set<int>[n_];
    for (int i = edges_idx_[t_s]; i < edges_idx_[t_e + 1]; ++i) {
        nbr[edges_t_[i].first].insert(edges_t_[i].second);
    }
    // printf("====== query: [%d, %d]\n", t_s, t_e);
    for (int u = 0; u < n_; ++u) {
        for (const auto& v : nbr[u]) {
            if (query(u, v, t_s, t_e, k)) {
                ++res;
                // printf("e = (%d, %d), k = %d\n", u, v, k);
            }
        }
    }
    delete[] nbr;
    return res;
}

int Graph::online_query(int t_s, int t_e, int k, int& snapshot_m) {
    // 排序去重
    vector<pair<int, int>> edges_arr(edges_t_.begin() + edges_idx_[t_s], edges_t_.begin() + edges_idx_[t_e + 1]);
    sort(edges_arr.begin(), edges_arr.end());
    edges_arr.erase(unique(edges_arr.begin(), edges_arr.end()), edges_arr.end());

    snapshot_m = edges_arr.size();

    unordered_map<int, int>* nbr = new unordered_map<int, int>[n_];
    for (int i = 0; i < edges_arr.size(); ++i) {
        nbr[edges_arr[i].first].insert(make_pair(edges_arr[i].second, i));
        nbr[edges_arr[i].second].insert(make_pair(edges_arr[i].first, i));
    }
    int m = edges_arr.size();
    int res = m;

    queue<int> q;

    vector<int>* cn = new vector<int>[m] {};
    vector<bool> del(m, false);
    for (int eid = 0; eid < m; ++eid) {
        int u = edges_arr[eid].first;
        int v = edges_arr[eid].second;
        if (nbr[u].size() > nbr[v].size()) swap(u, v);
        for (const auto& p : nbr[u]) {
            int w = p.first;
            if (nbr[v].find(w) != nbr[v].end()) {
                cn[eid].emplace_back(w);
            }
        }
        sort(cn[eid].begin(), cn[eid].end());
        if (cn[eid].size() < k - 2) {
            q.push(eid);
        }
    }

    while (!q.empty()) {
        int eid = q.front();
        q.pop();

        if (del[eid]) continue;
        del[eid] = true;
        --res;
        int u = edges_arr[eid].first;
        int v = edges_arr[eid].second;
        for (int w : cn[eid]) {
            int e1 = nbr[u][w];
            int e2 = nbr[v][w];

            cn[e1].erase(lower_bound(cn[e1].begin(), cn[e1].end(), v));
            if (cn[e1].size() < k - 2) q.push(e1);

            cn[e2].erase(lower_bound(cn[e2].begin(), cn[e2].end(), u));
            if (cn[e2].size() < k - 2) q.push(e2);
        }
    }

    delete[] nbr;
    delete[] cn;
    return res;
}

bool cmp1(const pair<int, int>& a, const pair<int, int>& b) {
    return a.first < b.first;
}

bool cmp2(const pair<int, int>& a, const pair<int, int>& b) {
    return a.second < b.second;
}
