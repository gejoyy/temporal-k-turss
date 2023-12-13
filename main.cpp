#include "graph.hpp"

#include <random>

using namespace std;

int main(int argc, char* argv[]) {
    // 1:-idx  2:grapg_path  3:index_path
    if (strcmp(argv[1], "-idx") == 0) {
        string graph_path(argv[2]);
        string idx_path(argv[3]);

        auto* g = new Graph();
        g->init_log("./log.txt");
        g->load(graph_path);

        g->index();
        g->write_idx(idx_path);
        // g->write_idx_txt(idx_path);
        delete g;
    }

    // 1:-idx-bl  2:grapg_path  3:index_path
    if (strcmp(argv[1], "-idx-bl") == 0) {
        string graph_path(argv[2]);
        string idx_path(argv[3]);

        auto* g = new Graph();
        g->init_log("./bl-log.txt");
        g->load(graph_path);

        g->index_baseline();
        g->write_idx(idx_path);
        // g->write_idx_txt(idx_path);
        delete g;
    }

    // -q   ./dataset/  ./index/  ./ktruss-log.txt  40  4
    if (strcmp(argv[1], "-q") == 0) {
        string graph_path(argv[2]);           // 图
        string idx_path(argv[3]);             // 索引位置
        string log_path("./ktruss-log.txt");  // 日志位置
        int t_range = atoi(argv[4]);          // 时间范围 百分比
        int cnt = atoi(argv[5]);              // 查几次

        auto* g = new Graph();

        g->load(graph_path);
        g->load_idx(idx_path);
        printf("max k-truss = %d\n", g->k_max_);

        int t_gap = g->t_ * t_range / 100;
        printf("Prepare %d queries\n", cnt);

        vector<pair<int, int>> queries;
        int t_s_max = g->t_ - t_gap - 1;

        random_device rd;
        uniform_int_distribution<int> dist(0, t_s_max - 1);

        for (int i = 0; i < cnt; ++i) {
            int t_s = 0;
            if (t_s_max != 0) t_s = dist(rd);  // 随机
            int t_e = t_s + t_gap;
            queries.emplace_back(make_pair(t_s, t_e));
        }

        FILE* log_f = fopen(log_path.c_str(), "a");
        fprintf(log_f, "\n\n======================= Query Processing =======================\n");
        time_t now = time(0);  // 把 now 转换为字符串形式
        fprintf(log_f, "%s\n", ctime(&now));
        fprintf(log_f, "Index path:%s\n", idx_path.c_str());
        fprintf(log_f, "t_range: %d%%, t_max = %d, k_max = %d,\n", t_range, g->t_, g->k_max_);

        for (int k_range = 2; k_range < 10; k_range += 2) {
            int k = max(g->k_max_ * k_range / 10, 3);
            fprintf(log_f, "\n---------- k_range = %d0%%,k = %d ----------\n", k_range, k);
            printf("\n-------- k_range = %d0%%,k = %d --------\n", k_range, k);
            vector<int> result_record;

            int empty_q = 0;
            int result_size = 0;

#ifdef _LINUX_
            struct timeval t_start, t_end;
            gettimeofday(&t_start, NULL);
#endif

            for (auto& ti : queries) {
                int res = g->query_all(ti.first, ti.second, k);
                result_size += res;
                if (res <= 0) ++empty_q;
                // printf("ts=%d, te=%d, k=%d, res=%d\n", ti.first, ti.second, k, res);
                result_record.emplace_back(res);
            }
#ifdef _LINUX_
            gettimeofday(&t_end, NULL);
            long long t_msec = (t_end.tv_sec - t_start.tv_sec) * 1000 * 1000 + (t_end.tv_usec - t_start.tv_usec);

            float average_result_size = float(result_size) / cnt;
            printf("Average query time: %lld *e-6 s(us)\n", t_msec / cnt);
            printf("Empty result queries: %d, average result size: %.2f\n", empty_q, average_result_size);
            fprintf(log_f, "Average query time: %lld *e-6 s(us)\n", t_msec / cnt);
            fprintf(log_f, "Empty result queries: %d, average result size: %.2f\n", empty_q, average_result_size);
#endif

            empty_q = 0;
            result_size = 0;

            unsigned long long average_snapshot_m = 0;

#ifdef _LINUX_
            gettimeofday(&t_start, NULL);
#endif

            for (auto& ti : queries) {
                int snapshot_edges;
                int res = g->online_query(ti.first, ti.second, k, snapshot_edges);
                average_snapshot_m += snapshot_edges;
                result_size += res;
                if (res <= 0) ++empty_q;
                //    printf("online query ts=%d, te=%d, k=%d, res=%d\n", ti.first, ti.second, k, res);
                result_record.emplace_back(res);
            }

#ifdef _LINUX_
            gettimeofday(&t_end, NULL);
            t_msec = (t_end.tv_sec - t_start.tv_sec) * 1000 * 1000 + (t_end.tv_usec - t_start.tv_usec);

            average_result_size = float(result_size) / cnt;
            printf("\nAverage online query time: %lld *e-6 s(us)\n", t_msec / cnt);
            printf("Empty result queries: %d, average result size: %.2f\n", empty_q, average_result_size);
            fprintf(log_f, "\nAverage online query time: %lld *e-6 s(us)\n", t_msec / cnt);
            fprintf(log_f, "Empty result queries: %d, average result size: %.2f\n", empty_q, average_result_size);

            average_snapshot_m /= cnt;
            printf("\nAverage snapshot edges number: %llu\n", average_snapshot_m);
            fprintf(log_f, "\nAverage snapshot edges number: %llu\n", average_snapshot_m);
#endif

            for (int i = 0; i < queries.size(); ++i) {
                cout << "query_all: " << result_record[i] << "     online_query: " << result_record[i + queries.size()] << endl;
            }
        }

        delete g;
        fclose(log_f);
    }

    // -txt  ./dataset/  ./index/.bin  ./index/.txt
    if (strcmp(argv[1], "-txt") == 0) {
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string write_path(argv[4]);
        auto* g = new Graph();

        g->load(graph_path);
        g->load_idx(idx_path);
        g->write_idx_txt(write_path);
        delete g;
    }
    return 0;
}