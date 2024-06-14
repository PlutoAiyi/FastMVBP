#pragma once
#include <math.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "LinearHeap.hpp"
#include "MBitSet.hpp"
#include "RandList.hpp"
struct SubGraph {
    vector<udi> left, right;

    SubGraph() {}

    SubGraph(const vector<udi>& _left, const vector<udi>& _right)
        : left(_left), right(_right) {}

    inline udi size() const { return left.size() + right.size(); }

    inline bool operator<(const SubGraph& oth) const {
        return size() > oth.size();
    }
    bool operator==(const SubGraph& rhs) const {
        return right == rhs.right && left == rhs.left;
    }
    bool contains(SubGraph& rhs) { 
        bool flag = true;
        for (int i = 0; i < rhs.left.size(); i++) {
            if (rhs.left[i] != left[i]) {
                flag = false;
                break;
            }
        }
        if (flag && right.size() >= rhs.right.size()) {
            for (int i = 0; i < rhs.right.size(); i++) {
                if (rhs.right[i] != right[i]) {
                    flag = false;
                    break;
                }
            }
        }
        return flag;
    }
    void sort() noexcept {
        std::sort(left.begin(), left.end());
        std::sort(right.begin(), right.end());
    }

    void print() const noexcept {
        using namespace std;
        cout << "|V|: " << size() <<  " |L|: " << left.size() <<  " |R|: " << right.size()<< endl;
        cout << "L: ";
        for (udi u : left) {
            cout << u << " ";
        }
        cout << endl;
        cout << "R: ";
        for (udi u : right) {
            cout << u << " ";
        }
        cout << endl;
    }
};

class Graph {
   public:
    udi vertex_tot;
    udi edge_tot;
    udi L_size;
    MBitSet in_L;
    std::vector<udi> original_vertex_number;
    std::vector<std::vector<udi>::const_iterator> iter;
    std::vector<udi> edge;
    std::vector<MBitSet> bit_g;

    // 用于遍历某个顶点的所有的边
    struct Vertex {
        udi id;
        const std::vector<std::vector<udi>::const_iterator>& iter;

        std::vector<udi>::const_iterator begin() const noexcept { return iter[id]; }

        std::vector<udi>::const_iterator end() const noexcept { return iter[id + 1]; }

        udi size() const noexcept { return iter[id + 1] - iter[id]; }
    };
    Vertex operator[](const udi& id) const noexcept { return Vertex{id, iter}; }

    Graph() : vertex_tot(0), edge_tot(0), L_size(0){};
    Graph(const Graph&) = delete;
    Graph& operator=(const Graph &) = delete;
    Graph(Graph&&) = default;
    Graph& operator=(Graph&&) = default;


    Graph(std::string dataset) {
        using namespace std;
        ifstream read(dataset);
        if (!read) {
            throw runtime_error("Error: File '" + dataset + "' not found.");
        }
        read >> vertex_tot >> L_size >> edge_tot;
        string temp;
        getline(read, temp);

        in_L = MBitSet(vertex_tot);
        for (udi i = 0; i < L_size; i++) {
            in_L.set(i);
        }
        edge.reserve(edge_tot * 2);
        original_vertex_number.reserve(vertex_tot);
        vector<udi> deg(vertex_tot, 0);
        udi u,v;
        while (!read.eof()) {
            getline(read, temp);
            istringstream line(temp);
            udi u, v;
            line >> u;
            while (line >> v) {
                deg[u]++;
                edge.emplace_back(v);
            }
        }
        iter.reserve(2 * edge_tot);
        udi sum = 0;
        for (udi i = 0; i < vertex_tot; i++) {
            iter.emplace_back(edge.begin() + sum);
            sum += deg[i];
        }
        iter.emplace_back(edge.end());
        for (udi u = 0; u < vertex_tot; u++) {
            original_vertex_number.emplace_back(u);
        }
    }

    void print() const noexcept {
        using namespace std;
        cout << "|V|:" << vertex_tot << endl;
        cout << "|L|:" << L_size << endl;
        cout << "|E|:" << edge_tot << endl;
        for (udi i = 0; i < edge.size(); i++) {
            cout << i << ' ';
            for (int v : (*this)[i]) {
                cout << v << ' ';
            }
            cout << endl;
        }
    }

    inline bool is_neighbor(udi u, udi v) const noexcept { return bit_g[u][v]; }

    inline bool is_L(udi u, bool exchanged) const noexcept {
        return in_L[u] ^ exchanged;
    }

    bool check(SubGraph sg, udi k) const noexcept {
        for (udi v : sg.left) {
            udi cnt = 0;
            for (udi u : sg.right) {
                if (!is_neighbor(v, u)) {
                    cnt++;
                }
            }
            if (cnt > k) {
                return false;
            }
        }
        for (udi v : sg.right) {
            udi cnt = 0;
            for (udi u : sg.left) {
                if (!is_neighbor(v, u)) {
                    cnt++;
                }
            }
            if (cnt > k) {
                return false;
            }
        }

        return true;
    }

    void reduction(udi lb_L, udi lb_R) {}
    Graph _2_hop_degeneracy_ordering_and_reduction(udi lb_L, udi lb_R, udi k,
                                                   RandList& delete_v,
                                                   bool exchange) const {
        std::vector<int> deg(vertex_tot);
        std::vector<udi> visit(vertex_tot);
        for (udi i = 0; i < vertex_tot; i++) { 
            deg[i] = (*this)[i].size();
            visit[i] = 1;
        }
        bool flag = true;
        while (flag) {
            flag = false;
            for (udi i = 0; i < vertex_tot; i++) {
                if (visit[i] == 0) {
                    continue;
                } 
                if ((in_L[i] && deg[i] + k >= lb_R) ||
                    (!in_L[i] && deg[i] + k >= lb_L)) {
                    continue;
                }
                flag = true;
                visit[i] = 0;
                // is_del.set(i);
                for (udi v : (*this)[i]) {
                    deg[v]--;
                }
            }
        }
        LinearHeap heap(vertex_tot, vertex_tot);
        std::vector<udi> temp(vertex_tot, vertex_tot);
        std::vector<udi> r_temp(vertex_tot, vertex_tot);
        // std::fill(nd.begin(), nd.end(), 0);
        Graph tg;
        udi new_num=0;
            for(udi i = 0; i < vertex_tot; i++){
                if(visit[i]){
                    temp[i]=new_num;
                    r_temp[new_num++]=i;
                }
            }
            tg.vertex_tot=new_num;
        tg.bit_g.resize(tg.vertex_tot, MBitSet(tg.vertex_tot));
        udi tot_L=0;
        udi tot_R=0;
        long long edges=0;
                std::vector<udi> nd(new_num,0);
        std:: vector<int> mapping_L(new_num,0);
        std:: vector<int> mapping_R(new_num,0);
        for (udi i = 0; i < new_num; i++) {
            udi u = r_temp[i];
            if(in_L[u]) {mapping_L[i]=tot_L++;}
            else mapping_R[i]=tot_R++;
            for (udi v : (*this)[u]) {
                if (visit[v]) {
                    tg.edge.emplace_back(temp[v]);
                    tg.bit_g[i].set(temp[v]);
                    nd[i]++;
                    edges++;
                }
            }
        }
        udi sum_t=0;
        for (udi i = 0; i < new_num; i++) {
            // udi u = r_temp[i];
            tg.iter.emplace_back(tg.edge.begin() + sum_t);
            sum_t += nd[i];
        }
        tg.iter.emplace_back(tg.edge.end());
        MBitSet is_cal(new_num);
        std::vector<int> DG(new_num, 0);
        // std::vector<int> DG_o2(new_num, 0);
        bool exact=true;
        if(new_num !=0)
        // std::cout<<2*edges/new_num<<std::endl;
        if(new_num !=0 && (2*edges/new_num) >100 && new_num>10000) exact=false;
        if(new_num !=0 && (2*edges/new_num) >20 && new_num>30000) exact=false;
        if(new_num !=0 && (2*edges/new_num) >10 && new_num>100000) exact=false;
        std::vector<std::vector<int>> con_L;
        std::vector<std::vector<int>> con_R;
        exact=true;
        if(exact){
                con_L.resize(tot_L, std::vector<int>(tot_L, 0));
                con_R.resize(tot_R, std::vector<int>(tot_R, 0));
        }
        if(exact){
            for(int i=0;i<new_num;i++){
                for (auto it1 = tg[i].begin(); it1 != tg[i].end(); ++it1) {
                    for (auto it2 = std::next(it1); it2 != tg[i].end(); ++it2) {
                        udi u = *it1;
                        udi v = *it2;
                        if(in_L[r_temp[u]]){
                            con_L[mapping_L[u]][mapping_L[v]]++; 
                            con_L[mapping_L[v]][mapping_L[u]]++;
                        }
                        else{
                            con_R[mapping_R[v]][mapping_R[u]]++;
                            con_R[mapping_R[u]][mapping_R[v]]++;
                        }
                    }
                }
            }
        }
        if(exact){
            for(int i=0;i<new_num;i++){
                int threshold=0;
                if(in_L[r_temp[i]]){
                    threshold=lb_R-2*k;
                }else{
                    threshold=lb_L-2*k;
                }
                is_cal.reset_all();
                is_cal.set(i);
                    // int j = 0;
                    DG[i]+=nd[i];
                    for (udi u : tg[i]) {
                            for (udi v : tg[u]) {
                                if ( in_L[r_temp[i]]&& !is_cal[v] && con_L[mapping_L[i]][mapping_L[v]]>= threshold) {
                                    DG[i]++;
                                    is_cal.set(v);
                                }else if(!in_L[r_temp[i]]&& !is_cal[v] && con_R[mapping_R[v]][mapping_R[i]]>= threshold){
                                    DG[i]++;
                                    is_cal.set(v);
                                }
                            }
                    }
                    heap.add(DG[i], i);
            }
        }else{
            for(int i=0;i<new_num;i++){
                    heap.add(nd[i], i);
            }
        }
        Graph ng;
        std::vector<udi> to(vertex_tot,  vertex_tot);
        int maxDG2L = 0;
        int maxDG2R = 0;
        int max_L = 0;
        MBitSet is_vi(new_num);
        MBitSet is_w(new_num);
        int num_L=0;
        int num_R=0;
        if(exact){
            while (!heap.empty()) {
                udi u = heap.pop_min();
                if(in_L[r_temp[u]]) num_L++;
                if(!in_L[r_temp[u]]) num_R++;
                if(tot_L+1-num_L<lb_L || tot_R+1-num_R<lb_R) {
                    u=r_temp[u];
                    to[u] = ng.vertex_tot++;
                    ng.original_vertex_number.emplace_back(original_vertex_number[u]);
                    continue;
                }
                if(!DG[u]) {
                u=r_temp[u];
                to[u] = ng.vertex_tot++;
                ng.original_vertex_number.emplace_back(original_vertex_number[u]);
                continue;
                }
                //
                if (in_L[r_temp[u]] && DG[u] > maxDG2L) {
                    maxDG2L = DG[u];
                }
                if (!in_L[r_temp[u]] && DG[u] > maxDG2R) maxDG2R = DG[u];
                is_vi.set(u);
                is_cal.reset_all();
                is_cal.set(u);
                is_w.reset_all();
                    for (udi w: tg[u]) {
                        if (!is_vi[w]) {
                            is_w.set(w);
                                heap.del(DG[w], w);
                                DG[w]--;
                                for(udi v:tg[u]){
                                    if(!is_w[v] && !is_vi[v]){
                                        int threshold=1;
                                        if(in_L[r_temp[w]]){
                                            threshold=lb_R-2*k;
                                        }else{
                                            threshold=lb_L-2*k;
                                        }
                                        if(in_L[r_temp[w]] && con_L[mapping_L[w]][mapping_L[v]]--==threshold){
                                            DG[w]--;
                                            heap.del(DG[v], v);
                                            DG[v]--;
                                            heap.add(DG[v], v);
                                        }else if(!in_L[r_temp[w]] && con_R[mapping_R[v]][mapping_R[w]]--==threshold){
                                            DG[w]--;
                                            heap.del(DG[v], v);
                                            DG[v]--;
                                            heap.add(DG[v], v);
                                        }
                                        if(in_L[r_temp[w]]){
                                            con_L[mapping_L[v]][mapping_L[w]]--;
                                        }else{
                                            con_R[mapping_R[w]][mapping_R[v]]--;
                                        }
                                    }
                                }
                                heap.add(DG[w], w);
                            int threshold=0;
                            if(in_L[r_temp[u]]){
                                threshold=lb_R-2*k;
                            }else{
                                threshold=lb_L-2*k;
                            }
                            for (udi v : tg[w]) {
                                if (in_L[r_temp[u]] && !is_vi[v] && con_L[mapping_L[u]][mapping_L[v]]>=threshold) {
                                    if (!is_cal[v]) {
                                        is_cal.set(v);
                                        heap.del(DG[v], v); 
                                        DG[v]--;
                                        heap.add(DG[v], v);
                                    }
                                }else if(!in_L[r_temp[u]] && !is_vi[v] && con_R[mapping_R[v]][mapping_R[u]]>=threshold){
                                    if (!is_cal[v]) {
                                        is_cal.set(v);
                                        heap.del(DG[v], v); 
                                        DG[v]--;
                                        heap.add(DG[v], v);
                                    }
                                }
                            }
                        }
                    }
                u=r_temp[u];
                to[u] = ng.vertex_tot++;
                ng.original_vertex_number.emplace_back(original_vertex_number[u]);
            }
        }else{
            while (!heap.empty()) {
                udi u = heap.pop_min();
                if(in_L[r_temp[u]]) num_L++;
                if(!in_L[r_temp[u]]) num_R++;
                if(tot_L+1-num_L<lb_L || tot_R+1-num_R<lb_R) {
                    u=r_temp[u];
                    to[u] = ng.vertex_tot++;
                    ng.original_vertex_number.emplace_back(original_vertex_number[u]);
                    continue;
                }
                if(!nd[u]) {
                u=r_temp[u];
                to[u] = ng.vertex_tot++;
                ng.original_vertex_number.emplace_back(original_vertex_number[u]);
                continue;
                }
                if (in_L[r_temp[u]] && nd[u] > maxDG2L) {
                    maxDG2L = nd[u];
                }
                if (!in_L[r_temp[u]] && nd[u] > maxDG2R) maxDG2R = nd[u];
                is_vi.set(u);
                    for (udi w: tg[u]) {
                        if (!is_vi[w]) {
                                heap.del(nd[w], w);
                                nd[w]--;
                                heap.add(nd[w], w);
                        }
                    }
                u=r_temp[u];
                to[u] = ng.vertex_tot++;
                ng.original_vertex_number.emplace_back(original_vertex_number[u]);
            }
        }
        std::vector<udi> new_deg(ng.vertex_tot);
        std::vector<udi> high_deg(ng.vertex_tot);
        std::fill(new_deg.begin(), new_deg.end(), 0);
        ng.in_L = MBitSet(ng.vertex_tot);
        ng.bit_g.resize(ng.vertex_tot, MBitSet(ng.vertex_tot));
        udi maxdegree = 0;
        udi mindegree = ng.vertex_tot;
        udi L = 0;
        for (udi i = 0; i < ng.vertex_tot; i++) {
            udi u = ng.original_vertex_number[i];
            if (in_L[u]) {
                ng.in_L.set(i);
                L++;
            }
            for (udi v : (*this)[u]) {
                if (to[v] < vertex_tot) {
                    ng.edge.emplace_back(to[v]);
                    ng.bit_g[i].set(to[v]);
                    new_deg[i]++;
                    if (to[v] > i) {
                        high_deg[i]++;
                    }
                }
            }
            if (new_deg[i] > maxdegree) maxdegree = new_deg[i];
            if (new_deg[i] < mindegree) mindegree = new_deg[i];
        }
        if (ng.vertex_tot > 2 * L) {
            exchange = true;
        }
        ng.edge_tot = ng.edge.size();
        udi sum = 0;
        for (udi i = 0; i < ng.vertex_tot; i++) {
            ng.iter.emplace_back(ng.edge.begin() + sum);
            sum += new_deg[i];
        }
        ng.iter.emplace_back(ng.edge.end());
        delete_v.reinit(ng.vertex_tot);
        bool flag_b = true;
        for (udi i = 0; i < ng.vertex_tot; i++) {
            if ((ng.in_L[i] && high_deg[i] + k < lb_R) ||
                (!ng.in_L[i] && high_deg[i] + k < lb_L)) {
                delete_v.add(i);
            }
        }
        return ng;
    }
    inline void get_original_vertex(SubGraph& sg) const {
        for (udi& u : sg.left) {
            u = original_vertex_number[u];
        }
        for (udi& u : sg.right) {
            u = original_vertex_number[u];
        }
    }
};