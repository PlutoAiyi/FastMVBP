#pragma once
#include <cassert>
#include <complex>
#include <queue>
#include <set>
#include <vector>
#include "Graph.hpp"
#include "MBitSet.hpp"
#include "RandList.hpp"
namespace Framework {
using namespace std;
const udi MAX=4294967295;
udi k, top_k;
udi num_e;
priority_queue<SubGraph> pq_k;
thread_local udi max_lvl=0;
udi global_ubL;
udi global_ubR;
udi current_pivot;
void record(SubGraph sg,const Graph &g,bool exchanged) {
        num_e++;
        if (pq_k.size() < top_k) {
            pq_k.push(sg);
        } else if (pq_k.top().size() < sg.size()) {
            pq_k.pop();
            pq_k.push(sg);
        }
}
udi lb_global() {
    if (pq_k.size() < top_k) {
        return 0;
    }
    return pq_k.top().size();
}
struct BiKplex {
    bool exchanged;
    udi lb_L, lb_R, ub_L, ub_R;
    RandList P_L, P_R, C_L, C_R, X_L, X_R;
    vector<udi> deg_in_P, deg_in_PC;
    bool flag_r,flag_l,flag_r2,flag_l2,is_line;
    vector<udi> reduction_removeL,reduction_removeR;
    udi olvl_r,olvl_l;
    BiKplex(const Graph& g, udi lb_L, udi lb_R, udi ub_L, udi ub_R)
        : exchanged(false),
          lb_L(lb_L),
          lb_R(lb_R),
          ub_L(ub_L),
          ub_R(ub_R),
          P_L(g.vertex_tot),
          P_R(g.vertex_tot),
          C_L(g.vertex_tot),
          C_R(g.vertex_tot),
          X_L(g.vertex_tot),
          X_R(g.vertex_tot),
        deg_in_P(g.vertex_tot),
        deg_in_PC(g.vertex_tot) {
        flag_r=false;
        flag_l=false;
        flag_r2=false;
        flag_l2=false;
        olvl_r=MAX;
        olvl_l=MAX;
        reduction_removeL.reserve(g.vertex_tot);
        reduction_removeR.reserve(g.vertex_tot);
    };
    BiKplex(){

    }
    void reInit(const BiKplex& other){
            lb_L = other.lb_L;
            lb_R = other.lb_R;
            ub_L = other.ub_L;
            ub_R = other.ub_R;
            exchanged=other.exchanged;
            flag_r=false;
            flag_l=false;
            flag_r2=false;
            flag_l2=false;
            olvl_r=MAX;
            olvl_l=MAX;
            P_L=other.P_L;
            P_R=other.P_R;
            C_L=other.C_L;
            C_R=other.C_R;
            X_L=other.X_L;
            X_R=other.X_R;
            deg_in_P=other.deg_in_P;
            deg_in_PC=other.deg_in_PC;
            reduction_removeL=other.reduction_removeL;
            reduction_removeR=other.reduction_removeR;
          }
    void reInit(const BiKplex& other,vector<udi>& nei3,bool flag){
            lb_L = other.lb_L;
            lb_R = other.lb_R;
            ub_L = other.ub_L;
            ub_R = other.ub_R;
            exchanged=other.exchanged;
            flag_r=false;
            flag_l=false;
            flag_r2=false;
            flag_l2=false;
            olvl_r=MAX;
            olvl_l=MAX;
            deg_in_P=other.deg_in_P;
            deg_in_PC=other.deg_in_PC;
            P_L=other.P_L;
            P_R=other.P_R;
            C_L=other.C_L;
            C_R=other.C_R;
            X_L=other.X_L;
            X_R=other.X_R;
            reduction_removeL=other.reduction_removeL;
            reduction_removeR=other.reduction_removeR;
          }
    void init(const BiKplex& other){
            lb_L = other.lb_L;
            lb_R = other.lb_R;
            ub_L = other.ub_L;
            ub_R = other.ub_R;
            exchanged=other.exchanged;
            flag_r=false;
            flag_l=false;
            flag_r2=false;
            flag_l2=false;
            olvl_r=MAX;
            olvl_l=MAX;
            P_L=other.P_L;
            P_R=other.P_R;
            C_L=other.C_L;
            C_R=other.C_R;
            X_L=other.X_L;
            X_R=other.X_R;
            deg_in_P=other.deg_in_P;
            deg_in_PC=other.deg_in_PC;
            reduction_removeL=other.reduction_removeL;
            reduction_removeR=other.reduction_removeR;
          }
    BiKplex& operator=(const BiKplex& other) {
        lb_L = other.lb_L;
        lb_R = other.lb_R;
        ub_L = other.ub_L;
        ub_R = other.ub_R;
        exchanged=other.exchanged;
        deg_in_P=other.deg_in_P;
        deg_in_PC=other.deg_in_PC;
        P_L=other.P_L;
        P_R=other.P_R;
        C_L=other.C_L;
        C_R=other.C_R;
        X_L=other.X_L;
        X_R=other.X_R;
        olvl_l=other.olvl_l;
        olvl_r=other.olvl_r;
        reduction_removeL=other.reduction_removeL;
        reduction_removeR=other.reduction_removeR;
        flag_r=other.flag_r;
        flag_l=other.flag_l;
        flag_r2=other.flag_r2;
        flag_l2=other.flag_l2;
        return *this;
    };
    void swap(RandList& lhs, RandList& rhs) noexcept {
        lhs.swap(rhs);
    }
    inline void exchange() {
        exchanged = !exchanged;
        std::swap(lb_L, lb_R);
        std::swap(ub_L, ub_R);
        swap(P_L, P_R);
        swap(C_L, C_R);
        swap(X_L, X_R);
        std::swap(reduction_removeL,reduction_removeR);
        std::swap(olvl_l,olvl_r);
        std::swap(flag_l2,flag_r2);
        std::swap(flag_l,flag_r);
        // swap(unlim_L, unlim_R);
    }
    inline void exchange(vector<udi>& remove_CL,vector<udi>& remove_CR,vector<udi>& remove_XL,vector<udi>& remove_XR) {
        exchanged = !exchanged;
        std::swap(lb_L, lb_R);
        std::swap(ub_L, ub_R);
        swap(P_L, P_R);
        swap(C_L, C_R);
        swap(X_L, X_R);
        // swap(unlim_L, unlim_R);
        std::swap(remove_CL, remove_CR);
        std::swap(remove_XL, remove_XR);
        std::swap(reduction_removeL,reduction_removeR);
        std::swap(olvl_l,olvl_r);
        std::swap(flag_l2,flag_r2);
        std::swap(flag_l,flag_r);
    }
    SubGraph get_subgraph() const {
        if (!exchanged) {
            return SubGraph(P_L.vlist, P_R.vlist);
        } else {
            return SubGraph(P_R.vlist, P_L.vlist);
        }
    }
    SubGraph get_subgraph_left() const {
        vector<udi> l = P_L.vlist;
        vector<udi> r = P_R.vlist;
        l.insert(l.end(), C_L.vlist.begin(), C_L.vlist.end());
        if (!exchanged) {
            return SubGraph(l, r);
        } else {
            return SubGraph(r, l);
        }
    }
    SubGraph get_subgraph_right() const {
        vector<udi> l = P_L.vlist;
        vector<udi> r = P_R.vlist;
        r.insert(r.end(), C_R.vlist.begin(), C_R.vlist.end());
        if (!exchanged) {
            return SubGraph(l, r);
        } else {
            return SubGraph(r, l);
        }
    }
    SubGraph get_subgraph_all() const {
        vector<udi> l = P_L.vlist, r = P_R.vlist;
        l.insert(l.end(), C_L.vlist.begin(), C_L.vlist.end());
        r.insert(r.end(), C_R.vlist.begin(), C_R.vlist.end());
        if (!exchanged) {
            return SubGraph(l, r);
        } else {
            return SubGraph(r, l);
        }
    }
    inline void P_L_add(const Graph& g, udi u) noexcept {
        P_L.add(u);
        // P_bit.set(u);
        // PC_bit.set(u);
        for (udi v : g[u]) {
            deg_in_P[v]++;
            deg_in_PC[v]++;
        }
        if (deg_in_P[u] + k == P_R.vnum) {
            auto pred = [&](udi v) { return !g.is_neighbor(u, v); };
            C_R_remove_if(g, pred);
            X_R_remove_if(g, pred);
        } else {
            // unlim_L.add(u);
        }
        for (int i = P_R.vnum - 1; i >= 0; i--) {
            udi v = P_R.vlist[i];
            if (deg_in_P[v] + k == P_L.vnum) {
                auto pred = [&](udi u) { return !g.is_neighbor(v, u); };
                C_L_remove_if(g, pred);
                X_L_remove_if(g, pred);
                // if(unlim_R.contains(v))
                // unlim_R.remove(v);
            }
        }
    }
    inline void P_L_add(const Graph& g, udi u,vector<udi> &remove_CR,vector<udi> &remove_CL,vector<udi> &remove_XR,vector<udi> &remove_XL)noexcept {
        P_L.add(u);

        for (int v : g[u]) {
            deg_in_P[v]++;
            deg_in_PC[v]++;
        }
        if (deg_in_P[u] + k == P_R.vnum) {
            auto pred = [&](int v) { return !g.is_neighbor(u, v); };
            C_R_remove_if_re(g, pred,remove_CR);
            X_R_remove_if_re(g, pred,remove_XR);
        }else{
            // unlim_L.add(u);
        }
        for (int i = P_R.vnum - 1; i >= 0; i--) {
            udi v = P_R.vlist[i];
            if (deg_in_P[v] + k == P_L.vnum) {
                auto pred = [&](udi u) { return !g.is_neighbor(v, u); };
                C_L_remove_if_re(g, pred,remove_CL);
                X_L_remove_if_re(g, pred,remove_XL);
                // if(unlim_R.contains(v))
                // unlim_R.remove(v);
            }
        }
    }

    inline void P_L_remove(const Graph& g, int u) noexcept{
        P_L.remove(u);
        for (int v : g[u]) {
            deg_in_P[v]--;
            deg_in_PC[v]--;
        }
    }
    inline void P_R_add(const Graph& g, udi v,vector<udi> &remove_CR,vector<udi> &remove_CL,vector<udi> &remove_XR,vector<udi> &remove_XL)noexcept {
        P_R.add(v);
        for (int u : g[v]) {
            deg_in_P[u]++;
            deg_in_PC[u]++;
        }
        if (deg_in_P[v] + k == P_L.vnum) {
            auto pred = [&](int u) { return !g.is_neighbor(v, u); };
                C_L_remove_if_re(g, pred,remove_CL);
                X_L_remove_if_re(g, pred,remove_XL);
        } else {
            // unlim_R.add(v);
        }
        for (int i = P_L.vnum - 1; i >= 0; i--) {
            int u = P_L.vlist[i];
            if(g.is_neighbor(v,u)){
                continue;
            }
            if (deg_in_P[u] + k == P_R.vnum) {
                auto pred = [&](int v) { return !g.is_neighbor(u, v); };
                C_R_remove_if_re(g, pred,remove_CR);
                X_R_remove_if_re(g, pred,remove_XR);
                // if(unlim_L.contains(u))
                // unlim_L.remove(u);
            }
        }
    }
inline void CL2PL(const Graph& g, udi v,vector<udi> &remove_CR,vector<udi> &remove_CL,vector<udi> &remove_XR,vector<udi> &remove_XL) noexcept{
        C_L.remove(v);
        P_L.add(v);
        for (int u : g[v]) {
            deg_in_P[u]++;
        }
        if (deg_in_P[v] + k == P_R.vnum) {
            auto pred = [&](udi u) { return !g.is_neighbor(v, u); };
                C_R_remove_if_re(g, pred,remove_CR);
                X_R_remove_if_re(g, pred,remove_XR);
        }
        for (int i = P_R.vnum - 1; i >= 0; i--) {
            udi u = P_R.vlist[i];
            if(g.is_neighbor(v,u)){
                continue;
            }
            if (deg_in_P[u] + k == P_L.vnum) {
                auto pred = [&](udi v) { return !g.is_neighbor(u, v); };
                C_L_remove_if_re(g, pred,remove_CL);
                X_L_remove_if_re(g, pred,remove_XL);
            }
        }
    }
inline void CL2PL(const Graph& g, udi v) noexcept{
        C_L.remove(v);
        P_L.add(v);
        for (int u : g[v]) {
            deg_in_P[u]++;
        }
        if (deg_in_P[v] + k == P_R.vnum) {
            auto pred = [&](udi u) { return !g.is_neighbor(v, u); };
                C_R_remove_if(g, pred);
                X_R_remove_if(g, pred);
        }
        for (int i = P_R.vnum - 1; i >= 0; i--) {
            udi u = P_R.vlist[i];
            if(g.is_neighbor(v,u)){
                continue;
            }
            if (deg_in_P[u] + k == P_L.vnum) {
                auto pred = [&](udi v) { return !g.is_neighbor(u, v); };
                C_L_remove_if(g, pred);
                X_L_remove_if(g, pred);
            }
        }
    }
inline void PL2CL(const Graph& g, udi v)noexcept {
        P_L.remove(v);
        C_L.add(v);
        for (int u : g[v]) {
            deg_in_P[u]--;
        }
    }
inline void PR2CR(const Graph& g, udi v) noexcept{
        P_R.remove(v);
        C_R.add(v);
        for (int u : g[v]) {
            deg_in_P[u]--;
        }
    }

inline void CR2PR(const Graph& g, udi v,vector<udi> &remove_CR,vector<udi> &remove_CL,vector<udi> &remove_XR,vector<udi> &remove_XL) noexcept{
        C_R.remove(v);
        P_R.add(v);
        for (int u : g[v]) {
            deg_in_P[u]++;
        }
        if (deg_in_P[v] + k == P_L.vnum) {
            auto pred = [&](int u) { return !g.is_neighbor(v, u); };
                C_L_remove_if_re(g, pred,remove_CL);
                X_L_remove_if_re(g, pred,remove_XL);
        }
        for (int i = P_L.vnum - 1; i >= 0; i--) {
            udi u = P_L.vlist[i];
            if(g.is_neighbor(v,u)){
                continue;
            }
            if (deg_in_P[u] + k == P_R.vnum) {
                auto pred = [&](udi v) { return !g.is_neighbor(u, v); };
                C_R_remove_if_re(g, pred,remove_CR);
                X_R_remove_if_re(g, pred,remove_XR);
            }
        }
    }
    inline void P_L_remove_if(const Graph& g,
                              const auto& pred) noexcept{
        for (int i = P_L.vnum - 1; i >= 0; i--) {
            if (pred(P_L.vlist[i])) {
                P_L_remove(g, P_L.vlist[i]);
            }
        }
    }

    inline void P_R_add(const Graph& g, udi v) noexcept{
        P_R.add(v);
        // P_bit.set(v);
        // PC_bit.set(v);
        for (udi u : g[v]) {
            deg_in_P[u]++;
            deg_in_PC[u]++;
        }
        if (deg_in_P[v] + k == P_L.vnum) {
            auto pred = [&](udi u) { return !g.is_neighbor(v, u); };
            C_L_remove_if(g, pred);
            X_L_remove_if(g, pred);
        }
        for (int i = P_L.vnum - 1; i >= 0; i--) {
            udi u = P_L.vlist[i];
            if(g.is_neighbor(v,u)){
                continue;
            }
            if (deg_in_P[u] + k == P_R.vnum) {
                auto pred = [&](udi v) { return !g.is_neighbor(u, v); };
                C_R_remove_if(g, pred);
                X_R_remove_if(g, pred);
                // if(unlim_L.contains(u))
                // unlim_L.remove(u);
            }
        }
    }

    inline bool C_L_add_and_check(const Graph& g)noexcept {
        udi all = P_L.vnum + C_L.vnum;
        for (int i = C_L.vnum - 1; i >= 0; i--) {
            udi u = C_L.vlist[i];
            C_L_remove(g, u);
            P_L_add(g, u);
            if (P_L.vnum + C_L.vnum < all) {
                return false;
            }
        }
        return true;
    }
    inline bool C_R_add_and_check(const Graph& g) noexcept{
        udi all = P_R.vnum + C_R.vnum;
        for (int i = C_R.vnum - 1; i >= 0; i--) {
            udi u = C_R.vlist[i];
            C_R_remove(g, u);
            P_R_add(g, u);
            if (P_R.vnum + C_R.vnum < all) {
                return false;
            }
        }
        return true;
    }
    inline void P_R_remove(const Graph& g, udi u)noexcept {
        P_R.remove(u);
        // P_bit.reset(u);
        // PC_bit.reset(u);
        for (udi v : g[u]) {
            deg_in_P[v]--;
            deg_in_PC[v]--;
        }
    }

    inline void P_R_remove_if(const Graph& g,
                              const auto& pred)noexcept {
        for (int i = P_R.vnum - 1; i >= 0; i--) {
            if (pred(P_R.vlist[i])) {
                P_R_remove(g, P_R.vlist[i]);
            }
        }
    }
    inline void C_L_add(const Graph& g, udi u)noexcept {
        C_L.add(u);
        // C_bit.set(u);
        // PC_bit.set(u);
        for (udi v : g[u]) {
            deg_in_PC[v]++;
        }
    }

    inline void C_L_remove(const Graph& g, udi u)noexcept {
        C_L.remove(u);
        // C_bit.reset(u);
        // PC_bit.reset(u);
        for (udi v : g[u]) {
            deg_in_PC[v]--;
        }
    }

    inline void C_L_remove_if(const Graph& g,
                              const auto& pred) noexcept{
        for (int i = C_L.vnum - 1; i >= 0; i--) {
            if (pred(C_L.vlist[i])) {
                C_L_remove(g, C_L.vlist[i]);
            }
        }
    }

    inline void C_R_add(const Graph& g, udi u)noexcept {
        C_R.add(u);
        // C_bit.set(u);
        // PC_bit.set(u);
        for (udi v : g[u]) {
            deg_in_PC[v]++;
        }
    }

    inline void C_R_remove(const Graph& g, udi u) noexcept{
        C_R.remove(u);
        // C_bit.reset(u);
        // PC_bit.reset(u);
        for (udi v : g[u]) {
            deg_in_PC[v]--;
        }
    }

    inline void C_R_remove_if(const Graph& g,
                              const auto& pred) noexcept{
        for (int i = C_R.vnum - 1; i >= 0; i--) {
            if (pred(C_R.vlist[i])) {
                C_R_remove(g, C_R.vlist[i]);
            }
        }
    }
    inline void C_R_remove_if_re(const Graph& g,
                              const auto& pred,vector<udi> &remove_CR)noexcept {
        for (int i = C_R.vnum - 1; i >= 0; i--) {
            if (pred(C_R.vlist[i])) {
                remove_CR.emplace_back(C_R.vlist[i]);
                C_R_remove(g, C_R.vlist[i]);
            }
        }
    }
    // inline std::vector<udi> C_R_remove_if_re(const Graph& g,
    //                           std::function<bool(udi)> pred) {
    //     vector<udi> C;
    //     C.reserve(C_R.vnum);
    //     for (int i = C_R.vnum - 1; i >= 0; i--) {
    //         if (pred(C_R.vlist[i])) {
    //             C.push_back(C_R.vlist[i]);
    //             C_R_remove(g, C_R.vlist[i]);
    //         }
    //     }
    //     return C;
    // }
    inline void C_L_remove_if_re(const Graph& g,
                              const auto& pred,vector<udi> &remove_CL) noexcept{
        for (int i = C_L.vnum - 1; i >= 0; i--) {
            if (pred(C_L.vlist[i])) {
                remove_CL.emplace_back(C_L.vlist[i]);
                C_L_remove(g, C_L.vlist[i]);
            }
        }
    }
    inline void X_L_remove_if_re(const Graph& g,
                              const auto& pred,vector<udi> &remove_XL) noexcept{
        for (int i = X_L.vnum - 1; i >= 0; i--) {
            if (pred(X_L.vlist[i])) {
                remove_XL.emplace_back(X_L.vlist[i]);
                X_L_remove(g, X_L.vlist[i]);
            }
        }
    }
    inline void X_R_remove_if_re(const Graph& g,
                              const auto& pred,vector<udi> &remove_XR) noexcept{
        for (int i = X_R.vnum - 1; i >= 0; i--) {
            if (pred(X_R.vlist[i])) {
                remove_XR.emplace_back(X_R.vlist[i]);
                X_R_remove(g, X_R.vlist[i]);
            }
        }
    }
    inline void X_L_add(const Graph& g, udi u) { X_L.add(u); }

    inline void X_L_remove(const Graph& g, udi u) { X_L.remove(u); }

    inline void X_L_remove_if(const Graph& g,
                              const auto& pred) noexcept{
        for (int i = X_L.vnum - 1; i >= 0; i--) {
            if (pred(X_L.vlist[i])) {
                X_L_remove(g, X_L.vlist[i]);
            }
        }
    }

    inline std::vector<udi> X_R_remove_if_re(const Graph& g,
                              const auto& pred) noexcept{
        vector<udi> X;
        X.reserve(C_R.vnum);
        for (int i = X_R.vnum - 1; i >= 0; i--) {
            if (pred(X_R.vlist[i])) {
                X.push_back(X_R.vlist[i]);
                C_R_remove(g, X_R.vlist[i]);
            }
        }
        return X;
    }

    inline void X_R_add(const Graph& g, udi u) noexcept{ X_R.add(u); }

    inline void X_R_remove(const Graph& g, udi u) noexcept{ X_R.remove(u); }

    inline void X_R_remove_if(const Graph& g,
                              const auto& pred) noexcept{
        for (int i = X_R.vnum - 1; i >= 0; i--) {
            if (pred(X_R.vlist[i])) {
                X_R_remove(g, X_R.vlist[i]);
            }
        }
    }

    void reduction(const Graph& g) noexcept{
        auto predl = [&](udi u) {
            return deg_in_P[u] + k < P_R.vnum || deg_in_PC[u] + k < lb_R;
        };
        C_L_remove_if(g, predl);
        X_L_remove_if(g, predl);

        auto predr = [&](udi v) {
            return deg_in_P[v] + k < P_L.vnum || deg_in_PC[v] + k < lb_L;
        };
        C_R_remove_if(g, predr);
        X_R_remove_if(g, predr);
    }
    void reduction(const Graph& g,vector<udi> &remove_CR,vector<udi> &remove_CL,vector<udi> &remove_XR,vector<udi> &remove_XL) noexcept{
        auto predl = [&](udi u) {
            return deg_in_P[u] + k < P_R.vnum || deg_in_PC[u] + k < lb_R;
        };
        C_L_remove_if_re(g, predl,remove_CL);
        X_L_remove_if_re(g, predl,remove_XL);

        auto predr = [&](udi v) {
            return deg_in_P[v] + k < P_L.vnum || deg_in_PC[v] + k < lb_L;
        };
        C_R_remove_if_re(g, predr,remove_CR);
        X_R_remove_if_re(g, predr,remove_XR);
    }
    bool is_need_search(){
    for(udi i=0;i<X_L.vnum;i++){
        if(deg_in_PC[X_L.vlist[i]]==C_R.vnum+P_R.vnum){
            return false;
        }
    }
    for(udi i=0;i<X_R.vnum;i++){
        if(deg_in_PC[X_R.vlist[i]]==C_L.vnum+P_L.vnum){
            return false;
        }
    }
    return true;
}
};
//BiKplex ns;
void new_BK(const Graph& g,udi lvl,BiKplex &ns) {
    if(!ns.exchanged){
        if(ns.P_L.vnum >= global_ubL || ns.P_R.vnum >=global_ubR) return;
    }else{
        if(ns.P_L.vnum >= global_ubR || ns.P_R.vnum >=global_ubL) return;
    }
    if (global_ubL + global_ubR < lb_global()+1) {
        return;
    }
    if(!ns.exchanged){
        if (global_ubL <= ns.lb_L || global_ubR <= ns.lb_R) {
            return;
        }
    }else{
        if (global_ubL <= ns.lb_R || global_ubR <= ns.lb_L) {
            return;
        }
    }
    udi L_size = ns.P_L.vnum + ns.C_L.vnum;
    udi R_size = ns.P_R.vnum + ns.C_R.vnum;
    ns.ub_L=min(ns.ub_L,L_size);
    ns.ub_R=min(ns.ub_R,R_size);
    if (ns.ub_L + ns.ub_R <= lb_global()) {
        return;
    }
    if (ns.ub_L < ns.lb_L || ns.ub_R < ns.lb_R) {
        return;
    }
    if (R_size == ns.lb_R && ns.P_L.vnum>k && !ns.C_R.empty()) {
        for(int i=0;i<ns.C_R.vnum;i++){
            if(ns.P_L.vnum > k + ns.deg_in_P[ns.C_R.vlist[i]]){
                return;
            }
        }
        for(int i=0;i<ns.P_L.vnum;i++){
            if(ns.P_R.vnum > k + ns.deg_in_PC[ns.P_L.vlist[i]]){
                return;
            }
        }
    }
    if (L_size == ns.lb_L && ns.P_R.vnum>k && !ns.C_L.empty() ) {
        for(int i=0;i<ns.C_L.vnum;i++){
            if(ns.P_R.vnum > k + ns.deg_in_P[ns.C_L.vlist[i]]){
                return;
            }
        }
        for(int i=0;i<ns.P_R.vnum;i++){
            if(ns.P_L.vnum > k + ns.deg_in_PC[ns.P_R.vlist[i]]){
                return;
            }
        }
    }
    if (ns.P_L.vnum > k && ns.P_L.vnum <ns.lb_L && ns.P_R.vnum > 0 && (ns.P_R.vnum == ns.ub_R || ns.C_R.empty() || R_size == ns.lb_R)) {
        int num_i=0;
        for(udi u:ns.C_L.vlist){
            if(ns.deg_in_PC[u]==R_size){
                num_i++;
            }
        }
        for(udi u:ns.P_L.vlist){
            if(ns.deg_in_PC[u]==R_size){
                num_i++;
            }
        }
        // int non=0;
        // for(udi u:ns.C_R.vlist){
        //     non+=(k+ns.deg_in_P[u]-ns.P_L.vnum);
        // }
        // for(udi u:ns.P_R.vlist){
        //     non+=(k+ns.deg_in_P[u]-ns.P_L.vnum);
        // }
        // ns.ub_L=min(ns.ub_L,ns.P_L.vnum+non+num_i);
        ns.ub_L=min(ns.ub_L,num_i+k*R_size);
        if (ns.ub_L < ns.lb_L) {
            return;
        }

        int edges = 0;
        if (ns.lb_L + ns.P_R.vnum * ns.P_L.vnum >
            ns.P_L.vnum + edges + ns.P_R.vnum * k) {
            for (int i = 0; i < ns.P_R.vnum; i++) {
                edges += ns.deg_in_P[ns.P_R.vlist[i]];
            }
            if(!ns.C_R.empty() && R_size == ns.lb_R){
                for (int i = 0; i < ns.C_R.vnum; i++) {
                    edges += ns.deg_in_P[ns.C_R.vlist[i]];
                }
            }
            vector<int> con;
            con.resize(g.vertex_tot,0);
            int num=0;
            if (!ns.C_R.empty() && R_size == ns.lb_R){
                for(udi u:ns.C_L.vlist){
                    if(ns.deg_in_PC[u]==R_size){
                        num++;
                    }
                }
                for(udi u:ns.P_L.vlist){
                    if(ns.deg_in_PC[u]==R_size){
                        num++;
                    }
                }
            if (ns.lb_L + R_size * ns.P_L.vnum >
                num + ns.P_L.vnum + edges + R_size * k) {
                        // cout<<"workR ";
                return;
            }
            }else{
                for(udi u:ns.C_L.vlist){
                    if(ns.deg_in_P[u]==ns.P_R.vnum){
                        num++;
                    }
                }
                for(udi u:ns.P_L.vlist){
                    if(ns.deg_in_P[u]==ns.P_R.vnum){
                        num++;
                    }
                }
                if (ns.lb_L + ns.P_R.vnum * ns.P_L.vnum >
                    num + ns.P_L.vnum + edges + ns.P_R.vnum * k) {
                        // cout<<"workR ";
                    return;
                }
            }
        }
    }
    if (ns.P_R.vnum > k && ns.P_R.vnum <ns.lb_R&& ns.P_L.vnum > 0 && (ns.P_L.vnum == ns.ub_L || ns.C_L.empty() || L_size == ns.lb_L)) {
        int num_i=0;
        for(udi u:ns.C_R.vlist){
            if(ns.deg_in_PC[u]==L_size){
                num_i++;
            }
        }
        for(udi u:ns.P_R.vlist){
            if(ns.deg_in_PC[u]==L_size){
                num_i++;
            }
        }
        ns.ub_R=min(ns.ub_R,num_i+k*L_size);
        if (ns.ub_R < ns.lb_R) {
            return;
        }
        int edges = 0;
        if (ns.lb_R + ns.P_L.vnum * ns.P_R.vnum >
            ns.P_R.vnum + edges + ns.P_L.vnum * k) {
        for (int i = 0; i < ns.P_L.vnum; i++) {
            edges += ns.deg_in_P[ns.P_L.vlist[i]];
        }
        if(!ns.C_L.empty() && L_size == ns.lb_L){
            for (int i = 0; i < ns.C_L.vnum; i++) {
                edges += ns.deg_in_P[ns.C_L.vlist[i]];
            }
        }
            int num=0;
            if (!ns.C_L.empty() && L_size == ns.lb_L){
                for(udi u:ns.C_R.vlist){
                    if(ns.deg_in_PC[u]==L_size){
                        num++;
                    }
                }
                for(udi u:ns.P_R.vlist){
                    if(ns.deg_in_PC[u]==L_size){
                        num++;
                    }
                }
                if (ns.lb_R + L_size * ns.P_R.vnum >
                    num + ns.P_R.vnum + edges + L_size * k) {
                    return;
                }
            }else{
                for(udi u:ns.C_R.vlist){
                    if(ns.deg_in_P[u]==ns.P_L.vnum){
                        num++;
                    }
                }
                for(udi u:ns.P_R.vlist){
                    if(ns.deg_in_P[u]==ns.P_L.vnum){
                        num++;
                    }
                }
                if (ns.lb_R + ns.P_L.vnum * ns.P_R.vnum >
                    num + ns.P_R.vnum + edges + ns.P_L.vnum * k) {
                    return;
                }
            }
    }
    }
    L_size = ns.P_L.vnum + ns.C_L.vnum;
    R_size = ns.P_R.vnum + ns.C_R.vnum;
    if (ns.C_L.empty() && ns.C_R.empty()) {
        if (ns.X_L.empty() && ns.X_R.empty()) {
            SubGraph sg = ns.get_subgraph();
            g.get_original_vertex(sg);
            record(sg,g,ns.exchanged);
        }
        return;
    }
    udi min_deg = g.vertex_tot;
    udi pivot = MAX;
    udi max_ndeg = 0;
    udi min_d_n = MAX;
    for(udi u:ns.P_L.vlist){
        if(R_size - ns.deg_in_PC[u]>k && R_size-ns.deg_in_PC[u] > max_ndeg){
                pivot = u;
                max_ndeg=R_size-ns.deg_in_PC[u];
        }
        if (ns.deg_in_PC[u] < min_deg) {
            min_deg = ns.deg_in_PC[u];
            min_d_n = u;
        }
    }
        if (min_deg + k < ns.lb_R) {
            return;
        } 
        if(min_deg + k < ns.ub_R)
        ns.ub_R = min_deg + k;
    if (lb_global() + 1 > ns.ub_R + ns.lb_L) ns.lb_L = lb_global() + 1 - ns.ub_R;
    udi min_deg_L=min_deg;
    udi min_d_n_L=min_d_n;
    min_deg = g.vertex_tot;
        for (udi v : ns.P_R.vlist) {
            if(L_size - ns.deg_in_PC[v]>k && L_size-ns.deg_in_PC[v] > max_ndeg){
                    pivot = v;
                    max_ndeg = L_size-ns.deg_in_PC[v];
            }
            if (ns.deg_in_PC[v] < min_deg) {
                min_deg = ns.deg_in_PC[v];
                min_d_n = v;
            }
        }
        if (min_deg + k < ns.lb_L) {
            return;
        }
        if(min_deg + k < ns.ub_L)
        ns.ub_L = min_deg + k;
        if (lb_global() + 1 > ns.ub_L + ns.lb_R) ns.lb_R = lb_global() + 1 - ns.ub_L;
    L_size = ns.P_L.vnum + ns.C_L.vnum;
    R_size = ns.P_R.vnum + ns.C_R.vnum;
    if (ns.ub_L + ns.ub_R <= lb_global()) {
        return;
    }
    if (ns.ub_L < ns.lb_L || ns.ub_R < ns.lb_R) {
        return;
    }
    if (L_size + R_size <= lb_global()) {
        return;
    }
    if (L_size < ns.lb_L || R_size < ns.lb_R) {
        return;
    }
    if (pivot != MAX && max_ndeg>k) {
        if(L_size+R_size<=lb_global()+1){
            return;
        }
        vector<udi> remove_CL,remove_CR,remove_XL,remove_XR;
        remove_CL.reserve(L_size);
        remove_CR.reserve(R_size);
        remove_XL.reserve(L_size);
        remove_XR.reserve(R_size);
        if (!g.is_L(pivot, ns.exchanged)) {
            ns.exchange(remove_CL,remove_CR,remove_XL,remove_XR);
        }
    if(lvl<=ns.olvl_l){
        ns.flag_l2=false;
        ns.olvl_l=MAX;
    }
    if(min_deg_L + k==ns.lb_R && !ns.flag_l2){
        ns.olvl_l=lvl;
        auto pred=[&](udi u){
            int cn=0;
            for(udi v:g[min_d_n_L]){
                if(g.is_neighbor(u,v)){
                    cn++;
                }
            }
            if(k+cn+ns.deg_in_P[min_d_n_L] + ns.deg_in_P[u]< ns.P_R.vnum + ns.deg_in_PC[min_d_n_L]){
                return true;
            }
            return false;
        };
        //ns.C_L_remove_if(g,pred);
        ns.C_L_remove_if_re(g,pred,remove_CL);
        ns.flag_l2=true;
    }
    if(lvl<=ns.olvl_r){
        ns.flag_r2=false;
        ns.olvl_r=MAX;
    }
    if(min_deg + k==ns.lb_L && !ns.flag_r2){
        ns.olvl_r=lvl;
        auto pred=[&](udi u){
            int cn=0;
            for(udi v:g[min_d_n]){
                if(g.is_neighbor(u,v)){
                    cn++;
                }
            }
            if(k+cn+ns.deg_in_P[min_d_n] + ns.deg_in_P[u]< ns.P_L.vnum + ns.deg_in_PC[min_d_n]){
                return true;
            }
            return false;
        };
        ns.C_R_remove_if_re(g,pred,remove_CR);
        ns.flag_r2=true;
    }
        vector<udi> doing;
        doing.reserve(max_ndeg);
        for (udi v : ns.C_R.vlist) {
            if (!g.is_neighbor(pivot, v)) {
                doing.emplace_back(v);
            }
        }
        udi bunch_num = ns.deg_in_P[pivot] + k - ns.P_R.vnum;
        udi idx = 0;
        for (idx = 0; idx < doing.size(); idx++) {
            udi v = doing[idx];
            if (!ns.C_R.contains(v)) {
                continue;
            }
            ns.C_R_remove(g, v);
            ns.X_R_add(g, v);
            udi  ubl=ns.ub_L;udi ubr =ns.ub_R;
            udi  lbl=ns.lb_L;udi lbr =ns.lb_R;
            if(ns.is_need_search() && ns.P_R.vnum + ns.C_R.vnum >= ns.lb_R&& ns.P_L.vnum + ns.C_L.vnum >= ns.lb_L){
            lvl++;
            new_BK(g,lvl,ns);
            lvl--;
            }
            if (!g.is_L(pivot, ns.exchanged)) {
                ns.exchange();
            }
            ns.ub_L=ubl;ns.ub_R=ubr;
            ns.lb_L=lbl;ns.lb_R=lbr;
            if (lb_global() + 1 > ns.ub_L + ns.lb_R) ns.lb_R = lb_global() + 1 - ns.ub_L;
            if (lb_global() + 1 > ns.ub_R + ns.lb_L) ns.lb_L = lb_global() + 1 - ns.ub_R;
            auto pred=[&](udi u){
                udi num = 0;
                for (udi x : g[v]) {
                    if ((g.is_neighbor(u, x)&& (ns.C_L.contains(x) || ns.P_L.contains(x)))) {
                        num++;
                    }
                }
                if (num + 2 * k < ns.lb_L) return true;
                return false;
            };
            ns.C_R_remove_if_re(g,pred,remove_CR);
            ns.X_R_remove(g, v);
            ns.P_R_add(g, v,remove_CR,remove_CL,remove_XR,remove_XL);
            ns.reduction(g,remove_CR,remove_CL,remove_XR,remove_XL);
            if(ns.P_R.vnum + ns.C_R.vnum <ns.lb_R || ns.P_L.vnum + ns.C_L.vnum <ns.lb_L){
                break;
            }
            bunch_num--;
            if (bunch_num == 0) {
                break;
            }
        }

        for (; idx < doing.size(); idx++) {
            udi v = doing[idx];
            if (ns.C_R.contains(v)) {
                ns.C_R_remove(g, v);
            }
        }
        udi  ubl=ns.ub_L;udi ubr =ns.ub_R;
        udi  lbl=ns.lb_L;udi lbr =ns.lb_R;
        if(ns.is_need_search() && ns.P_R.vnum + ns.C_R.vnum >= ns.lb_R&& ns.P_L.vnum + ns.C_L.vnum >= ns.lb_L){
        lvl++;
        new_BK(g,lvl,ns);
        lvl--;
        }
        if (!g.is_L(pivot, ns.exchanged)) {
            ns.exchange();
        }
        ns.ub_L=ubl;ns.ub_R=ubr;
        ns.lb_L=lbl;ns.lb_R=lbr;
        if (lb_global() + 1 > ns.ub_L + ns.lb_R) ns.lb_R = lb_global() + 1 - ns.ub_L;
        if (lb_global() + 1 > ns.ub_R + ns.lb_L) ns.lb_L = lb_global() + 1 - ns.ub_R;
        for(udi u:doing){
            if(ns.P_R.contains(u)){
                // ns.PR2CR(g,u);
                ns.P_R_remove(g,u);
                ns.C_R_add(g,u);
            }
            if(!ns.C_R.contains(u)){
                ns.C_R_add(g,u);
            }
        }
        for(udi u:remove_CR){
            if(!ns.C_R.contains(u))
            ns.C_R_add(g,u);
        }
        for(udi u:remove_CL){
            if(!ns.C_L.contains(u))
            ns.C_L_add(g,u);
        }
        for(udi u:remove_XR){
            if(!ns.X_R.contains(u))
            ns.X_R_add(g,u);
        }
        for(udi u:remove_XL){
            if(!ns.X_L.contains(u))
            ns.X_L_add(g,u);
        }
    } else {
        pivot=MAX;
        min_d_n=MAX;
        min_d_n_L=MAX;
        min_deg=MAX;
        min_deg_L=MAX;
        for (udi u : ns.C_L.vlist) {
            if(R_size - ns.deg_in_PC[u] > k && R_size - ns.deg_in_PC[u] > max_ndeg){
                max_ndeg=R_size-ns.deg_in_PC[u];
                pivot = u;
            }
            if (ns.deg_in_PC[u] < min_deg) {
                min_deg = ns.deg_in_PC[u];
                min_d_n = u;
            }
        }
        min_d_n_L=min_d_n;
        min_deg_L=min_deg;
        for (udi u : ns.C_R.vlist) {
            if (L_size - ns.deg_in_PC[u] > k && L_size - ns.deg_in_PC[u] > max_ndeg) {
                    max_ndeg=L_size-ns.deg_in_PC[u];
                    pivot = u;
            }
            if (ns.deg_in_PC[u] < min_deg) {
                min_deg = ns.deg_in_PC[u];
                min_d_n = u;
            }
        }
        if (max_ndeg <= k || pivot==MAX) {
            if(ns.C_L.vnum+ns.P_L.vnum > ns.ub_L || ns.C_R.vnum+ns.P_R.vnum >ns.ub_R) return;
            if(!ns.exchanged){
                if(ns.C_L.vnum+ns.P_L.vnum >= global_ubL || ns.C_R.vnum+ns.P_R.vnum >=global_ubR) return;
            }else{
                if(ns.C_L.vnum+ns.P_L.vnum >= global_ubR || ns.C_R.vnum+ns.P_R.vnum >=global_ubL) return;
            }
            if (ns.X_L.empty() && ns.X_R.empty()) {
                    SubGraph sg = ns.get_subgraph_all();
                    g.get_original_vertex(sg);
                    record(sg,g,ns.exchanged);
                    return;
            }
            if(!ns.is_need_search()) return;
            vector<udi> lim_L,lim_R;
            lim_L.reserve(ns.P_L.vnum);
            lim_R.reserve(ns.P_R.vnum);
            for(udi u:ns.P_L.vlist){
                if(ns.deg_in_PC[u]==R_size - k){
                    lim_L.emplace_back(u);
                }
            }
            for(udi u:ns.C_L.vlist){
                if(ns.deg_in_PC[u]==R_size - k){
                    lim_L.emplace_back(u);
                }
            }
            for(udi u:ns.X_R.vlist){
                if(ns.deg_in_PC[u]>=L_size - k){
                bool flag=true;
                    for(udi v:lim_L){
                        if(!g.is_neighbor(u,v)){
                            flag=false;
                            break;
                        }
                    }
                    if(flag) return;
                }
            }
            for(udi u:ns.P_R.vlist){
                if(ns.deg_in_PC[u]==L_size - k){
                    lim_R.emplace_back(u);
                }
            }
            for(udi u:ns.C_R.vlist){
                if(ns.deg_in_PC[u]==L_size - k){
                    lim_R.emplace_back(u);
                }
            }
            for(udi u:ns.X_L.vlist){
                if(ns.deg_in_PC[u]>=R_size - k){
                    bool flag=true;
                    for(udi v:lim_R){
                        if(!g.is_neighbor(u,v)){
                            flag=false;
                            break;
                        }
                    }
                    if(flag) return;        
                }
            }
            SubGraph sg = ns.get_subgraph_all();
            g.get_original_vertex(sg);
            record(sg,g,ns.exchanged);
            return;
        }
        if(L_size+R_size<=lb_global()+1){
            return;
        }
        if (!g.is_L(pivot, ns.exchanged)) {
            ns.exchange();
        }
        ns.C_L_remove(g, pivot);
        ns.X_L_add(g, pivot);
        udi  ubl=ns.ub_L;udi ubr =ns.ub_R;
        udi  lbl=ns.lb_L;udi lbr =ns.lb_R;
        if(ns.is_need_search() && ns.P_R.vnum + ns.C_R.vnum >= ns.lb_R&& ns.P_L.vnum + ns.C_L.vnum >= ns.lb_L){
        lvl++;
        new_BK(g,lvl,ns);
        lvl--;
        }
        //ns.exchanged=O_E;
        if (!g.is_L(pivot, ns.exchanged)) {
            ns.exchange();
        }
        ns.ub_L=ubl;ns.ub_R=ubr;
        ns.lb_L=lbl;ns.lb_R=lbr;
        if (lb_global() + 1 > ns.ub_L + ns.lb_R) ns.lb_R = lb_global() + 1 - ns.ub_L;
        if (lb_global() + 1 > ns.ub_R + ns.lb_L) ns.lb_L = lb_global() + 1 - ns.ub_R;
        vector<udi> remove_CL,remove_CR,remove_XL,remove_XR;
        remove_CL.reserve(L_size);
        remove_CR.reserve(R_size);
        remove_XL.reserve(L_size);
        remove_XR.reserve(R_size);
        ns.X_L_remove(g, pivot);
        ns.P_L_add(g, pivot,remove_CR,remove_CL,remove_XR,remove_XL);
        ns.reduction(g,remove_CR,remove_CL,remove_XR,remove_XL);
        ubl=ns.ub_L;ubr =ns.ub_R;
        lbl=ns.ub_L;lbr =ns.ub_R;
        if(ns.is_need_search() && ns.P_R.vnum + ns.C_R.vnum >= ns.lb_R && ns.P_L.vnum + ns.C_L.vnum >= ns.lb_L){
        lvl++;
        new_BK(g,lvl,ns);
        lvl--;
        }
        //ns.exchanged=O_E;
        if (!g.is_L(pivot, ns.exchanged)) {
            ns.exchange();
        }
        ns.ub_L=ubl;ns.ub_R=ubr;
        ns.lb_L=lbl;ns.lb_R=lbr;
        if (lb_global() + 1 > ns.ub_L + ns.lb_R) ns.lb_R = lb_global() + 1 - ns.ub_L;
        if (lb_global() + 1 > ns.ub_R + ns.lb_L) ns.lb_L = lb_global() + 1 - ns.ub_R;
        for(udi u:remove_CR){
            ns.C_R_add(g,u);
        }
        for(udi u:remove_CL){
            ns.C_L_add(g,u);
        }
        for(udi u:remove_XR){
            ns.X_R_add(g,u);
        }
        for(udi u:remove_XL){
            ns.X_L_add(g,u);
        }
        ns.P_L_remove(g,pivot);
        ns.C_L_add(g,pivot);
        // ns.PL2CL(g,pivot);
    }
}

void EnOneHop(const Graph& g,const BiKplex& os,BiKplex& ns,vector<udi>&subset,udi start,vector<udi> &all,udi k,udi num){
    if(subset.size()<=num){
            ns=os;

            bool flag=false;
            for(udi v:subset){
            int nei=0;
            for(udi u:ns.P_L.vlist){
                for(udi x:ns.C_R.vlist){
                    if(g.is_neighbor(x,u) && g.is_neighbor(x,v)){
                        nei++;
                    }
                }
            }
            if(nei<ns.lb_L-2*k){
                flag=true;
                break;
            }
            ns.CL2PL(g,v);

            }
            for(udi u:ns.C_R.vlist){
                if(ns.deg_in_P[u] + k < ns.P_L.vnum){
                    return;
                }
            }
            if(ns.C_R.vnum + ns.P_R.vnum<ns.lb_R || ns.C_L.vnum + ns.P_L.vnum<ns.lb_L) return;
            if(subset.size()<=num && !flag)
            for(udi v:all){
                if(ns.C_L.contains(v)){
                    ns.C_L_remove(g,v);
                    if(ns.deg_in_PC[v] + k >= ns.lb_R)
                    ns.X_L_add(g,v);
                    ns.X_L_add(g,v);
                }
            }
            if(!flag)
            ns.reduction(g);
            if(!flag &&ns.is_need_search()){
                if(ns.C_R.vnum + ns.P_R.vnum>=ns.lb_R && ns.C_L.vnum + ns.P_L.vnum>=ns.lb_L){
                    new_BK(g,0,ns);
                }
            }
    }
    for(udi i=start;i<all.size() && subset.size()<num;i++){
        subset.emplace_back(all[i]);
        EnOneHop(g,os,ns,subset,i+1,all,k,num);
        subset.pop_back();
    }
}
void generateSubsetsForL(const Graph& g,const BiKplex& os,BiKplex& ns,vector<udi>&subset,udi start,vector<udi> &all,udi k){
    if(subset.size()<=k){
            ns=os;
            bool flag=false;
            for(udi v:subset){
            for(udi u:ns.P_R.vlist){
                            int nei=0;
                for(udi x:ns.C_L.vlist){
                    if(g.is_neighbor(x,u) && g.is_neighbor(x,v)){
                        nei++;
                    }
                }
            if(nei<ns.lb_R-2*k){
                flag=true;
                break;
            }
            }
            if(flag) break;
            ns.P_R_add(g, v);
            }
            if(subset.size()<=k && !flag)
            for(udi v:all){
                if(!ns.P_R.contains(v)){
                    if(ns.deg_in_PC[v] + k >= ns.lb_L)
                    ns.X_R_add(g,v);
                }
            }
            if(ns.lb_L==ns.C_L.vnum+ns.P_L.vnum && !flag){
                for(udi u:ns.C_L.vlist){
                    if(ns.deg_in_P[u] + k < ns.P_R.vnum){
                        return;
                    }
                }
            }
            if(!flag)
            ns.reduction(g);
            if(!flag && ns.is_need_search() ){
                if(ns.C_R.vnum + ns.P_R.vnum>=ns.lb_R && ns.C_L.vnum + ns.P_L.vnum>=ns.lb_L){
                    new_BK(g,0,ns);
                }
            }
    }
    for(udi i=start;i<all.size() && subset.size()<k;i++){
        subset.emplace_back(all[i]);
        generateSubsetsForL(g,os,ns,subset,i+1,all,k);
        subset.pop_back();
    }

}
void generateSubsets(const Graph& g,const BiKplex& os,BiKplex& ns,vector<udi>&subset,udi start,vector<udi> &all,udi k){
    if(subset.size()<=k){
            ns=os;
            bool flag=false;
            for(udi v:subset){
            for(udi u:ns.P_L.vlist){
                            int nei=0;
                for(udi x:ns.C_R.vlist){
                    if(g.is_neighbor(x,u) && g.is_neighbor(x,v)){
                        nei++;
                    }
                }
            if(nei<ns.lb_L-2*k){
                flag=true;
                break;
            }
            }
            if(flag) break;
            ns.P_L_add(g, v);
            }
            if(subset.size()<=k && !flag)
            for(udi v:all){
                if(!ns.P_L.contains(v)){
                    if(ns.deg_in_PC[v] + k >= ns.lb_R)
                    ns.X_L_add(g,v);
                }
            }
            if(ns.lb_R==ns.C_R.vnum+ns.P_R.vnum && !flag){
                for(udi u:ns.C_R.vlist){
                    if(ns.deg_in_P[u] + k < ns.P_L.vnum){
                        return;
                    }
                }
            }
            if(!flag)
            ns.reduction(g);
            if(!flag && ns.is_need_search() ){
                if(ns.C_R.vnum + ns.P_R.vnum>=ns.lb_R && ns.C_L.vnum + ns.P_L.vnum>=ns.lb_L){
                    if(ns.P_L.vnum<k){
                        new_BK(g,0,ns);
                    }else{
                        new_BK(g,0,ns);
                    }

                }
            }
    }
    for(udi i=start;i<all.size() && subset.size()<k;i++){
        subset.emplace_back(all[i]);
        generateSubsets(g,os,ns,subset,i+1,all,k);
        subset.pop_back();
    }

}
void decomposition(const Graph& g,BiKplex& s,RandList & delete_v) {
    udi num_L = 0, num_R = 0;
    num_e = 0;
    int flag=0;
    int skip_num=0;
    int max_num=0,max_num_r=0;
    if(s.exchanged){
        s.exchange();
    }
    global_ubL=s.ub_L+1;
    global_ubR=s.ub_R+1;
    for (int i=0;i<g.vertex_tot;i++) {
        if(delete_v.contains(i)) {
            continue;
        }
        if(g.is_L(i,s.exchanged)){
            if(g[i].size()<s.lb_R-k){
                continue;
            }
        }else{
            if(g[i].size()<s.lb_L-k){
                continue;
            }
        }
        int ui=i;
        current_pivot=ui;
        max_lvl=0;
        BiKplex ns=s;
        if (g.is_L(ui, s.exchanged)) {
            ns.P_L_add(g,ui);
            for(udi v:g[ui]){
                    int count=0;
                    for(udi x:g[v]){
                        if(x>=ui){
                            count++;
                        }
                    }
                    if(count<ns.lb_L - k){
                        continue;
                    }
                    if(v>ui){
                        ns.C_R_add(g,v);
                    }else if(v<ui){
                        ns.X_R_add(g,v);
                    }
            }
            vector<udi> con_ui(g.vertex_tot,0);
            for(udi v:ns.C_R.vlist){
                for(udi x:g[v]){
                    if(ns.C_L.contains(x) || ns.X_L.contains(x)) continue;
                    if(x>ui && ns.deg_in_PC[x]>=ns.lb_R-2*k){
                        con_ui[x]=ns.deg_in_PC[x];
                        ns.C_L_add(g,x);
                    }else if(x<ui && ns.deg_in_PC[x]>=ns.lb_R-2*k){
                        ns.X_L_add(g,x);
                    }
                }
            }
            for(int i=0;i<ns.C_L.vnum;i++){
                for(udi v:g[ns.C_L.vlist[i]]){
                    if(v>ui){
                        if(ns.C_R.contains(v)){
                            if(ns.deg_in_PC[v]<ns.lb_L-k){
                                ns.C_R_remove(g,v);
                                for(udi r:g[v]){
                                    con_ui[r]--;
                                }
                            }
                        }else{
                            if(ns.deg_in_PC[v]>=ns.lb_L-k){
                                ns.C_R_add(g,v);
                            }
                        }
                    }else if(v<ui){
                        if(ns.X_R.contains(v)){
                            if(ns.deg_in_PC[v]<ns.lb_L-k){
                                ns.X_R_remove(g,v);
                            }
                        }else{
                            if(ns.deg_in_PC[v]>=ns.lb_L-k){
                                ns.X_R_add(g,v);
                            }
                        }
                    }
                }
            }
            udi num;
                RandList three_hop;
                three_hop.reinit(ns.C_R.vnum);
                RandList two_hop;
                two_hop.reinit(ns.C_L.vnum);
                int num_de=0;
            while(true){
                num=ns.C_L.vnum+ns.C_R.vnum+ns.X_L.vnum+ns.X_R.vnum;
                for(int i=ns.C_L.vnum-1;i>=0;i--){
                    if(con_ui[ns.C_L.vlist[i]]<ns.lb_R-2*k || ns.deg_in_PC[ns.C_L.vlist[i]]<ns.lb_R - k){
                        if(con_ui[ns.C_L.vlist[i]]<ns.lb_R-2*k) num_de++;
                        ns.C_L_remove(g,ns.C_L.vlist[i]);
                    }
                }
                for(int i=ns.C_R.vnum-1;i>=0;i--){
                    if(ns.deg_in_PC[ns.C_R.vlist[i]]<ns.lb_L - k){
                        udi v=ns.C_R.vlist[i];
                        if(g.is_neighbor(ui,v)){
                                for(udi r:g[v]){
                                    con_ui[r]--;
                                }
                        }
                        ns.C_R_remove(g,v);
                    }
                }
                auto predl = [&] (udi u){
                    if(ns.deg_in_PC[u]<ns.lb_R - k){
                        return true;
                    }else{
                        return false;
                    }
                };
                ns.X_L_remove_if(g,predl);
                auto predr = [&] (udi u){
                    if(ns.deg_in_PC[u]<ns.lb_L - k){
                        return true;
                    }else{
                        return false;
                    }
                };
                ns.X_R_remove_if(g,predr);
                if(num==ns.C_L.vnum+ns.C_R.vnum+ns.X_L.vnum+ns.X_R.vnum){
                    break;
                }
            }
            bool is_search=true;
            for(int i=0;i<ns.X_L.vnum;i++){
                if(ns.deg_in_PC[ns.X_L.vlist[i]]==ns.C_R.vnum){
                    is_search=false;
                    break;
                }
            }
            for(int i=0;i<ns.X_R.vnum;i++){
                if(ns.deg_in_PC[ns.X_R.vlist[i]]==ns.C_L.vnum+1){
                    is_search=false;
                    break;
                }
            }
            ns.ub_L=ns.ub_L > (ns.C_L.vnum+1) ? (ns.C_L.vnum+1):ns.ub_L;
            ns.ub_R=ns.ub_R > (ns.C_R.vnum) ? (ns.C_R.vnum):ns.ub_R;
            vector<udi> nei3;
            int one=0;
            nei3.reserve(ns.C_R.vnum);
            for(int i=ns.C_R.vnum-1;i>=0;i--){
                if(g.is_neighbor(ui,ns.C_R.vlist[i])){
                    one++;
                }else{
                    nei3.emplace_back(ns.C_R.vlist[i]);
                    ns.C_R_remove(g,ns.C_R.vlist[i]);
                }
            }
             if(ns.C_R.vnum + k < ns.lb_R || ns.C_L.vnum + 1  < ns.lb_L) continue;
            if(is_search){
                    if(ns.lb_L==ns.C_L.vnum+ns.P_L.vnum){
                    int num_i=0;
                    for(udi u:ns.C_R.vlist){
                        if(ns.deg_in_PC[u]==ns.C_L.vnum+ns.P_L.vnum){
                            num_i++;
                        }
                    }
                    if(num_i + ns.lb_L*k + k< ns.lb_R) {
                         continue;
                    }else{
                        ns.ub_R=min(ns.ub_R,num_i + ns.lb_L*k+k);
                    }
                }
                    for(udi u:nei3){
                        ns.C_R_add(g,u);
                    }
                new_BK(g,0,ns);
            }
        } else {
            ns.P_R_add(g,ui);
            for(udi v:g[ui]){
                    int count=0;
                    for(udi x:g[v]){
                        if(x>=ui){
                            count++;
                        }
                    }
                    if(count<ns.lb_R - k){
                        continue;
                    }
                    if(v>ui){
                        ns.C_L_add(g,v);
                    }else if(v<ui){
                        ns.X_L_add(g,v);
                    }
            }
            vector<udi> con_ui(g.vertex_tot,0);
            for(udi v:ns.C_L.vlist){
                for(udi x:g[v]){
                    if(ns.C_R.contains(x) || ns.X_R.contains(x)) continue;
                    if(x>ui && ns.deg_in_PC[x]>=ns.lb_L-2*k){
                        con_ui[x]=ns.deg_in_PC[x];
                        ns.C_R_add(g,x);
                    }else if(x<ui &&ns.deg_in_PC[x]>=ns.lb_L-2*k){
                        ns.X_R_add(g,x);
                    }
                }
            }
            for(int i=0;i<ns.C_R.vnum;i++){
                for(udi v:g[ns.C_R.vlist[i]]){
                    if(v>ui){
                        if(ns.C_L.contains(v)){
                            if(ns.deg_in_PC[v]<ns.lb_R-k){
                                ns.C_L_remove(g,v);
                                for(udi u:g[v]){
                                    con_ui[u]--;
                                }
                            }
                        }else{
                            if(ns.deg_in_PC[v]>=ns.lb_R-k){
                                ns.C_L_add(g,v);
                            }
                        }
                    }else if(v<ui){
                        if(ns.X_L.contains(v)){
                            if(ns.deg_in_PC[v]<ns.lb_R-k){
                                ns.X_L_remove(g,v);
                            }
                        }else{
                            if(ns.deg_in_PC[v]>=ns.lb_R-k){
                                ns.X_L_add(g,v);
                            }
                        }
                    }
                }
            }
            udi num;
                RandList three_hop;
                three_hop.reinit(ns.C_R.vnum);
                RandList two_hop;
                two_hop.reinit(ns.C_L.vnum);
            while(true){
                num=ns.C_L.vnum+ns.C_R.vnum+ns.X_L.vnum+ns.X_R.vnum;
                for(int i=ns.C_R.vnum-1;i>=0;i--){
                    if(ns.deg_in_PC[ns.C_R.vlist[i]]<ns.lb_L - k || con_ui[ns.C_R.vlist[i]]<ns.lb_L-2*k){
                        ns.C_R_remove(g,ns.C_R.vlist[i]);
                    }
                }
                for(int i=ns.C_L.vnum-1;i>=0;i--){
                    if(ns.deg_in_PC[ns.C_L.vlist[i]]<ns.lb_R - k){
                        udi v=ns.C_L.vlist[i];
                        ns.C_L_remove(g,v);
                        if(g.is_neighbor(ui,v)){
                            for(udi r:g[v]){
                                con_ui[r]--;
                            }
                        }
                    } 
                }
                auto predr = [&] (udi u){
                    if(ns.deg_in_PC[u]<ns.lb_L - k){
                        return true;
                    }else{
                        return false;
                    }
                };
                ns.X_R_remove_if(g,predr);
                auto predl = [&] (udi u){
                    if(ns.deg_in_PC[u]<ns.lb_R - k){
                        return true;
                    }else{
                        return false;
                    }
                };
                ns.X_L_remove_if(g,predl);
                if(num==ns.C_L.vnum+ns.C_R.vnum+ns.X_L.vnum+ns.X_R.vnum){
                    break; 
                }
            }
            bool is_search=true;
            for(int i=0;i<ns.X_R.vnum;i++){
                if(ns.deg_in_PC[ns.X_R.vlist[i]]==ns.C_L.vnum){
                    is_search=false;
                    break;
                }
            }
            for(int i=0;i<ns.X_L.vnum;i++){
                if(ns.deg_in_PC[ns.X_L.vlist[i]]==ns.C_R.vnum+1){
                    is_search=false;
                    break;
                }
            }
            int one=0;
            vector<udi> nei3;
            bool flag=false;
            nei3.reserve(ns.C_L.vnum);
            for(int i=ns.C_L.vnum-1;i>=0;i--){
                if(ns.deg_in_PC[ns.C_L.vlist[i]] + k ==ns.lb_R){
                    flag=true;
                }
                if(g.is_neighbor(ui,ns.C_L.vlist[i])){
                    one++;
                }else{
                    nei3.emplace_back(ns.C_L.vlist[i]);
                    ns.C_L_remove(g,ns.C_L.vlist[i]);
                }
            }
            if(ns.C_R.vnum + 1 < ns.lb_R || ns.C_L.vnum + k < ns.lb_L) continue;
            ns.ub_L=ns.ub_L > (ns.C_L.vnum+ nei3.size()) ? (ns.C_L.vnum+ nei3.size()):ns.ub_L;
            ns.ub_R=ns.ub_R > (ns.C_R.vnum+1) ? (ns.C_R.vnum+1):ns.ub_R;
            if(is_search && nei3.size()>10*(ns.C_L.vnum+ns.C_R.vnum)){
                vector<udi> subset;
                subset.reserve(k+1);
                if(ns.lb_R==ns.C_R.vnum+ns.P_R.vnum){
                    for(udi u:ns.C_L.vlist){
                        if(ns.deg_in_PC[u]==ns.C_R.vnum+ns.P_R.vnum){
                            num++;
                        }
                    }
                    if(num + ns.lb_R*k + k< ns.lb_L) {
                         continue;
                    }else{
                        ns.ub_L=min(ns.ub_L,num + ns.lb_R*k+k);
                    }
                }
                BiKplex os;
                os.reInit(ns);
                generateSubsets(g,os,ns,subset,0,nei3,k);
            }else if(is_search){
                int num=0;
                    if(ns.lb_R==ns.C_R.vnum+ns.P_R.vnum){
                    for(udi u:ns.C_L.vlist){
                        if(ns.deg_in_PC[u]==ns.C_R.vnum+ns.P_R.vnum){
                            num++;
                        }
                    }
                    if(num + ns.lb_R*k + k< ns.lb_L) {
                         continue;
                    }else{
                        ns.ub_L=min(ns.ub_L,num + ns.lb_R*k+k);
                    }
                }
                for(udi u:nei3){
                    ns.C_L_add(g,u);
                }
                    new_BK(g,0,ns);
            }
        }
    }
}
void IB(const Graph& g, udi theta_L, udi theta_R) {
    udi _lb_L = theta_L, _ub_L = k;
    udi _lb_R = theta_R, _ub_R = k;
    vector<int> l;
    vector<int> r;
    int min_l, min_r;
    for (udi i = 0; i < g.vertex_tot; i++) {
        if (g.in_L[i]) {
            l.push_back(g[i].size());
        } else {
            r.push_back(g[i].size());
        }
    }
    nth_element(r.begin(),r.begin()+theta_R-1, r.end(),greater<int>());
    _ub_L = r[theta_R - 1] + k + 1;

    _lb_L = _ub_L;
    _lb_R = _ub_R;
    _lb_L = _ub_L >> 3;
    if (_lb_L <= theta_L * 2) _lb_L = _ub_L;
    bool move = false;
    bool is_run=false;
    while (true) {
        if (pq_k.size() < top_k) {
            if(is_run){
                _ub_L=_lb_L;
            }
            _lb_L = max(theta_L, _lb_L / 2);
        } else {
            if(lb_global()+1 > 2*_ub_L){
                break;
            }else{
                _ub_L=_lb_L;
                _lb_L = max(theta_L, _lb_L / 2);
            }
        }
        if (lb_global() + 1> _ub_L) {
            _lb_R = max(theta_R, lb_global() + 1 - _ub_L);
        } else {
            _lb_R = theta_R;
        }
        if (_lb_L == _ub_L) {
            break;
        }
        RandList delete_v;
        bool exchange=false;
        Graph ng = g._2_hop_degeneracy_ordering_and_reduction(_lb_L, _lb_R, k,delete_v,exchange);
            nth_element(l.begin(),l.begin()+_lb_L-1, l.end(),greater<int>());
            _ub_R = l[_lb_L - 1] + k + 1;
             if(_ub_R<_lb_R) continue;
        cout << "_lb_L:" << _lb_L << ",_ub_L" << _ub_L <<",_lb_R:" << _lb_R
             << ",_ub_R" << _ub_R << ",graph_size:" << ng.vertex_tot << endl;
        BiKplex os(ng, _lb_L, _lb_R, _ub_L-1 , _ub_R-1);
        os.exchanged=exchange;
        decomposition(ng, os,delete_v);
        is_run=true;
        delete_v.clear();
    }
    if (pq_k.size() == top_k && _lb_L !=theta_L) {
        _ub_L = _lb_L;
        _lb_L = theta_L;
        if (lb_global() + 1 > _ub_L) {
            _lb_R = max(theta_R, lb_global() + 1 - _ub_L);
        } else {
            _lb_R = theta_R;
        }
        RandList delete_v;
        nth_element(l.begin(),l.begin()+_lb_L-1, l.end(),greater<int>());
        _ub_R = l[_lb_L - 1] + k + 1;
        if(_ub_R<_lb_R) {
            cout<<"_ub_R:"<<_ub_R<<",_lb_R:"<<_lb_R<<endl;
            cout<<"skip"<<endl;
            return;
        }
        bool exchange=false;
        Graph ng = g._2_hop_degeneracy_ordering_and_reduction(_lb_L, _lb_R, k,delete_v,exchange);
        cout << "_lb_L:" << _lb_L << ",_ub_L" << _ub_L << ",_lb_R:" << _lb_R
             << ",_ub_R" << _ub_R << ",graph_size:" << ng.vertex_tot << endl;
        BiKplex os(ng, _lb_L, _lb_R, _ub_L-1 , _ub_R-1 );
        os.exchanged=exchange;
        decomposition(ng, os,delete_v);
        delete_v.clear();
    }
}
};  // namespace Framework