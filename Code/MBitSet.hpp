#pragma once
#include <algorithm>
#include <bit>
#include <vector>
//using std::udi;
using std::vector;

class MBitSet {
   public:
    udi n;
    vector<udi> vec;
    const udi B = sizeof(udi) * 8;

    MBitSet(udi _n = 0) {
        n = _n;
        vec.resize((n + B - 1) / B);
    }

    void operator=(const MBitSet& oth) {
        n = oth.n;
        vec = oth.vec;
    }

    inline void reset_all() { std::fill(vec.begin(), vec.end(), 0); }

    inline void flip() {
        for (udi& v : vec) {
            v = ~v;
        }
    }

    inline void set(udi x) { vec[x / B] |= 1ull << (x % B); }

    inline void reset(udi x) { vec[x / B] &= ~(1ull << (x % B)); }

    inline bool get(udi x) const { return vec[x / B] & (1ull << (x % B)); }

    inline bool operator[](udi x) const { return vec[x / B] & (1ull << (x % B)); }

    inline udi lowbit() const {
        for (udi i = 0; i < vec.size(); i++) {
            if (vec[i]) {
                return i * B + std::countr_zero(vec[i]);
            }
        }
    }

    inline bool empty() const {
        for (const udi& v : vec) {
            if (v) {
                return false;
            }
        }
        return true;
    }

	inline int count_ones() {
        int count = 0;
        for (int i = 0; i < vec.size(); ++i) {
            if(vec[i])
            count += std::popcount(vec[i]);
        }
        return count;
    }
};