#pragma once
#include <vector>
//using std::udi;
using udi = unsigned int;
class LinearHeap {
   public:
    udi max_val;
    udi max_len;

    udi cnt;

    udi min_val;
    std::vector<udi> head;
    std::vector<udi> pre, next;

    LinearHeap(udi max_val, udi max_len)
        : max_val(max_val),
          max_len(max_len),
          cnt(0),
          min_val(max_val),
          head(max_val + 1, max_len),
          pre(max_len + 1, max_len),
          next(max_len + 1, max_len) {}

    inline void link(udi l, udi r) {
        next[l] = r;
        pre[r] = l;
    }

    inline void add(udi val, udi id) {
        link(id, head[val]);
        head[val] = id;
        if (val < min_val) {
            min_val = val;
        }
        cnt++;
    }

    inline void del(udi val, udi id) {
        if (head[val] == id) {
            head[val] = next[id];
            pre[next[id]] = max_len;
        } else {
            link(pre[id], next[id]);
        }
        pre[id] = max_len;
        next[id] = max_len;
        cnt--;
    }

    inline udi pop_min() {
        while (min_val != max_val && head[min_val] == max_len) {
            min_val++;
        }
        if (min_val == max_val) {
            return max_len;
        }
        udi id = head[min_val];
        del(min_val, id);
        return id;
    }

    inline bool empty() {
        return !cnt;
    }
};