//
// Created by 李曌珩 on 2019-06-16.
//

#ifndef HOMEWORK3_VECTOR_H
#define HOMEWORK3_VECTOR_H

#include <cstdio>
#include <string>
#include <vector>
#include <set>
#include <cassert>
#include <cmath>
#include <queue>
#include <algorithm>

class Vector {
private:
    double a[4];
    int _size;

public:
    Vector() { _size = 0; }
    Vector(int n) { _size = n; }
    Vector(int n, double x) {
        _size = n;
        for (int i = 0; i < n; i++)
            a[i] = x;
    }

    void swap(Vector &b) {
        assert(b.size() == _size);
        for (int i = 0; i < _size; i++)
            std::swap(a[i], b.a[i]);
    }

    int size() const { return _size; }
    void clear() { _size = 0; }

    double& operator [] (int x) { return a[x]; }
    const double& operator [] (int x) const { return a[x]; }

    void push_back(double x) { a[_size++] = x; }
    void pop_back() { _size--; }
};

#endif //HOMEWORK3_VECTOR_H
