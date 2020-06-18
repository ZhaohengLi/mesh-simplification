//
// Created by 李曌珩 on 2019-06-16.
//

#ifndef HOMEWORK3_CAL_H
#define HOMEWORK3_CAL_H

#include <cstdio>
#include <string>
#include <vector>
#include <set>
#include <cassert>
#include <cmath>
#include <queue>
#include <algorithm>
#include "Vector.h"

#define EPSLOOSE 0.1
#define EPS 1e-8

using std::min;
using std::max;
using std::make_pair;

typedef std::vector<Vector> Matrix;



void printVector(const Vector &v) {
    for (int i = 0; i < v.size(); i++)
        printf("%.4lf\t", v[i]);
    printf("\n");
}

void printMatrx(const Matrix &m) {
    for (auto v : m)
        printVector(v);
}

double norm(const Vector &v) {
    double t = 0;
    for (int i = 0; i < v.size(); i++) t += v[i] * v[i];
    return sqrt(t);
}

Vector crossProduct(const Vector &a, const Vector &b) {
    assert(a.size() == 3 && b.size() == 3);
    Vector c(3);
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}

double innerProduct(const Vector &a, const Vector &b) {
    assert(a.size() == b.size());
    double c = 0;
    for (int i = 0; i < a.size(); i++)
        c += a[i] * b[i];
    return c;
}

Matrix outerProduct(const Vector &a, const Vector &b) {
    Matrix c(a.size(), Vector(b.size(), 0));
    for (int i = 0; i < a.size(); i++)
        for (int j = 0; j < b.size(); j++)
            c[i][j] = a[i] * b[j];
    return c;
}

void outerProductFast(const Vector &a, const Vector &b, Matrix &c) {
    assert(a.size() == c.size());
    if (a.size() == 0) return;
    assert(b.size() == c[0].size());
    for (int i = 0; i < a.size(); i++)
        for (int j = 0; j < b.size(); j++)
            c[i][j] += a[i] * b[j];
}

Vector innerProduct(const Vector &a, const Matrix &b) {
    assert(a.size() == b.size());
    if (a.size() == 0) return Vector();
    Vector c(b[0].size(), 0);
    for (int i = 0; i < b.size(); i++)
        for (int j = 0; j < b[0].size(); j++)
            c[j] += a[i] * b[i][j];
    return c;
}

Vector operator + (const Vector &a, const Vector &b) {
    assert(a.size() == b.size());
    Vector c(a.size());
    for (int i = 0; i < a.size(); i++)
        c[i] = a[i] + b[i];
    return c;
}

Matrix operator + (const Matrix &a, const Matrix &b) {
    assert(a.size() == b.size());
    Matrix c(a.size());
    for (int i = 0; i < a.size(); i++)
        c[i] = a[i] + b[i];
    return c;
}

Vector operator - (const Vector &a, const Vector &b) {
    assert(a.size() == b.size());
    Vector c(a.size());
    for (int i = 0; i < a.size(); i++)
        c[i] = a[i] - b[i];
    return c;
}

Vector operator * (const double &a, const Vector &b) {
    Vector c(b.size());
    for (int i = 0; i < b.size(); i++)
        c[i] = a * b[i];
    return c;
}

Vector operator / (const Vector &a, const double &b) {
    assert(b != 0);
    Vector c(a.size());
    for (int i = 0; i < a.size(); i++)
        c[i] = a[i] / b;
    return c;
}

Vector solveEquation(Matrix m, int n) {
    assert(m.size() >= n);
    if (m.size() == 0) return Vector();
    assert(m[0].size() > n);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++)
            if (fabs(m[i][i]) < fabs(m[j][i])) m[i].swap(m[j]);
        if (fabs(m[i][i]) < EPSLOOSE) throw 200;
        m[i] = m[i] / m[i][i];
        for (int j = i + 1; j < n; j++)
            m[j] = m[j] - m[j][i] * m[i];
    }
    Vector v(n);
    for (int i = n - 1; i >= 0; i--) {
        assert(fabs(m[i][i] - 1) < EPS);
        v[i] = -m[i][n];
        for (int j = i + 1; j < n; j++) {
            v[i] -= m[i][j] * v[j];
        }
    }

    return v;
}


#endif //HOMEWORK3_CAL_H
