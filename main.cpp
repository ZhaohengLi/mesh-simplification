
#include "Vector.h"
#include "cal.h"

#define INFD 1e8
#define BUFFER_SIZE 1024
#define TOLERATE 2.0

typedef std::pair<int, int> Edge;

class Model {
    std::vector< Vector > vertex;
    std::vector<bool> removed;
    std::vector< std::set<Edge> > face;
    std::set<Edge> edge;
    std::priority_queue< std::pair<double, Edge> > heap;
    size_t faceN;

    double edgeLen(const Edge &e) { return norm(vertex[e.first] - vertex[e.second]); }

public:
    void clear() {
        vertex.clear();
        removed.clear();
        face.clear();
        edge.clear();
        faceN = 0;
    }

    size_t getEdgeN() { return edge.size(); }
    size_t getVertexN() { return vertex.size(); }
    size_t getFaceN() { return faceN; }

    void loadFromFile(std::string filename) {
        clear();

        char buffer[BUFFER_SIZE];
        FILE *file = fopen(filename.c_str(), "r");
        std::vector<std::string> vertexIn;
        std::vector<std::string> faceIn;
        while (fgets(buffer, BUFFER_SIZE, file) != NULL) {
            int ptr = 0;
            while (buffer[ptr] != 0 && buffer[ptr] != 'v' && buffer[ptr] != 'f' && buffer[ptr] != '#') ptr++;
            if (buffer[ptr] == 'v') vertexIn.push_back(std::string(buffer));
            if (buffer[ptr] == 'f') faceIn.push_back(std::string(buffer));
        }
        fclose(file);

        size_t vertexN = vertexIn.size();
        vertex.resize(vertexN, Vector(3, 0));
        removed.resize(vertexN, false);
        face.resize(vertexN);
        faceN = faceIn.size();

        for (int i = 0; i < vertexN; i++) {
            sscanf(vertexIn[i].c_str(), "%*s%lf%lf%lf", &vertex[i][0], &vertex[i][1], &vertex[i][2]);
        }

        for (const auto &f : faceIn) {
            int v[3];
            sscanf(f.c_str(), "%*s%d%d%d", v, v + 1, v + 2);
            v[0] --; v[1] --; v[2] --;
            face[v[0]].insert(make_pair(v[1], v[2]));
            face[v[1]].insert(make_pair(v[2], v[0]));
            face[v[2]].insert(make_pair(v[0], v[1]));
            std::sort(v, v + 3);
            assert(0 <= v[0] && v[0] < v[1] && v[1] < v[2] && v[2] < vertexN);
            edge.insert(make_pair(v[0], v[1]));
            edge.insert(make_pair(v[1], v[2]));
            edge.insert(make_pair(v[0], v[2]));
        }
    }

    void saveToFile(std::string filename) {
        FILE *file = fopen(filename.c_str(), "w");
        size_t vertexN = vertex.size();
        std::vector<int> vertexID(vertexN, 0);
        int vertexReal = 0;

        for (int i = 0; i < vertexN; i++) {
            if (removed[i]) continue;
            vertexID[i] = ++vertexReal;
            fprintf(file, "v %.8lf %.8lf %.8lf\n", vertex[i][0], vertex[i][1], vertex[i][2]);
        }

        for (int i = 0; i < vertexN; i++) {
            if (removed[i]) continue;
            for (const auto &f : face[i]) {
                assert(!removed[f.first] && !removed[f.second]);
                assert(vertexID[f.first] && vertexID[f.second] && vertexID[i]);
                if (i < f.first && i < f.second) {
                    fprintf(file, "f %d %d %d\n", vertexID[i], vertexID[f.first], vertexID[f.second]);
                }
            }
        }
    }

    std::pair<Vector, double> getPosition(const Edge &e) {
        Matrix q(4, Vector(4, 0));
        for (const auto &f : face[e.first]) {
            auto n = crossProduct(vertex[f.first] - vertex[e.first], vertex[f.second] - vertex[e.first]);
            n = n / norm(n);
            n.push_back(-innerProduct(vertex[e.first], n));
            outerProductFast(n, n, q);
        }
        for (const auto &f : face[e.second]) {
            auto n = crossProduct(vertex[f.first] - vertex[e.second], vertex[f.second] - vertex[e.second]);
            n = n / norm(n);
            n.push_back(-innerProduct(vertex[e.second], n));
            outerProductFast(n, n, q);
        }

        Vector v;
        try {
            v = solveEquation(q, 3);
        } catch(...) {
            v = (vertex[e.first] + vertex[e.second]) / 2;
        }
        if (norm(v - vertex[e.first]) + norm(v - vertex[e.second]) > TOLERATE * norm(vertex[e.first] - vertex[e.second])) {
            v = (vertex[e.first] + vertex[e.second]) / 2;
        }
        v.push_back(1);
        double cost = innerProduct(innerProduct(v, q), v);
        assert(cost > -EPS);
        v.pop_back();
        return make_pair(v, cost);
    }

    std::pair<Edge, Vector> selectEdge(double threshold) {
        Edge idx = make_pair(-1, -1);
        Vector pos;
        std::pair<double, Edge> tmp;
        while (!heap.empty()) {
            tmp = heap.top();
            heap.pop();
            if (edge.find(tmp.second) == edge.end()) continue;
            if (removed[tmp.second.first] || removed[tmp.second.second]) continue;
            if (edgeLen(tmp.second) > threshold) continue;
            auto act = getPosition(tmp.second);
            if (fabs(act.second + tmp.first) > EPS) continue;
            idx = tmp.second;
            pos = act.first;
            break;
        }
        printf("%lf %d %d", -tmp.first, idx.first, idx.second);
        return std::make_pair(idx, pos);
    }

    bool faceReverse(const Edge &e, const Vector &v1, const Vector &v2) {
        const auto &x = vertex[e.first];
        const auto &y = vertex[e.second];
        return innerProduct(crossProduct(x - v1, y - v1), crossProduct(x - v2, y - v2)) < 0;
    }

    void addToHeap(const Edge &e, double threshold) {
        if (edgeLen(e) > threshold) return;
        auto pos = getPosition(e);
        heap.push(make_pair(-pos.second, e));
    }

    void updateNeighborEdge(int v, double threshold) {
        std::set<int> neighbor;
        for (const auto &f : face[v]) {
            neighbor.insert(f.first);
            neighbor.insert(f.second);
        }
        for (auto x : neighbor) {
            addToHeap(make_pair(min(x, v), max(x, v)), threshold);
        }
    }

    void removeEdge(const Edge &e, const Vector &v, double threshold) {
        std::vector<Edge> toRev;
        for (const auto &f : face[e.first]) {
            if (f.first == e.second || f.second == e.second) continue;
            auto reverse = faceReverse(f, vertex[e.first], v);
            if (!reverse) continue;
            toRev.push_back(f);
            assert(face[f.second].find(make_pair(e.first, f.first)) != face[f.second].end());
            face[f.second].erase(make_pair(e.first, f.first));
            face[f.second].insert(make_pair(f.first, e.first));

            assert(face[f.first].find(make_pair(f.second, e.first)) != face[f.first].end());
            face[f.first].erase(make_pair(f.second, e.first));
            face[f.first].insert(make_pair(e.first, f.second));
        }
        for (const auto &f : toRev) {
            face[e.first].erase(f);
            face[e.first].insert(make_pair(f.second, f.first));
        }


        for (const auto &f : face[e.second]) {
            assert(face[f.second].find(make_pair(e.second, f.first)) != face[f.second].end());
            face[f.second].erase(make_pair(e.second, f.first));
            auto reverse = faceReverse(f, vertex[e.second], v);
            if (f.first != e.first && f.second != e.first) {
                if (reverse) {
                    face[f.second].insert(make_pair(f.first, e.first));
                } else {
                    face[f.second].insert(make_pair(e.first, f.first));
                }
            }

            assert(face[f.first].find(make_pair(f.second, e.second)) != face[f.first].end());
            face[f.first].erase(make_pair(f.second, e.second));
            if (f.first != e.first && f.second != e.first) {
                if (reverse) {
                    face[f.first].insert(make_pair(e.first, f.second));
                } else {
                    face[f.first].insert(make_pair(f.second, e.first));
                }
            }

            if (f.first == e.first || f.second == e.first)
                faceN--;
            else {
                if (reverse) {
                    face[e.first].insert(make_pair(f.second, f.first));
                } else {
                    face[e.first].insert(f);
                }
            }

            auto tmp = make_pair(min(e.second, f.first), max(e.second, f.first));
            if (edge.find(tmp) != edge.end())
                edge.erase(tmp);
            tmp = make_pair(min(e.second, f.second), max(e.second, f.second));
            if (edge.find(tmp) != edge.end())
                edge.erase(tmp);
            if (f.first != e.first && f.second != e.first) {
                edge.insert(make_pair(min(e.first, f.first), max(e.first, f.first)));
                edge.insert(make_pair(min(e.first, f.second), max(e.first, f.second)));
            }
        }

        edge.erase(e);
        vertex[e.first] = v;
        vertex[e.second].clear();
        removed[e.second] = true;
        face[e.second].clear();

        std::set<int> neighbor;
        for (const auto &f: face[e.first]) {
            neighbor.insert(f.first);
            neighbor.insert(f.second);
        }
        for (auto nb : neighbor) {
            updateNeighborEdge(nb, threshold);
        }
    }

    void buildHeap(double threshold) {
        while (!heap.empty()) heap.pop();
        for (const auto &e : edge) {
            addToHeap(e, threshold);
        }
    }

    void simplify(size_t target, double threshold) {
        buildHeap(threshold);
        while (faceN > target) {
            printf("%c%zu ", 13, faceN);
            auto e = selectEdge(threshold);
            if (e.first != make_pair(-1, -1))
                removeEdge(e.first, e.second, threshold);
            else {
                printf("%cERROR: No enough edges under threshold.\n", 13);
                printf("Warning: Current result will be save.\n");
                return;
            }
            fflush(stdout);
        }
    }

};

int main(int argc, char **argv) {
    Model model;
    if (argc != 4) {
        printf("Usage:\n ./main [Input Object] [Output Object] [Simplify Rate] [Threshold Value]");
        return 0;
    }
    std::string inputModelFileName(argv[1]);
    std::string outputModelFileName(argv[2]);
    double simplifyRate = atof(argv[3]);

    printf("inputModelFileName: %s\n", inputModelFileName.c_str());
    printf("outputModelFileName: %s\n", outputModelFileName.c_str());
    printf("simplifyRate: %.4lf\n", simplifyRate);
    printf("threshold: %.4lf\n", INFD);
    printf("------------------------------------\n");

    time_t start = time(0);

    model.loadFromFile(inputModelFileName);

    size_t all = model.getFaceN();
    size_t simple = all * simplifyRate;

    printf("vertex: %zu\n", model.getVertexN());
    printf("edge: %zu\n", model.getEdgeN());
    printf("simple / all = %zu / %zu\n", simple, all);
    model.simplify(simple, INFD);

    model.saveToFile(outputModelFileName);
    time_t end = time(0);
    printf("%cSave to [%s] successfully. Time %ld sec.\n", 13, outputModelFileName.c_str(), end - start);
    return 0;
}