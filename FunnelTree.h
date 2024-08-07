#ifndef FUNNEL_TREE_H
#define FUNNEL_TREE_H
#include <string>
#include <functional>
#include <array>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <utility>
#include <memory>
using namespace std;


class Point_3 {
public:
    const double x, y, z;
    const int name;
    Point_3(const double &x, const double &y, const double &z, const int &name) : x(x), y(y), z(z), name(name) {}
    bool operator==(const Point_3 &other) const { return x == other.x && y == other.y && z == other.z; }
    string str() const { return "p(" + to_string(name) + '|' + to_string(x) + ' ' + to_string(y) + ' ' + to_string(z) + ')'; }
    string shortstr() const { return to_string(name); }
};

typedef array<const Point_3*, 2> Edge;
bool operator==(const Edge &a, const Edge &b);
string to_string(const Edge& e);

typedef array<const Point_3*, 3> Triangle;
string to_string(const Triangle &t);

struct Hasher {
    size_t operator()(const Point_3 &p) const {
        return ((hash<double>()(p.x) ^ (hash<double>()(p.y) << 1)) >> 1) ^ (hash<double>()(p.z) << 1);
    }
    size_t operator()(const Edge &e) const { return Hasher()(*e[0]) ^ Hasher()(*e[1]); }
    size_t operator()(const Triangle &t) const { return ((Hasher()(*t[0]) ^ (Hasher()(*t[1]) << 1)) >> 1) ^ (Hasher()(*t[2]) << 1); }
};

typedef unordered_map<Point_3, vector<int>, Hasher> PointDict;
typedef unordered_map<Edge, vector<int>, Hasher> EdgeDict;

class Funnel {
private:
    array<Funnel*, 2> child;
    shared_ptr<vector<int>> triangles;
    const Point_3 *x;
    bool removed;

public:
    const Point_3 &p, &q;
    double length_sp, length_pq, spq_angle, psq_angle, psw_angle, topright_angle;

    Funnel(const Point_3 &p, const Point_3 &q, const Point_3 &x, const shared_ptr<vector<int>> &sequence,
            const double length_sp = 0, const double &length_pq = 0, const double &spq_angle = 0,
            const double &psq_angle = 0, const double &psw_angle = 0, const double &topright_angle = 0)
            : p(p), q(q), x(&x), triangles(sequence), child({nullptr, nullptr}), length_sp(length_sp), length_pq(length_pq),
            spq_angle(spq_angle), psq_angle(psq_angle), psw_angle(psw_angle), topright_angle(topright_angle), removed(0) {}


    void set_x(const Point_3& x) { this->x = &x; }
    bool getRemoved() const { return removed; }
    void addChildren(const bool &i, const Funnel &f);
    const auto& getChild(const int &i) { return child[i]; }
    void delchild(const bool &i);
    const auto& get_x() const { return *x; }
    const auto& getSequence() const { return triangles; }
    void add_to_sequence(const int &i) { triangles->push_back(i); }
    string strshort() const { return "(" + p.shortstr() + ' ' + q.shortstr() + ' ' + x->shortstr() + ' ' + to_string(removed) + ')'; }
    string str() const;
    ~Funnel() { for(auto i : child) if (i != nullptr) delete i; }
};

class TriangleMesh {
private:
    int countFace;
    vector<Triangle> dict_triangle;
    PointDict dict_vertices;
    EdgeDict dict_edges;

public:
    TriangleMesh(ifstream &input);
    vector<Triangle> getFaces() const { return dict_triangle; }
    auto getVertices() const { return dict_vertices; }
    auto getEdges() const { return dict_edges; }
    int size() const { return countFace; }
};

vector<vector<Funnel*>> FunnelTree(const Point_3 &s, const TriangleMesh& mesh);

#endif