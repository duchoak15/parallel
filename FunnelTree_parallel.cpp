#define _USE_MATH_DEFINES
#include "FunnelTree_parallel.h"
#include <cmath>
#include <algorithm>

const double M_PI2 = M_PI - 1e-12;

double cal_distance(const Point_3 &a, const Point_3 &b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

bool operator==(const Edge &a, const Edge &b) { return (*a[0] == *b[0] && *a[1] == *b[1]) || (*a[0] == *b[1] && *a[1] == *b[0]); }
string to_string(const Edge& e) { return "(" + e[0]->shortstr() + ' ' + e[1]->shortstr() + ')'; }
string to_string(const Triangle &t) { return "tri(" + t[0]->shortstr() + ' ' + t[1]->shortstr() + ' ' + t[2]->shortstr() + ')'; }

TriangleMesh::TriangleMesh(ifstream &input) {
    int v, temp;
    input >> v >> countFace >> temp;

    vector<const Point_3*> pVec;

    for (int i = 0; i < v; i++) {
        double x, y, z;
        input >> x >> y >> z;
        auto temp2 = dict_vertices.emplace(Point_3(x, y, z, i), vector<int>());
        pVec.push_back(&temp2.first->first);
    }

    for (int i = 0; i < countFace; i++) {
        int x, y, z;
        input >> temp >> x >> y >> z;

        const array<const Point_3*, 3> point = {pVec[x], pVec[y], pVec[z]};
        for (auto j : point) dict_vertices[*j].push_back(i);
        dict_triangle.emplace_back(point);

        for (auto j : array<Edge, 3>{Edge{pVec[x], pVec[y]}, Edge{pVec[y], pVec[z]}, Edge{pVec[z], pVec[x]}}) {
            dict_edges.emplace(j, vector<int>());
            dict_edges[j].push_back(i);
        }

    }
}

void Funnel::delchild(const bool &i) {
    if (child[i] == nullptr) return;
    child[i]->removed = 1;
    child[i]->delchild(0);
    child[i]->delchild(1);
}

void Funnel::addChildren(const bool &i, const Funnel &f) {
    if (child[i] != nullptr) delete child[i];
    child[i] = new Funnel(f);
}

pair<const Point_3*, size_t> determine_v(const Point_3 &x, const Point_3 &q, const vector<int> &S,
                                   const EdgeDict &dict_edges, const vector<Triangle> &triangles) {
    auto temp = dict_edges.find(Edge{&x, &q});

    if (temp != dict_edges.end() && !(x == q))
        for (int f : temp->second) {
            if (f == S.back()) continue;
            if (find(S.begin(), S.end(), f) != S.end()) break;
            for (auto ve : triangles[f]) 
                if (!(*ve == x || *ve == q)) return make_pair(ve, f);
        }

    return make_pair(nullptr, -1);
}

double cal_angle_with_length(const double &vx, const double &vp, const double &xp) {
    return acos((vx * vx + vp * vp - xp * xp) / (2 * vx * vp));
}

double cal_pv(const double &pqv_angle, const double &pq, const double &vq) {
    return sqrt(pq * pq + vq * vq - 2 * pq * vq * cos(pqv_angle));
}

double cal_angle3(const Point_3 &a, const Point_3 &b, const Point_3 &c) {
    const double ab = cal_distance(a, b), bc = cal_distance(b, c), ca = cal_distance(c, a);
    return acos((ab * ab + bc * bc - ca * ca) / (2 * ab * bc));
}

void Clip_off_funnels(const double &sv, const double &pvs, Funnel &funnel, Funnel &funnel_old) {
    const Funnel *const funnel_pv_old = funnel_old.getChild(0), *const funnel_vq_old = funnel_old.getChild(1);
    const double l2 = funnel_vq_old->length_sp, pvz2 = asin(funnel_pv_old->length_sp * sin(funnel_pv_old->spq_angle) / l2);

    if (l2 > sv) {
        if (pvs > pvz2) funnel_old.delchild(1);
        else if (pvz2 > pvs) funnel_old.delchild(0);
    }
    else if (sv > l2) {
        if (pvs > pvz2) funnel.delchild(0);
        else if (pvz2 > pvs) funnel.delchild(1);
    }
    else if (pvs > pvz2) {
        funnel.delchild(0);
        funnel_old.delchild(1);
    }
    else if (pvz2 > pvs) {
        funnel.delchild(1);
        funnel_old.delchild(0);
    }
}

vector<vector<Funnel*>> FunnelTree(const Point_3 &s, const TriangleMesh &mesh) {
    const PointDict dict_vertices = mesh.getVertices();
    if (dict_vertices.find(s) == dict_vertices.end()) return {};

    
    const EdgeDict edges_dict = mesh.getEdges();
    const vector<Triangle> dict_triangle = mesh.getFaces();
    const int size = mesh.size();

    vector<vector<vector<Funnel*>>> tree;
    int sz_parallel = 0;
    for (int face : dict_vertices.at(s)) {
            const Triangle f = dict_triangle[face];
            const Point_3 *p, *q;
            if (s == *f[0]) {
                q = f[1];
                p = f[2];
            } else if (s == *f[1]) {
                q = f[2];
                p = f[0];
            } else {
                q = f[0];
                p = f[1];
            }

            const double spq_angle = cal_angle3(s, *p, *q), psw_angle = cal_angle3(*p, s, *q);
            Funnel *const new_funnel =  new Funnel(*p, *q, *p, make_shared<vector<int>>(vector<int>{face}),
                                                  cal_distance(s, *p), cal_distance(*p, *q), spq_angle, psw_angle, psw_angle);
            tree.push_back({{new_funnel}});
            sz_parallel++;
        }

    vector<int> depth(sz_parallel);
    vector<bool> flag(sz_parallel);
    bool sflag;
    unordered_map<Triangle, Funnel*, Hasher> funnel_has_2_children;

    for (int i = 0; i < size - 1; ++i)
    {   
        sflag = 1;
        #pragma omp parallel for
        for (int id = 0; id < sz_parallel; ++id) {
            if (flag[id]) continue; // using continue because openmp forbids break clause
            vector<Funnel*> l_f_iplus1;
            for (Funnel *funnel : tree[id][i]) {
                if (funnel->getRemoved()) continue;

                const Point_3 &p = funnel->p, &q = funnel->q, *x, *v;
                double spv, length_pv, length_vq;
                vector<int> S_new = *funnel->getSequence();
                short int sign;
                while (1) {
                    x = &funnel->get_x();

                    auto temp = determine_v(*x, q, S_new, edges_dict, dict_triangle);
                    v = temp.first;
                    if (v == nullptr) break;
                    const int &f = temp.second;
                    if (find(S_new.begin(), S_new.end(), f) == S_new.end()) S_new.push_back(f);

                    length_vq = cal_distance(*v, q);
                    double topright_angle = funnel->topright_angle + cal_angle3(*x, q, *v);
                    sign = 1;
                    if (topright_angle >= M_PI2) {
                        topright_angle = M_PI * 2 - topright_angle;
                        sign = -1;
                    }

                    length_pv = cal_pv(topright_angle, funnel->length_pq, length_vq);
                    const double vpq = cal_angle_with_length(length_pv, funnel->length_pq, length_vq);
                    spv = funnel->spq_angle + vpq * sign;

                    if (spv < M_PI2) break;
                    funnel->set_x(*v);
                    funnel->add_to_sequence(f);
                    funnel->topright_angle = topright_angle;
                }

                if (v == nullptr) continue;

                const double length_sv = cal_pv(spv, funnel->length_sp, length_pv),
                             psv = cal_angle_with_length(funnel->length_sp, length_sv, length_pv),
                             pvq = cal_angle_with_length(length_pv, length_vq, funnel->length_pq),
                             psw_1 = min(funnel->psw_angle, psv),
                             top_right_new = max(cal_angle3(*x, *v, q) - pvq * sign, 0.0);
                
                
                funnel->addChildren(0, Funnel(p, *v, *x, make_shared<vector<int>>(S_new),
                                              funnel->length_sp, length_pv, spv, psv, psw_1, top_right_new));
                l_f_iplus1.push_back(funnel->getChild(0));

                if (psv < funnel->psw_angle) {
                    const double pvs = cal_angle_with_length(length_pv, length_sv, funnel->length_sp),
                                 vsq = funnel->psq_angle - psv,
                                 vsw = funnel->psw_angle - psv;

                    funnel->addChildren(1, Funnel(*v, q, *v, funnel->getChild(0)->getSequence(), length_sv, length_vq, pvq - pvs, vsq, vsw));
                    l_f_iplus1.push_back(funnel->getChild(1));

                    #pragma omp critical
                    {   // emplace is not thread-safe
                        auto temp = funnel_has_2_children.emplace(Triangle{&p, v, &q}, funnel);
                        if (!temp.second) Clip_off_funnels(length_sv, pvs, *funnel, *temp.first->second);
                    }
                }
            }
            if (l_f_iplus1.empty()) {flag[id] = 1; continue;}
            tree[id].push_back(l_f_iplus1);
            depth[id]++; sflag = 0;
        }

        if (sflag) break;
    }

    int max_depth = *max_element(depth.begin(), depth.end());
    vector<vector<Funnel*>> combine_tree(max_depth + 1);

    for (int id = 0; id < sz_parallel; ++id)
        for (int i = 0; i <= depth[id]; ++i)
            combine_tree[i].insert(combine_tree[i].end(), tree[id][i].begin(), tree[id][i].end());

    return combine_tree;
}

string Funnel::str() const {
    if (child[1] != nullptr) return strshort() + " - " + child[0]->strshort() + ' ' + child[1]->strshort();
    if (child[0] != nullptr) return strshort() + " - " + child[0]->strshort();
    return strshort();
}