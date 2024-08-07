/*
compile:            g++ main main.cpp FunnelTree.cpp
run with file.txt:  ./main file.txt
*/
#include "FunnelTree.h"
#include <iostream>
#include <chrono>

size_t size(const vector<vector<Funnel*>> &tree) {
    size_t n = 0;
    for (auto i : tree) n += i.size();
    return n;
}

void testmesh(const char* filename) {
    ifstream file(filename);
    TriangleMesh mesh(file);

    // point - indexes of faces including point
    // point format: p(index|x y z)
    for (auto i : mesh.getVertices()) {
        cout << i.first.str();
        for (auto j : i.second) cout << ' ' << j;
        cout << '\n';
    }
    cout << endl;

    // edge - indexes of faces including edge
    // edge format: vertices' indexes
    for (auto i : mesh.getEdges()) {
        cout << to_string(i.first);
        for (auto j : i.second) cout << ' ' << j;
        cout << '\n';
    }
    cout << endl;

    // triangle format: tri([vertices' indexes])
    for (auto i : mesh.getFaces()) cout << to_string(i) << '\n';
}

void run(const char* filename) {
    cout << "File: " << filename << '\n';
    ifstream file(filename);
    TriangleMesh mesh(file);
    cout << "There are " << mesh.getVertices().size() << " vertices, "
                         << mesh.getFaces().size() << " faces and "
                         << mesh.getEdges().size() << " edges.\n";

    int temp;
    double x, y, z;
    file.clear();
    file.seekg(0);
    file >> temp >> temp >> temp >> x >> y >> z;
    Point_3 p(x, y, z, 0);

    
    auto start = chrono::high_resolution_clock::now();
    vector<vector<Funnel*>> tree = FunnelTree(p, mesh);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start);
    cout << "Funnel tree initialized with " << size(tree) << " nodes in " << duration.count() << " seconds.\n";
}

int main(int argc, char *argv[]) {
    freopen("out.txt", "w", stdout);
    run(argv[1]);
}
