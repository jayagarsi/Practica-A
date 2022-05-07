#include <iostream>
#include <vector>
#include <set>
#include <stack>
using namespace std;

using Adjacencia = set<int>;
using Graf = vector<Adjacencia>;

bool isMinimalPositiveInfluenceDominatingSet(const Graf& G, const set<int>& D, vector<bool>& visitat) {

    bool esMPIDS = false;

    std::set<int>::iterator it = D.begin();

    while (it != D.end() and not esMPIDS) {                                         // per cada vertex u de D mirem que passa si el treiem
        set<int> Adjacents = G[*it];                                            // adjacencies del vertex u

        for (int aresta : Adjacents) {                                          // per tota adjacencia v de u
            int nAdjacenciesD = 0;
            bool marcat = false;
            for (int a : G[aresta]) {                                           // mirem les adjacencies d'v tret el vertex u

                std::set<int>::iterator i = D.find(a);
                if (i != D.end() and *i != *it) ++nAdjacenciesD;
                //else if (i != D.end() and *i == *it and G[aresta].size() == 1) marcat = true;
            }

            int nAdjacenciesTotals = G[aresta].size() - 1;

            float proporcioAdjaD = proporcioAdjaD = float(nAdjacenciesD)/float(nAdjacenciesTotals);
            if (proporcioAdjaD < 0.5) esMPIDS = true;

        }

        ++it;

    }

    return esMPIDS;
}

bool isPositiveInfluenceDominatingSet(const Graf& G, const set<int>& D, vector<bool>& visitat) {

    stack<int> S;

    S.push(0);
    bool esPIDS = true;                             // Comencem assumint que ho es

    while (not S.empty() and esPIDS) {              // mentre la pila estigui plena i tots els vertexs
                                                    // visitats fins ara tinguin com a minim la meitat de les adjacencies a D
        int v = S.top();
        S.pop();
        int nAdjacenciesD = 0;

        if (not visitat[v]) {

            visitat[v] = true;

            for (int aresta : G[v]) {

                std::set<int>::const_iterator it = D.find(aresta);
                if (it != D.end()) ++nAdjacenciesD;
                S.push(aresta);
                
            }
            int nAdjacenciesTotals = G[v].size();
            float proporcioAdjaD = float(nAdjacenciesD)/float(nAdjacenciesTotals);
            if (proporcioAdjaD < 0.5) esPIDS = false;
        }
    }
    return esPIDS;
}

int main() {

    cout << "Lectura del Graf G = (V,E)" << endl;
    cout << "Insereix el nombre de vertexs n = |V| i el nombre d'arestes m = |E|" << endl;

    int n, m;
    cin >> n >> m;

    Graf G(n);
    cout << "Introdueix m arestes de la forma v u on v,u son vertexs de V" << endl;

    int v, u;
    for (int i = 0; i < m; ++i) {
        cin >> v >> u;
        G[v-1].insert(u-1);
        G[u-1].insert(v-1);
    }

    cout << "Introdueix el nombre de vertexs del conjunt D i seguidament els vertexs que hi pertanyen" << endl;

    int nD;
    cin >> nD;
    set<int> D;

    int w;
    
    for (int i = 0; i < nD; ++i) {
        cin >> w;
        D.insert(w);
    }

    vector<bool> visitat(n, false);

    cout << endl;

    cout << "Comprovant si D es PIDS..." << endl;

    cout << endl;

    if (isPositiveInfluenceDominatingSet(G, D, visitat)) {
        cout << "El conjunt es d'influencia positiva" << endl;
        cout << endl;
        cout << "Comprovant si D es minimal..." << endl;
        cout << endl;
        if (isMinimalPositiveInfluenceDominatingSet(G, D, visitat))
            cout << "El conjunt es d'influencia positiva i es minimal" << endl;
        else cout << "El conjunt es d'influencia positiva pero no es minimal" << endl;
            
    }
    else {
        cout << "El conjunt no es d'influencia positiva" << endl;
    }
}