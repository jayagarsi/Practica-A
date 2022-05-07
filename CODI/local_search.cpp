/***************************************************************************
    local_search.cpp
    (C) 2021 by C.Blum & M.Blesa
    
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "./HEADERS/Timer.h"
#include "./HEADERS/Random.h"
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>
#include <set>
#include <limits>
#include <iomanip>

// global variables concerning the random number generator (in case needed)
time_t t;
Random* rnd;

// Data structures for the problem data
int n_of_nodes;
int n_of_arcs;
vector< set<int> > neighbors;

// string for keeping the name of the input file
string inputFile;

// number of applications of local search
int n_apps = 1;

// dummy parameters as examples for creating command line parameters -> 
// see function read_parameters(...)
int dummy_integer_parameter = 0;
int dummy_double_parameter = 0.0;


inline int stoi(string &s) {

  return atoi(s.c_str());
}

inline double stof(string &s) {

  return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFile = argv[++iarg];
        
        // reading the number of applications of local search 
        // from the command line (if provided)
        else if (strcmp(argv[iarg],"-n_apps")==0) n_apps = atoi(argv[++iarg]); 
        
        // example for creating a command line parameter param1 -> 
        // integer value is stored in dummy_integer_parameter
        else if (strcmp(argv[iarg],"-isol")==0) {
            dummy_integer_parameter = atoi(argv[++iarg]); 
        }
        // example for creating a command line parameter param2 -> 
        // double value is stored in dummy_double_parameter
        else if (strcmp(argv[iarg],"-param2")==0) {
            dummy_double_parameter = atof(argv[++iarg]);  
        }
        iarg++;
    }
}

void printSolution(const set<int> &solution) {
    for (int i : solution) cout << i << " ";
    cout << endl;
}

int count_neighbors(const vector<set<int>> &solution, int neighbor) {
    int actual_neighbors = 0;
    for (set<int> vertex : solution) if (vertex.find(neighbor) != vertex.end()) actual_neighbors++;
    return actual_neighbors;
}

bool is_solution(const vector<int> &neighborsSize, const vector<set<int>> &neighbors, int deletedVertex) {
    for (int neighbor : neighbors[deletedVertex]) {
        int needed_neighbors = round(neighbors[neighbor].size() / 2.0);
        if (neighborsSize[neighbor] - 1 < needed_neighbors) return false;
    }
    return true;
}

vector<pair<int, int>> generate_successors(const set<int> &solutionVertexs, const vector<int> &neighborsSize, pair<int, int> solutionInfo, const vector<set<int>> &neighbors) {
    vector<pair<int, int>> successors;
    for (int i : solutionVertexs) {
        pair<int, int> successor(solutionInfo.first, i);
        if(is_solution(neighborsSize, neighbors, i)) {
            successor.first--;
            successors.push_back(successor);
        }
    }
    return successors;
}

double getHeuristic(const vector<set<int>> &neighbors, pair<int, int> solutionInfo) {

    int solutionSize = solutionInfo.first;
    int lastDeletedVertex = solutionInfo.second;

    if (neighbors[lastDeletedVertex].size() != 0)
        return solutionSize - 1 / double(neighbors[lastDeletedVertex].size());
    else
        return solutionSize - 1;
}

void MergeSortedIntervals(vector<set<int> >& v, int s, int m, int e) {

    // temp is used to temporary store the vector obtained by merging
    // elements from [s to m] and [m+1 to e] in v
    vector<set<int> > temp;

    int i, j;
    i = s;
    j = m + 1;

    while (i <= m && j <= e) {

        if (v[i].size() <= v[j].size()) {
            temp.push_back(v[i]);
            ++i;
        }
        else {
            temp.push_back(v[j]);
            ++j;
        }

    }

    while (i <= m) {
        temp.push_back(v[i]);
        ++i;
    }

    while (j <= e) {
        temp.push_back(v[j]);
        ++j;
    }

    for (int i = s; i <= e; ++i)
        v[i] = temp[i - s];

}

// the MergeSort function
// Sorts the array in the range [s to e] in v using
// merge sort algorithm
void MergeSort(vector<set<int> >& v, int s, int e) {
    if (s < e) {
        int m = (s + e) / 2;
        MergeSort(v, s, m);
        MergeSort(v, m + 1, e);
        MergeSortedIntervals(v, s, m, e);
    }
}

int computeH(const set<int>& v, set<int>& S) {

    int compt = 0;
    for (int a : v) {                                   // per tot vei de a
        std::set<int>::iterator it = S.find(a);
        if (it != S.end()) ++compt;                     // Veins de v a S
    }

    int n = ceil(float(v.size())/2);                    // upper bound de grau(v)/2
    return n-compt;
}

int cover_degree(const vector<set<int> >& G, set<int>& S, set<int>& SC, int i) {

    int argMax = 0;
    int aux = 0;
    int nCoverMax = -1;

    for (int a : G[i]) {

        std::set<int>::iterator it = SC.find(a);
        int covered = 0;
        if (it != SC.end()) {
            for (int aresta : G[a]) {
                if (computeH(G[aresta], S) > 0) {
                    ++covered;
                    aux = aresta;
                }
            }

            if (covered > nCoverMax) {
                nCoverMax = covered;
                argMax = a;
            }
        }
    }
    return argMax;
}

void greedy(const vector<set<int> >& G, set<int>& S) {

    int n = G.size();
    set<int> SC;
    for (int i = 0; i < n; ++i) {
        SC.insert(i);
    }

    for (int i = 0; i < n; ++i) {

        int p = computeH(G[i], S);                  // mirem si el vertex i esta cobert
        if (p > 0) {
        
            for (int j = 0; j < p; ++j) {
                int argmax = cover_degree(G, S, SC, i);             // cover degree dels vertexs adjacents a i que estan en SC
                S.insert(argmax);
                SC.erase(argmax);
            }
        }
    }
}

/**********
Main function
**********/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);
    
    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // initializing the random number generator.
    // A random number in (0,1) is obtained with: double rnum = rnd->next();
    rnd = new Random((unsigned) time(&t));
    rnd->next();

    // vectors for storing the result and the computation time 
    // obtained by the <n_apps> applications of local search
    vector<double> results(n_apps, std::numeric_limits<int>::max());
    vector<double> times(n_apps, 0.0);

    // opening the corresponding input file and reading the problem data
    ifstream indata;
    indata.open(inputFile.c_str());
    if(not indata) { // file couldn't be opened
        cout << "Error: file could not be opened" << endl;
    }

    indata >> n_of_nodes;
    indata >> n_of_arcs;
    neighbors = vector< set<int> >(n_of_nodes);
    int u, v;
    while(indata >> u >> v) {
        neighbors[u - 1].insert(v - 1);
        neighbors[v - 1].insert(u - 1);
    }
    indata.close();

    // main loop over all applications of local search
    for (int na = 0; na < n_apps; ++na) {

        // the computation time starts now
        Timer timer;

        // Example for requesting the elapsed computation time at any moment: 
        // double ct = timer.elapsed_time(Timer::VIRTUAL);

        cout << "start application " << na + 1 << endl;

        // HERE GOES YOUR LOCAL SEARCH METHOD

        // The starting solution for local search may be randomly generated, 
        // or you may incorporate your greedy heuristic in order to produce 
        // the starting solution.

        set<int> solutionVertexs;

        if (dummy_integer_parameter == 0) {                     // solucio inicial amb tots els vertexs
            for (int i = 0; i < neighbors.size(); i++) solutionVertexs.insert(i);
        }

        else if (dummy_integer_parameter == 1) {
            MergeSort(neighbors, 0, neighbors.size()-1);            // ordenem primer els vertexs respecte el seu grau de manera ascendent
            greedy(neighbors, solutionVertexs);
        }

        else cout << "BAD USAE OF isol" << endl;

        pair<int, int> solutionInfo(solutionVertexs.size(), -1); // first --> size, second --> lastDeletedVertex

        vector<int> neighborsSize;
        for (int i = 0; i < neighbors.size(); i++) neighborsSize.push_back(neighbors[i].size());
        double HeurActual = getHeuristic(neighbors, solutionInfo);

        bool end = false;
        while (not end) {
            vector<pair<int, int>> successorsInfo = generate_successors(solutionVertexs, neighborsSize ,solutionInfo, neighbors);
            bool foundSuccessor = false;
            for (pair<int, int> successorInfo : successorsInfo) {
                double HeurPossible = getHeuristic(neighbors, successorInfo);
                if (HeurPossible < HeurActual) {
                    solutionInfo = successorInfo;
                    HeurActual = HeurPossible;
                    foundSuccessor = true;
                }
            }
            if (foundSuccessor) {
                solutionVertexs.erase(solutionInfo.second);
                for (int neighbor : neighbors[solutionInfo.second]) neighborsSize[neighbor]--;
            }
            else end = true;

            // cout << "New solution: " << endl;
            // printSolution(solution);
            cout << "SOLUTION SIZE: " << solutionInfo.first << endl;
        }

        double ct = timer.elapsed_time(Timer::VIRTUAL);

        // Whenever you move to a new solution, first take the computation 
        // time as explained above. Say you store it in variable ct.
        // Then, write the following to the screen: 
        //cout << "Solution: " << endl;
        //printSolution(solutionVertexs);

        cout << "value: " << solutionVertexs.size();
        cout << "\ttime: " << ct << endl;

        // When a local minimum is reached, store the value of the 
        // corresponding solution in vector results: 
        results[na] = solutionVertexs.size();
        
        // Finally store the needed computation time (that is, the time 
        // measured once the local minimum is reached) in vector times: 
        times[na] = ct;

        cout << "end application " << na + 1 << endl;
    }

    // calculating the average of the results and computation times, and 
    // their standard deviations, and write them to the screen
    double r_mean = 0.0;
    int r_best = std::numeric_limits<int>::max();
    double t_mean = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        r_mean = r_mean + results[i];
        if (int(results[i]) < r_best) r_best = int(results[i]);
        t_mean = t_mean + times[i];
    }
    r_mean = r_mean/double(results.size());
    t_mean = t_mean/double(times.size());
    double rsd = 0.0;
    double tsd = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        rsd = rsd + pow(results[i]-r_mean,2.0);
        tsd = tsd + pow(times[i]-t_mean,2.0);
    }
    rsd = rsd/double(results.size());
    if (rsd > 0.0) {
        rsd = sqrt(rsd);
    }
    tsd = tsd/double(results.size());
    if (tsd > 0.0) {
        tsd = sqrt(tsd);
    }
    cout << r_best << "\t" << r_mean << "\t" << rsd << "\t";
    cout << t_mean << "\t" << tsd << endl;
}

