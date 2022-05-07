#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Timer.h"
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
#include <ilcplex/ilocplex.h>

ILOSTLBEGIN

int n_of_nodes; 
int n_of_arcs;  
vector< set<int> > neighbors; 

string inputFile;

double time_limit = 3200.0;


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
		else if (strcmp(argv[iarg],"-t")==0) time_limit = atof(argv[++iarg]);
		iarg++;
	}
}

ILOSOLVECALLBACK4(loggingCallback,
		Timer&, timer,
		double&, time_stamp,
		double&, result,
		double&, gap) {

	if (hasIncumbent()) {
		IloNum nv = getIncumbentObjValue();
		double newTime = timer.elapsed_time(Timer::VIRTUAL);
		double newGap = 100.0*getMIPRelativeGap();
		if (result > double(nv)) {
			cout << "value " << nv << "\ttime " << newTime << "\tgap " << newGap << endl;
			result = double(nv);
			time_stamp = newTime;
			gap = newGap;
		}
	}
}


void run_cplex(Timer& timer) {

	double result = std::numeric_limits<int>::max();
	double time_stamp = 0.0;
	double gap = 0.0;

	IloEnv env;
	env.setOut(env.getNullStream());
	try{

		IloModel model(env);

		IloNumVarArray x(env, n_of_nodes, 0, 1, ILOINT);

		IloExpr obj(env);

		for (int i = 0; i < n_of_nodes; ++i) obj += x[i];

		model.add(IloMinimize(env, obj));
		obj.end();

		for (int i = 0; i < n_of_nodes; ++i) {
			IloExpr expr(env);
			for (set<int>::iterator sit = neighbors[i].begin(); sit != neighbors[i].end(); ++sit) expr += x[*sit];
			int value = ceil(float(neighbors[i].size()) / 2.0);
			model.add(expr >= value);
		}

		IloCplex cpl(model);

		cpl.setParam(IloCplex::TiLim, time_limit);
		cpl.setParam(IloCplex::EpGap, 0.0);
		cpl.setParam(IloCplex::EpAGap, 0.0);
		cpl.setParam(IloCplex::Threads, 1);
		cpl.setWarning(env.getNullStream());
		cpl.use(loggingCallback(env, timer, time_stamp, result, gap));
		cpl.solve();

		if (cpl.getStatus() == IloAlgorithm::Optimal or 
				cpl.getStatus() == IloAlgorithm::Feasible) 
		{
			double newTime = timer.elapsed_time(Timer::VIRTUAL);
			double lastVal = double(cpl.getObjValue());
			double lastGap = 100.0*cpl.getMIPRelativeGap();
			if (lastGap < 0.0) lastGap *= -1.0;
			cout << "value " << lastVal;
			cout << "\ttime " << newTime;
			cout << "\tgap " << lastGap << endl;
			if (cpl.getStatus() == IloAlgorithm::Optimal) {
				cout << "optimality proven" << endl;
			}
			cout << "nodes/vertices in the solution: (";
			bool first = true;
			int counter = 0;
			for (int i = 0; i < n_of_nodes; ++i) {
				IloNum xval = cpl.getValue(x[i]);
				if (xval > 0.9) {
					++counter;
					if (first) {
						cout << i;
						first = false;
					}
					else cout << "," << i;
				}
			}
			cout << ")" << endl;
			cout << "Total amount of nodes: " << counter << endl;
		}
		env.end();
	}
	catch(IloException& e) {
		cerr  << " ERROR: " << e << endl;
	}
	env.end();
}


// Main function

int main( int argc, char **argv ) {

	read_parameters(argc,argv);

	std::cout << std::setprecision(2) << std::fixed;

	ifstream indata;
	indata.open(inputFile.c_str());
	if (not indata) { 
		cout << "Error: file could not be opened" << endl;
	}

	indata >> n_of_nodes;
	indata >> n_of_arcs;
	neighbors = vector< set<int> >(n_of_nodes);
	int u, v;
	while (indata >> u >> v) {
		neighbors[u - 1].insert(v - 1);
		neighbors[v - 1].insert(u - 1);
	}
	indata.close();
	Timer timer;
	run_cplex(timer);

}

