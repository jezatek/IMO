#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <cstdlib>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <set>
#include <map>
#include <omp.h>

#include "lab1.h"
#include "lab2.h"
#include "lab3.h"
#include "lab4.h"

using namespace std;

/// @brief for x time each iteration deletes 30 points from cycles.
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be empty at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be empty at start
/// @param LS Using Local Search each iteration
float paraLNS(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, long long time, int &iterations, bool LS = true)
{
    float best_result = 999999999;
#pragma omp parallel
    {
        vector<int> pom_indexes_of_first_cycle;
        vector<int> pom_indexes_of_second_cycle;
        vector<int> local_first, local_second;
        two_regret_heuristics(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
        changeEdgeMemory(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);

        float localbest = resultFromCycles(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
        local_first = pom_indexes_of_first_cycle;
        local_second = pom_indexes_of_second_cycle;
        int iter = 0;
        auto start = chrono::steady_clock::now();
        while (chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count() < time)
        {
            iter++;
            pom_indexes_of_first_cycle = local_first;
            pom_indexes_of_second_cycle = local_second;
            shreadCycle(pom_indexes_of_first_cycle, 20);
            shreadCycle(pom_indexes_of_second_cycle, 20);

            two_regret_heuristics(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle, false);

            if (LS)
                changeEdgeCandidates(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle, true, 7);
            float result = resultFromCycles(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
            // cout << "best_result: " << best_result <<endl;
            // cout << "result: " << result <<endl;

            if (result < localbest)
            {
                localbest = result;
                local_first = pom_indexes_of_first_cycle;
                local_second = pom_indexes_of_second_cycle;
            }
        }
#pragma omp critical
        {
            if (localbest < best_result)
            {
                best_result = localbest;
                indexes_of_first_cycle = local_first;
                indexes_of_second_cycle = local_second;
                iterations += iter;
            }
        }
    }
    return best_result;
}
int shared_edges(const vector<int> &a1, const vector<int> &a2, const vector<int> &b1, const vector<int> &b2)
{
    auto edges = [](const vector<int> &cycle)
    {
        set<pair<int, int>> s;
        for (int i = 0; i < cycle.size(); i++)
        {
            int u = cycle[i], v = cycle[(i + 1) % cycle.size()];
            if (u > v)
                swap(u, v);
            s.insert({u, v});
        }
        return s;
    };

    auto ea = edges(a1);
    auto eb = edges(a2);
    ea.merge(edges(b1));
    eb.merge(edges(b2));

    int count = 0;
    for (auto &e : ea)
        if (eb.count(e))
            count++;
    return count;
}

int same_cycle_vertex_pairs(const vector<int> &a1, const vector<int> &a2, const vector<int> &b1, const vector<int> &b2)
{
    vector<int> group_a(200), group_b(200);
    for (int v : a1)
        group_a[v] = 0;
    for (int v : a2)
        group_a[v] = 1;
    for (int v : b1)
        group_b[v] = 0;
    for (int v : b2)
        group_b[v] = 1;

    int count = 0;
    for (int i = 0; i < 200; i++)
        for (int j = i + 1; j < 200; j++)
            if (group_a[i] == group_a[j] && group_b[i] == group_b[j])
                count++;
    return count;
}
struct Solution
{
    vector<int> f, s;
    float r;
};

void generateStats(vector<vector<double>> &dist, const string goodCycles, const string outputFile)
{
    int N = 1000;
    vector<Solution> o;
    for (int i = 0; i < N; i++)
    {
        vector<int> f, s;
        two_regret_heuristics(dist, f, s, true);
        changeEdge(dist, f, s, false);
        float cost = resultFromCycles(dist, f, s);
        o.push_back({f, s, cost});
    }
    Solution goodS;
    float pom;
    ifstream file(goodCycles);
    for (int i = 0; i < 100; i++)
    {
        file >> pom;
        goodS.f.push_back(--pom);
    }
    for (int i = 0; i < 100; i++)
    {
        file >> pom;
        goodS.s.push_back(--pom);
    }
    goodS.r = resultFromCycles(dist, goodS.f, goodS.s);
    file.close();
    ofstream fout(outputFile);
    for (int i = 0; i < N; i++)
    {
        auto &sol = o[i];

        int edge_sim_to_best = shared_edges(sol.f, sol.s, goodS.f, goodS.s);
        int cycle_sim_to_best = same_cycle_vertex_pairs(sol.f, sol.s, goodS.f, goodS.s);

        double avg_edge_sim = 0, avg_cycle_sim = 0;
        for (int j = 0; j < N; j++)
        {
            if (i == j)
                continue;
            avg_edge_sim += shared_edges(sol.f, sol.s, o[j].f, o[j].s);
            avg_cycle_sim += same_cycle_vertex_pairs(sol.f, sol.s, o[j].f, o[j].s);
        }
        avg_edge_sim /= (N - 1);
        avg_cycle_sim /= (N - 1);

        fout << sol.r << " "
             << edge_sim_to_best << " "
             << cycle_sim_to_best << " "
             << avg_edge_sim << " "
             << avg_cycle_sim << "\n";
    }
    cout << "JEstemTutej2" << endl;

    fout.close();
}