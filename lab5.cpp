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

#include "lab1.h"
// #include "lab2.h"
#include "lab3.h"
// #include "lab4.h"
using namespace std;
struct Solution
{
    double res;
    vector<int> tab1;
    vector<int> tab2;
    bool operator<(const Solution &m) const
    {
        return res < m.res;
    }
    Solution(vector<int> &t1, vector<int> &t2, vector<vector<double>> &distance_matrix)
    {
        tab1 = t1;
        tab2 = t2;
        res = resultFromCycles(distance_matrix, t1, t2);
    }
    /// @brief Rekombination from 2 results
    /// @param t1_1 firstCycle from 1 result
    /// @param t1_2 secondCycle from 1 result
    /// @param t2_1 firstCycle from 2 result
    /// @param t2_2 secondCycle from 2 result
    Solution(const Solution &s1, const Solution &s2, vector<vector<double>> &distance_matrix)
    {
        vector<int> newt1;
        vector<int> newt2;
        int n = s1.tab1.size();
        // for (int i = 0; i < t1_1.size(); ++i)
        // {
        //     if (find(t2_1.begin(), t2_1.end(), t1_1[i]) != t2_1.end())
        //         newt1.push_back(t1_1[i]);
        // }
        for (int i = 0; i < n; ++i)
        {
            int a1 = s1.tab1[i];
            int b1 = s1.tab1[(i + 1) % n];
            int a2 = s1.tab2[i];
            int b2 = s1.tab2[(i + 1) % n];

            for (int j = 0; j < n; ++j)
            {
                int c1 = s2.tab1[j];
                int d1 = s2.tab1[(j + 1) % n];
                int c2 = s2.tab2[j];
                int d2 = s2.tab2[(j + 1) % n];

                if ((a1 == c1 && b1 == d1) || (a1 == c2 && b1 == d2))
                {
                    if (newt1.empty() || newt1.back() != a1)
                        newt1.push_back(a1);
                    if (newt1[0] != b1)
                        newt1.push_back(b1);
                }

                if ((a2 == c1 && b2 == d1) || (a2 == c2 && b2 == d2))
                {
                    if (newt2.empty() || newt2.back() != a2)
                        newt2.push_back(a2);
                    if (newt2[0] != b2)
                        newt2.push_back(b2);
                }
            }
        }
        two_regret_heuristics(distance_matrix, newt1, newt2, false);
        tab1 = newt1;
        tab2 = newt2;
        res = resultFromCycles(distance_matrix, tab1, tab2);
    }
};
/// @brief Hybrid evolution alhorithm
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be empty at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be empty at start
float HAE(std::vector<std::vector<double>> &distance_matrix, std::vector<int> &indexes_of_first_cycle, std::vector<int> &indexes_of_second_cycle, long long time, int &iterations, bool LS = true)
{
    int MAX_SIZE = 20;
    iterations = 0;
    auto start = chrono::steady_clock::now();
    set<Solution> solutions;
    for (int i = 0; i < MAX_SIZE + 10; i++)
    {
        vector<int> idx1;
        vector<int> idx2;
        two_regret_heuristics(distance_matrix, idx1, idx2);
        Solution s(idx1, idx2, distance_matrix);
        solutions.insert(s);
    }

    while (solutions.size() > MAX_SIZE)
    {
        solutions.erase(--solutions.end());
    }
    while (chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count() < time)
    {
        iterations++;
        int r1 = get_random_number(0, solutions.size() - 1);
        int r2;
        do
        {
            r2 = get_random_number(0, solutions.size() - 1);
        } while (r1 == r2);
        auto it1 = solutions.begin();
        std::advance(it1, r1);
        auto it2 = solutions.begin();
        std::advance(it2, r2);
        Solution s(*it1, *it2, distance_matrix);
        if (LS)
        {
            changeEdgeMemory(distance_matrix, s.tab1, s.tab2);
        }
        solutions.insert(s);
        while (solutions.size() > MAX_SIZE)
        {
            solutions.erase(--solutions.end());
        }
    }
    indexes_of_first_cycle = solutions.begin()->tab1;
    indexes_of_second_cycle = solutions.begin()->tab2;
    return solutions.begin()->res;
}
