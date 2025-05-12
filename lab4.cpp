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
#include <unordered_map>

#include "lab1.h"
#include "lab2.h"
#include "lab3.h"
using namespace std;

/// @brief Performs 200 times rand -> edge swap with memory and finds max res
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be empty at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be empty at start
float MSLS(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle)
{
    int iterations = 200;
    float best_result = MAXFLOAT;
    for (int i = 0; i < iterations; i++)
    {
        // cout << "iteration: " << i <<endl;

        vector<int> pom_indexes_of_first_cycle;
        vector<int> pom_indexes_of_second_cycle;
        randomRes(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
        changeEdgeMemory(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
        float result = resultFromCycles(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
        if (result < best_result)
        {
            best_result = result;
            indexes_of_first_cycle = pom_indexes_of_first_cycle;
            indexes_of_second_cycle = pom_indexes_of_second_cycle;
        }
    }
    cout << best_result << endl;
    return best_result;
}

void randomLocalChangeIter(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, int number_of_changes = 5)
{
    for (int i = 0; i < number_of_changes; i++)
    {
        if (rand() % 2)
        {
            if (rand() % 2)
            {
                int a = get_random_number(0, indexes_of_first_cycle.size() - 1);
                int b = get_random_number(0, indexes_of_first_cycle.size() - 1);
                if (a > b)
                    swap(a, b);
                if (a < b)
                    reverse(indexes_of_first_cycle.begin() + a, indexes_of_first_cycle.begin() + b);
            }
            else
            {
                int a = get_random_number(0, indexes_of_second_cycle.size() - 1);
                int b = get_random_number(0, indexes_of_second_cycle.size() - 1);
                if (a > b)
                    swap(a, b);
                if (a < b)
                    reverse(indexes_of_second_cycle.begin() + a, indexes_of_second_cycle.begin() + b);
            }
        }
        else
        {
            vector<int>::iterator iter1;
            vector<int>::iterator iter2;
            if (rand() % 2)
                iter1 = indexes_of_first_cycle.begin() + get_random_number(0, indexes_of_first_cycle.size() - 1);
            else
                iter1 = indexes_of_second_cycle.begin() + get_random_number(0, indexes_of_second_cycle.size() - 1);
            if (rand() % 2)
                iter2 = indexes_of_first_cycle.begin() + get_random_number(0, indexes_of_first_cycle.size() - 1);
            else
                iter2 = indexes_of_second_cycle.begin() + get_random_number(0, indexes_of_second_cycle.size() - 1);
            swap(iter1, iter2);
        }
        // get_random_number
    }
}

/// @brief for x time does small perturbations than fixes them
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be empty at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be empty at start
float ILS(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, long long time)
{
    vector<int> pom_indexes_of_first_cycle;
    vector<int> pom_indexes_of_second_cycle;
    randomRes(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
    changeEdgeMemory(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);

    float best_result = resultFromCycles(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
    indexes_of_first_cycle = pom_indexes_of_first_cycle;
    indexes_of_second_cycle = pom_indexes_of_second_cycle;

    auto start = chrono::steady_clock::now();
    while (chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count() < time)
    {
        pom_indexes_of_first_cycle = indexes_of_first_cycle;
        pom_indexes_of_second_cycle = indexes_of_second_cycle;
        randomLocalChangeIter(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle, 5);
        changeEdgeMemory(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
        float result = resultFromCycles(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
        // cout << "best_result: " << best_result <<endl;
        // cout << "result: " << result <<endl;

        if (result < best_result)
        {
            best_result = result;
            indexes_of_first_cycle = pom_indexes_of_first_cycle;
            indexes_of_second_cycle = pom_indexes_of_second_cycle;
        }
    }
    return best_result;
}

void shreadCycle(vector<int> &cycle, int procent)
{
    int n = cycle.size();
    int r = get_random_number(0, n - 1);
    int remove_count = n * procent / 100;
    int end_index = r + remove_count - 1;
    if (end_index >= n)
    {
        cycle.erase(cycle.begin() + r, cycle.end());
        cycle.erase(cycle.begin(), cycle.begin() + (end_index - n + 1));
    }
    else
    {
        cycle.erase(cycle.begin() + r, cycle.begin() + end_index + 1);
    }
}

/// @brief for x time each iteration deletes 30 points from cycles.
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be empty at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be empty at start
/// @param LS Using Local Search each iteration
float LNS(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, long long time, bool LS = true)
{
    vector<int> pom_indexes_of_first_cycle;
    vector<int> pom_indexes_of_second_cycle;
    randomRes(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
    changeEdgeMemory(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);

    float best_result = resultFromCycles(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
    indexes_of_first_cycle = pom_indexes_of_first_cycle;
    indexes_of_second_cycle = pom_indexes_of_second_cycle;
    int iter = 0;
    auto start = chrono::steady_clock::now();
    while (chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count() < time)
    {
        iter++;
        pom_indexes_of_first_cycle = indexes_of_first_cycle;
        pom_indexes_of_second_cycle = indexes_of_second_cycle;
        shreadCycle(pom_indexes_of_first_cycle, 30);
        shreadCycle(pom_indexes_of_second_cycle, 30);

        two_regret_heuristics(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle, false);

        if (LS)
            changeEdgeMemory(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
        float result = resultFromCycles(distance_matrix, pom_indexes_of_first_cycle, pom_indexes_of_second_cycle);
        // cout << "best_result: " << best_result <<endl;
        // cout << "result: " << result <<endl;

        if (result < best_result)
        {
            best_result = result;
            indexes_of_first_cycle = pom_indexes_of_first_cycle;
            indexes_of_second_cycle = pom_indexes_of_second_cycle;
        }
    }
    return best_result;
}
