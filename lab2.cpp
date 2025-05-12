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

#include "lab1.h"

using namespace std;
vector<int>::iterator nextIter(vector<int>::iterator iter, vector<int> *tab)
{
    if (iter + 1 == tab->end())
    {
        return tab->begin();
    }
    else
        return (iter + 1);
}
vector<int>::iterator prevIter(vector<int>::iterator iter, vector<int> *tab)
{
    if (iter == tab->begin())
    {
        return tab->end() - 1;
    }
    else
        return (iter - 1);
}

/// @brief Performs changing of the vetexes
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be filled at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be filled at start
void changeWierzholek(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, bool steepest)
{
    int n = distance_matrix.size();
    bool improved = true;
    // int ananas = 99;
    while (improved)
    {
        // ananas++;
        improved = false;
        // if (ananas % 100 == 0)
        // cout << resultFromCycles(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle) << " " << ananas << " " << improved << endl;
        vector<int>::iterator besti;
        vector<int>::iterator bestj;
        float bestSoFar = 0;
        vector<pair<int, vector<int> *>> all_elements;
        for (int i = 0; i < indexes_of_first_cycle.size(); i++)
            all_elements.emplace_back(i, &indexes_of_first_cycle);
        for (int i = 0; i < indexes_of_second_cycle.size(); i++)
            all_elements.emplace_back(i, &indexes_of_second_cycle);
        vector<pair<pair<int, vector<int> *>, pair<int, vector<int> *>>> pairs;
        for (int i = 0; i < all_elements.size(); i++)
            for (int j = i + 1; j < all_elements.size(); j++)
                pairs.emplace_back(all_elements[i], all_elements[j]);
        shuffle(pairs.begin(), pairs.end(), mt19937(random_device()()));
        for (auto &pr : pairs)
        {
            vector<int> *tab1 = pr.first.second;
            vector<int> *tab2 = pr.second.second;
            auto iIter = tab1->begin() + pr.first.first;
            auto jIter = tab2->begin() + pr.second.first;
            float delta = 0;
            delta = -distance_matrix[*iIter][*nextIter(iIter, tab1)] - distance_matrix[*iIter][*prevIter(iIter, tab1)];
            delta -= (distance_matrix[*jIter][*nextIter(jIter, tab2)] + distance_matrix[*jIter][*prevIter(jIter, tab2)]);
            delta += (distance_matrix[*jIter][*nextIter(iIter, tab1)] + distance_matrix[*jIter][*prevIter(iIter, tab1)]);
            delta += (distance_matrix[*iIter][*nextIter(jIter, tab2)] + distance_matrix[*iIter][*prevIter(jIter, tab2)]);
            if (*nextIter(iIter, tab1) == *jIter)
                delta += 2 * distance_matrix[*iIter][*jIter];
            if (*prevIter(iIter, tab1) == *jIter)
                delta += 2 * distance_matrix[*iIter][*jIter];
            if (steepest)
            {
                if (delta < bestSoFar)
                {
                    improved = true;
                    besti = iIter;
                    bestj = jIter;
                    bestSoFar = delta;
                }
            }
            else
            {
                if (delta < 0)
                {
                    improved = true;
                    // int a = resultFromCycles(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
                    swap(*iIter, *jIter);
                    // if (resultFromCycles(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle) == a)
                    // {
                    //     cout << *iIter << " " << *jIter;
                    //     cin >> ananas;
                    // }
                    // if (ananas % 100 == 0)
                    // {
                    //     cout << a << " " << resultFromCycles(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle) << endl;
                    //     cout << *iIter << " co" << *jIter;
                    // }

                    break;
                    // continue;
                }
            }
        }
        if (steepest && improved)
        {
            swap(*besti, *bestj);
        }
    }
}

/// @brief Performs changing of the edges
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be filled at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be filled at start
void changeEdge(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, bool steepest)
{
    int n = distance_matrix.size();
    bool improved = true;
    // int ananas = 99;
    while (improved)
    {
        // ananas++;
        improved = false;
        // if (ananas % 100 == 0)
        // cout << resultFromCycles(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle) << " " << ananas << " " << improved << endl;
        vector<int>::iterator besti;
        vector<int>::iterator bestj;
        float bestSoFar = 0;
        vector<pair<int, vector<int> *>> all_elements;
        for (int i = 0; i < indexes_of_first_cycle.size(); i++)
            all_elements.emplace_back(i, &indexes_of_first_cycle);
        for (int i = 0; i < indexes_of_second_cycle.size(); i++)
            all_elements.emplace_back(i, &indexes_of_second_cycle);
        vector<pair<pair<int, vector<int> *>, pair<int, vector<int> *>>> pairs;
        for (int i = 0; i < all_elements.size(); i++)
            for (int j = i + 1; j < all_elements.size(); j++)
                pairs.emplace_back(all_elements[i], all_elements[j]);
        shuffle(pairs.begin(), pairs.end(), mt19937(random_device()()));
        for (auto &pr : pairs)
        {
            vector<int> *tab1 = pr.first.second;
            vector<int> *tab2 = pr.second.second;
            auto iIter = tab1->begin() + pr.first.first;
            auto jIter = tab2->begin() + pr.second.first;

            // if (tab1 != tab2)
            //     continue;

            float delta = 0;
            if (tab1 == tab2)
            {
                delta -= distance_matrix[*iIter][*nextIter(iIter, tab1)];
                delta -= distance_matrix[*jIter][*nextIter(jIter, tab2)];
                delta += distance_matrix[*iIter][*jIter];
                delta += distance_matrix[*nextIter(iIter, tab1)][*nextIter(jIter, tab2)];
            }
            else
            {
                delta = -distance_matrix[*iIter][*nextIter(iIter, tab1)] - distance_matrix[*iIter][*prevIter(iIter, tab1)];
                delta -= (distance_matrix[*jIter][*nextIter(jIter, tab2)] + distance_matrix[*jIter][*prevIter(jIter, tab2)]);
                delta += (distance_matrix[*jIter][*nextIter(iIter, tab1)] + distance_matrix[*jIter][*prevIter(iIter, tab1)]);
                delta += (distance_matrix[*iIter][*nextIter(jIter, tab2)] + distance_matrix[*iIter][*prevIter(jIter, tab2)]);
            }
            // cout << "i: " << i << endl;
            // cout << "j: " << j << endl;
            // cout << "delta: " << delta << endl;
            if (steepest)
            {
                if (delta < bestSoFar)
                {
                    improved = true;
                    besti = iIter;
                    bestj = jIter;
                    bestSoFar = delta;
                }
            }
            else
            {
                if (delta < 0)
                {
                    improved = true;

                    // cout << "tutaj1" << endl;
                    // printVector(*tab1);
                    if (tab1 == tab2)
                    {
                        if (iIter < jIter)
                            reverse(nextIter(iIter, tab1), jIter + 1);
                        else
                            reverse(nextIter(jIter, tab1), iIter + 1);
                    }
                    else
                    {
                        swap(*iIter, *jIter);
                    }
                    // printVector(*tab1);

                    break;
                    ;
                }
            }
        }
        if (steepest && improved)
        {
            auto found = find(indexes_of_first_cycle.begin(), indexes_of_first_cycle.end(), *besti);
            auto found2 = find(indexes_of_first_cycle.begin(), indexes_of_first_cycle.end(), *bestj);
            if (found != indexes_of_first_cycle.end() && found2 != indexes_of_first_cycle.end())
            {
                // cout << "tutaj2" << endl;
                // printVector(indexes_of_first_cycle);
                if (besti < bestj)
                    reverse(nextIter(besti, &indexes_of_first_cycle), bestj + 1);
                else
                    reverse(nextIter(bestj, &indexes_of_first_cycle), besti + 1);
                // printVector(indexes_of_first_cycle);
            }
            else if (found == indexes_of_first_cycle.end() && found2 == indexes_of_first_cycle.end())
            {
                // cout << "tutaj3" << endl;
                // printVector(indexes_of_second_cycle);
                if (besti < bestj)
                    reverse(nextIter(besti, &indexes_of_second_cycle), bestj + 1);
                else
                    reverse(nextIter(bestj, &indexes_of_second_cycle), besti + 1);
                // printVector(indexes_of_second_cycle);
            }
            else
            {
                swap(*besti, *bestj);
            }
        }
    }
}

void randomChange(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, long long time)
{
    auto start = chrono::steady_clock::now();
    while (chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count() < time)
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