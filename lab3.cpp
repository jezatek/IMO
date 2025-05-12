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
using namespace std;
struct Move
{
    double delta;
    vector<int> *tab1;
    vector<int> *tab2;
    int i;
    int j;
    int prev_i;
    int next_i;
    int prev_j;
    int next_j;
    bool operator<(const Move &m) const
    {
        return delta < m.delta;
    }
    Move(vector<int>::iterator &iIter, vector<int>::iterator &jIter, vector<int> &t1, vector<int> &t2, vector<vector<double>> &distance_matrix)
    {
        i = *iIter;
        j = *jIter;
        tab1 = &t1;
        tab2 = &t2;
        prev_i = *prevIter(iIter, &t1);
        prev_j = *prevIter(jIter, &t2);
        next_i = *nextIter(iIter, &t1);
        next_j = *nextIter(jIter, &t2);
        if (tab1 == tab2)
        {
            delta = 0;
            delta -= distance_matrix[*iIter][*nextIter(iIter, tab1)];
            delta -= distance_matrix[*jIter][*nextIter(jIter, tab2)];
            delta += distance_matrix[*iIter][*jIter];
            delta += distance_matrix[*nextIter(iIter, tab1)][*nextIter(jIter, tab2)];
        }
        else
        {
            delta = 0;
            delta = -distance_matrix[*iIter][*nextIter(iIter, tab1)] - distance_matrix[*iIter][*prevIter(iIter, tab1)];
            delta -= (distance_matrix[*jIter][*nextIter(jIter, tab2)] + distance_matrix[*jIter][*prevIter(jIter, tab2)]);
            delta += (distance_matrix[*jIter][*nextIter(iIter, tab1)] + distance_matrix[*jIter][*prevIter(iIter, tab1)]);
            delta += (distance_matrix[*iIter][*nextIter(jIter, tab2)] + distance_matrix[*iIter][*prevIter(jIter, tab2)]);
        }
    }
};
float recalcDelta(vector<int>::iterator &iIter, vector<int>::iterator &jIter, vector<int> &t1, vector<int> &t2, vector<vector<double>> &distance_matrix)
{
    float delta;
    if (t1 == t2)
    {
        delta = 0;
        delta -= distance_matrix[*iIter][*nextIter(iIter, &t1)];
        delta -= distance_matrix[*jIter][*nextIter(jIter, &t2)];
        delta += distance_matrix[*iIter][*jIter];
        delta += distance_matrix[*nextIter(iIter, &t1)][*nextIter(jIter, &t2)];
    }
    else
    {
        delta = 0;
        delta = -distance_matrix[*iIter][*nextIter(iIter, &t1)] - distance_matrix[*iIter][*prevIter(iIter, &t1)];
        delta -= (distance_matrix[*jIter][*nextIter(jIter, &t2)] + distance_matrix[*jIter][*prevIter(jIter, &t2)]);
        delta += (distance_matrix[*jIter][*nextIter(iIter, &t1)] + distance_matrix[*jIter][*prevIter(iIter, &t1)]);
        delta += (distance_matrix[*iIter][*nextIter(jIter, &t2)] + distance_matrix[*iIter][*prevIter(jIter, &t2)]);
    }
    return delta;
}
void updateLM(vector<pair<int, vector<int> *>> &els, vector<int>::iterator iter, vector<int> *tab, multiset<Move> &LM, vector<vector<double>> &distance_matrix)
{
    for (auto &el : els)
    {
        auto otherIter = el.second->begin() + el.first;
        if (iter == otherIter || nextIter(iter, tab) == otherIter || prevIter(iter, tab) == otherIter)
            continue;
        auto i1 = nextIter(iter, tab);
        auto i2 = prevIter(iter, tab);
        Move m(i1, otherIter, *tab, *el.second, distance_matrix);
        if (m.delta < 0)
        {
            LM.insert(m);
        }
        Move m1(iter, otherIter, *tab, *el.second, distance_matrix);
        if (m1.delta < 0)
        {
            LM.insert(m1);
        }
        Move m2(i2, otherIter, *tab, *el.second, distance_matrix);
        if (m2.delta < 0)
        {
            LM.insert(m2);
        }
    }
}
void fillLM(multiset<Move> &LM, multiset<Move> &dLM)
{
    for (auto &m : dLM)
    {
        LM.insert(m);
    }
    dLM.clear();
}

/// @brief Performs changing of the edges with memory to reduce time
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be filled at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be filled at start
void changeEdgeMemory(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle)
{
    vector<pair<int, vector<int> *>> els;
    for (int i = 0; i < indexes_of_first_cycle.size(); i++)
        els.emplace_back(i, &indexes_of_first_cycle);
    for (int i = 0; i < indexes_of_second_cycle.size(); i++)
        els.emplace_back(i, &indexes_of_second_cycle);

    multiset<Move> LM;
    multiset<Move> delayedLM;
    for (int p = 0; p < els.size(); ++p)
    {
        for (int q = p + 1; q < els.size(); ++q)
        {
            auto iIter = els[p].second->begin() + els[p].first;
            auto jIter = els[q].second->begin() + els[q].first;

            Move m(iIter, jIter, *els[p].second, *els[q].second, distance_matrix);
            if (m.delta < 0)
            {
                LM.insert(m);
            }
        }
    }
    // cout << "LM size: " << LM.size() << endl;
    bool improved = true;
    // int ananas = 0;
    while (!LM.empty())
    {
        improved = false;
        // ananas++;
        auto it = LM.begin();

        Move m = *it;
        LM.erase(it);

        auto iIter = find(m.tab1->begin(), m.tab1->end(), m.i);
        auto jIter = find(m.tab2->begin(), m.tab2->end(), m.j);

        if (iIter == m.tab1->end() || jIter == m.tab2->end())
            continue;

        int ni = *nextIter(iIter, m.tab1);
        int nj = *nextIter(jIter, m.tab2);
        int pi = *prevIter(iIter, m.tab1);
        int pj = *prevIter(jIter, m.tab2);

        if (m.tab1 == m.tab2)
        {
            bool ok1 = (ni == m.next_i);
            bool ok2 = (pi == m.next_i);
            bool ok3 = (nj == m.next_j);
            bool ok4 = (pj == m.next_j);

            if ((ok1 && ok3) || (ok2 && ok4))
            {
                if (recalcDelta(iIter, jIter, *m.tab1, *m.tab2, distance_matrix) >= 0)
                    continue;

                if (iIter < jIter)
                    reverse(nextIter(iIter, m.tab1), jIter + 1);
                else
                    reverse(nextIter(jIter, m.tab1), iIter + 1);

                improved = true;
                fillLM(LM, delayedLM);
                updateLM(els, iIter, m.tab1, LM, distance_matrix);
                updateLM(els, jIter, m.tab2, LM, distance_matrix);
            }
            else if ((ok1 && ok4) || (ok2 && ok3))
            {
                Move newM = m;
                delayedLM.insert(newM);
            }
        }
        else
        {
            bool cond1 = (ni == m.next_i && pi == m.prev_i);
            bool cond2 = (pi == m.next_i && ni == m.prev_i);
            bool cond3 = (nj == m.next_j && pj == m.prev_j);
            bool cond4 = (pj == m.next_j && nj == m.prev_j);

            if ((cond1 || cond2) && (cond3 || cond4))
            {
                if (recalcDelta(iIter, jIter, *m.tab1, *m.tab2, distance_matrix) >= 0)
                    continue;

                swap(*iIter, *jIter);

                improved = true;

                fillLM(LM, delayedLM);
                updateLM(els, iIter, m.tab1, LM, distance_matrix);
                updateLM(els, jIter, m.tab2, LM, distance_matrix);
            }
        }
    }
    // cout << "ananans" << ananas;
}

/// @brief Performs changing of the edges with NN
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be filled at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be filled at start
void changeEdgeCandidates(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, bool steepest, int numNN = 10)
{
    int n = distance_matrix.size();

    vector<vector<int>> candidates(n);
    for (int i = 0; i < n; i++)
    {
        vector<pair<double, int>> dists;
        for (int j = 0; j < n; j++)
        {
            if (i != j)
                dists.emplace_back(distance_matrix[i][j], j);
        }
        sort(dists.begin(), dists.end());
        for (int k = 0; k < numNN; k++)
            candidates[i].push_back(dists[k].second);
    }

    bool improved = true;
    while (improved)
    {
        improved = false;

        unordered_map<int, pair<vector<int> *, vector<int>::iterator>> point_map;
        for (auto it = indexes_of_first_cycle.begin(); it != indexes_of_first_cycle.end(); ++it)
            point_map[*it] = {&indexes_of_first_cycle, it};
        for (auto it = indexes_of_second_cycle.begin(); it != indexes_of_second_cycle.end(); ++it)
            point_map[*it] = {&indexes_of_second_cycle, it};

        float bestSoFar = 0;
        vector<int>::iterator besti, bestj;
        vector<int> *bestTab1 = nullptr, *bestTab2 = nullptr;

        for (auto &[p, info1] : point_map)
        {
            auto &[tab1, iIter] = info1;
            for (int q : candidates[p])
            {
                if (point_map.find(q) == point_map.end())
                    continue;

                auto &[tab2, jIter] = point_map[q];

                if (tab1 == tab2)
                {
                    auto pi = prevIter(iIter, tab1);
                    auto pj = prevIter(jIter, tab2);
                    float delta1 = recalcDelta(iIter, jIter, *tab1, *tab2, distance_matrix);
                    float delta2 = recalcDelta(pi, pj, *tab1, *tab2, distance_matrix);

                    if (steepest)
                    {
                        if (delta1 < bestSoFar)
                        {
                            improved = true;
                            besti = iIter;
                            bestj = jIter;
                            bestTab1 = tab1;
                            bestTab2 = tab2;
                            bestSoFar = delta1;
                        }
                        if (delta2 < bestSoFar)
                        {
                            improved = true;
                            besti = prevIter(iIter, tab1);
                            bestj = prevIter(jIter, tab2);
                            bestTab1 = tab1;
                            bestTab2 = tab2;
                            bestSoFar = delta2;
                        }
                    }
                    else
                    {
                        if (delta1 < 0)
                        {
                            improved = true;
                            if (iIter < jIter)
                                reverse(nextIter(iIter, tab1), jIter + 1);
                            else
                                reverse(nextIter(jIter, tab1), iIter + 1);
                            break;
                        }
                        if (delta2 < 0)
                        {
                            improved = true;
                            auto piIter = prevIter(iIter, tab1);
                            auto pjIter = prevIter(jIter, tab2);
                            if (piIter < pjIter)
                                reverse(nextIter(piIter, tab1), pjIter + 1);
                            else
                                reverse(nextIter(pjIter, tab1), piIter + 1);
                            break;
                        }
                    }
                }
                else
                {
                    auto ni = nextIter(iIter, tab1);
                    auto pj = prevIter(jIter, tab2);
                    float delta1 = recalcDelta(pj, iIter, *tab2, *tab1, distance_matrix);
                    float delta2 = recalcDelta(ni, jIter, *tab1, *tab2, distance_matrix);

                    if (steepest)
                    {
                        if (delta1 < bestSoFar)
                        {
                            improved = true;
                            besti = prevIter(jIter, tab2);
                            bestj = iIter;
                            bestTab1 = tab2;
                            bestTab2 = tab1;
                            bestSoFar = delta1;
                        }
                        if (delta2 < bestSoFar)
                        {
                            improved = true;
                            besti = nextIter(iIter, tab1);
                            bestj = jIter;
                            bestTab1 = tab1;
                            bestTab2 = tab2;
                            bestSoFar = delta2;
                        }
                    }
                    else
                    {
                        if (delta1 < 0)
                        {
                            improved = true;
                            swap(*prevIter(jIter, tab2), *iIter);
                            break;
                        }
                        if (delta2 < 0)
                        {
                            improved = true;
                            swap(*nextIter(iIter, tab1), *jIter);
                            break;
                        }
                    }
                }
            }
            if (improved && !steepest)
                break;
        }

        if (improved && steepest)
        {
            if (bestTab1 == bestTab2)
            {
                if (besti < bestj)
                    reverse(nextIter(besti, bestTab1), bestj + 1);
                else
                    reverse(nextIter(bestj, bestTab1), besti + 1);
            }
            else
            {
                swap(*besti, *bestj);
            }
        }
    }
}
