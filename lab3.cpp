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

            if ((cond1 && cond3) || (cond2 && cond4))
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
