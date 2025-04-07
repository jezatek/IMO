#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <cstdlib>
#include <algorithm>
#include <iterator>
#include <chrono>
#define nrOfTrials 100
using namespace std;

struct Node
{
    int id;
    int x;
    int y;
};

vector<Node> read_coordinates_from_file(const string &filename)
{
    ifstream file(filename);
    vector<Node> nodes;

    if (!file.is_open())
    {
        cerr << "Error opening file!" << endl;
        return nodes;
    }

    string line;
    bool reading_coords = false;

    while (getline(file, line))
    {
        if (line == "NODE_COORD_SECTION")
        {
            reading_coords = true;
            continue;
        }

        if (reading_coords)
        {
            if (line == "EOF")
                break;

            istringstream iss(line);
            Node node;
            if (iss >> node.id >> node.x >> node.y)
            {
                nodes.push_back(node);
            }
        }
    }

    file.close();
    return nodes;
}
double euclidean_distance(Node n1, Node n2)
{
    return floor(sqrt((n2.x - n1.x) * (n2.x - n1.x) + (n2.y - n1.y) * (n2.y - n1.y)) + 0.5);
}

vector<vector<double>> create_distance_matrix(vector<Node> nodes)
{
    vector<vector<double>> matrix;
    matrix.resize(nodes.size());
    for (int i = 0; i < matrix.size(); i++)
        matrix[i].resize(matrix.size());

    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix[i].size(); j++)
        {
            matrix[i][j] = euclidean_distance(nodes[i], nodes[j]);
        }
    }

    return matrix;
}

// int get_random_number(int min, int max) {
//     random_device rd;
//     default_random_engine generator(rd());
//     uniform_int_distribution<int> distribution(min, max);
//
//     return distribution(generator);
// }

int get_random_number(int min, int max)
{
    return min + rand() % (max - min + 1);
}

int get_index_of_furthest_node(vector<vector<double>> &distance_matrix, vector<bool> &used_nodes, int index)
{
    int furthest_node_index = -1;
    int furthest_distance = 0;

    for (int i = 0; i < distance_matrix[0].size(); i++)
    {
        if ((distance_matrix[index][i] > furthest_distance || furthest_node_index == -1) && used_nodes[i] == false)
        {
            furthest_node_index = i;
            furthest_distance = distance_matrix[index][i];
        }
    }
    return furthest_node_index;
}

int get_index_of_closest_node(vector<vector<double>> &distance_matrix, vector<bool> &used_nodes, int index)
{
    int closest_node_index = -1;
    int closest_distance = 99999999;

    for (int i = 0; i < distance_matrix[0].size(); i++)
    {
        if ((distance_matrix[index][i] < closest_distance || closest_node_index == -1) && used_nodes[i] == false)
        {
            closest_node_index = i;
            closest_distance = distance_matrix[index][i];
        }
    }

    return closest_node_index;
}

void nearest_neighbour(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle)
{
    indexes_of_first_cycle.resize(distance_matrix.size() - distance_matrix.size() / 2, -1);
    indexes_of_second_cycle.resize(distance_matrix.size() / 2, -1);

    vector<bool> used_nodes(distance_matrix.size(), false);

    int first_node_index = get_random_number(0, distance_matrix.size() - 1);
    indexes_of_first_cycle[0] = first_node_index;
    used_nodes[first_node_index] = true;

    int furthest_node_index = get_index_of_furthest_node(distance_matrix, used_nodes, first_node_index);
    indexes_of_second_cycle[0] = furthest_node_index;
    used_nodes[furthest_node_index] = true;

    int i = 1;

    while (true)
    {
        if (i > indexes_of_first_cycle.size() - 1)
            break;
        int closest_node_index = get_index_of_closest_node(distance_matrix, used_nodes, indexes_of_first_cycle[i - 1]);
        indexes_of_first_cycle[i] = closest_node_index;
        used_nodes[closest_node_index] = true;

        if (i > indexes_of_second_cycle.size() - 1)
            break;
        closest_node_index = get_index_of_closest_node(distance_matrix, used_nodes, indexes_of_second_cycle[i - 1]);
        indexes_of_second_cycle[i] = closest_node_index;
        used_nodes[closest_node_index] = true;

        i++;
    }
}

int get_index_closest_to_pair(vector<vector<double>> &distance_matrix, vector<bool> &used_nodes, int index1, int index2)
{
    int closest_node_index = -1;

    int closest_distance = 99999999;
    int distance;
    for (int i = 0; i < distance_matrix[0].size(); i++)
    {
        distance = distance_matrix[index1][i] + distance_matrix[index2][i];
        if ((distance < closest_distance || closest_node_index == -1) && used_nodes[i] == false)
        {
            closest_node_index = i;
            closest_distance = distance;
        }
    }
    return closest_node_index;
}

void randomRes(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle)
{
    int n = distance_matrix.size();
    vector<int> indexes;
    for (int i = 0; i < n; i++)
    {
        indexes.push_back(i);
    }
    shuffle(indexes.begin(), indexes.end(), mt19937(random_device()()));
    for (int i = 0; i < (n + 1) / 2; i++)
        indexes_of_first_cycle.push_back(indexes[i]);
    for (int i = (n + 1) / 2; i < n; i++)
        indexes_of_second_cycle.push_back(indexes[i]);
}

void stupidInit(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle)
{
    int n = distance_matrix.size();
    vector<int> indexes;
    for (int i = 0; i < n; i++)
    {
        indexes.push_back(i);
    }
    // shuffle(indexes.begin(), indexes.end(), mt19937(random_device()()));
    for (int i = 0; i < (n + 1) / 2; i++)
        indexes_of_first_cycle.push_back(indexes[i]);
    for (int i = (n + 1) / 2; i < n; i++)
        indexes_of_second_cycle.push_back(indexes[i]);
}

void greedy_cycle(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle)
{
    vector<bool> used_nodes(distance_matrix.size(), false);
    int first_cycle_full_size = distance_matrix.size() - distance_matrix.size() / 2;
    int second_cycle_full_size = distance_matrix.size() / 2;

    int first_node_index = get_random_number(0, distance_matrix.size() - 1);
    indexes_of_first_cycle.push_back(first_node_index);
    used_nodes[first_node_index] = true;

    int furthest_node_index = get_index_of_furthest_node(distance_matrix, used_nodes, first_node_index);
    indexes_of_second_cycle.push_back(furthest_node_index);
    used_nodes[furthest_node_index] = true;

    int first_node_pair_index = get_index_of_closest_node(distance_matrix, used_nodes, first_node_index);
    indexes_of_first_cycle.push_back(first_node_pair_index);
    used_nodes[first_node_pair_index] = true;

    int furthest_node_pair_index = get_index_of_closest_node(distance_matrix, used_nodes, furthest_node_index);
    indexes_of_second_cycle.push_back(furthest_node_pair_index);
    used_nodes[furthest_node_pair_index] = true;

    int i = 2;

    int best_insert_index;
    int best_node_index;
    int min_distance;

    while (true)
    {
        if (i > first_cycle_full_size - 1)
            break;

        best_insert_index = -1;
        best_node_index = -1;
        min_distance = 99999999;

        for (int j = 1; j < indexes_of_first_cycle.size(); j++)
        {
            int first_index = indexes_of_first_cycle[j];
            int second_index = indexes_of_first_cycle[j - 1];
            int node_index = get_index_closest_to_pair(distance_matrix, used_nodes, first_index, second_index);
            int distance = -distance_matrix[first_index][second_index] + distance_matrix[first_index][node_index] + distance_matrix[second_index][node_index];

            if (distance < min_distance || best_node_index == -1)
            {
                best_node_index = node_index;
                min_distance = distance;
                best_insert_index = j;
            }
        }

        indexes_of_first_cycle.insert(indexes_of_first_cycle.begin() + best_insert_index, best_node_index);
        used_nodes[best_node_index] = true;

        if (i > second_cycle_full_size - 1)
            break;

        best_insert_index = -1;
        best_node_index = -1;
        min_distance = 99999999;

        for (int j = 1; j < indexes_of_second_cycle.size(); j++)
        {
            int first_index = indexes_of_second_cycle[j];
            int second_index = indexes_of_second_cycle[j - 1];
            int node_index = get_index_closest_to_pair(distance_matrix, used_nodes, first_index, second_index);
            int distance = -distance_matrix[first_index][second_index] + distance_matrix[first_index][node_index] + distance_matrix[second_index][node_index];

            if (distance < min_distance || best_node_index == -1)
            {
                best_node_index = node_index;
                min_distance = distance;
                best_insert_index = j;
            }
        }

        indexes_of_second_cycle.insert(indexes_of_second_cycle.begin() + best_insert_index, best_node_index);
        used_nodes[best_node_index] = true;

        i++;
    }
}

double cost_of_insertion(vector<vector<double>> &distance_matrix, int indexA, int indexB, int inserted_indexC)
{
    double distanceAB = distance_matrix[indexA][indexB];
    double distanceAC = distance_matrix[indexA][inserted_indexC];
    double distanceBC = distance_matrix[indexB][inserted_indexC];

    return distanceAC + distanceBC - distanceAB;
}

double lowestCost(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_cycle, int inserted_node_index, int &insertion_index)
{
    vector<double> insertion_costs(indexes_of_cycle.size(), -1);

    for (int i = 0; i < indexes_of_cycle.size(); i++)
    {
        int insertion_cost = cost_of_insertion(distance_matrix, indexes_of_cycle[i], indexes_of_cycle[(i + 1) % indexes_of_cycle.size()], inserted_node_index);
        insertion_costs[i] = insertion_cost;
    }

    int lowest_cost_index = -1;
    double lowest_cost = -1;

    for (int i = 0; i < insertion_costs.size(); i++)
    {
        if (insertion_costs[i] < lowest_cost || lowest_cost_index == -1)
        {
            lowest_cost_index = i;
            lowest_cost = insertion_costs[i];
        }
    }

    insertion_index = lowest_cost_index + 1;
    return lowest_cost;
}

void two_regret_heuristics(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, float w1 = 1, float w2 = 1)
{
    vector<bool> used_nodes(distance_matrix.size(), false);
    int first_cycle_full_size = distance_matrix.size() - distance_matrix.size() / 2;
    int second_cycle_full_size = distance_matrix.size() / 2;
    int first_node_index = get_random_number(0, distance_matrix.size() - 1);
    indexes_of_first_cycle.push_back(first_node_index);
    used_nodes[first_node_index] = true;

    int furthest_node_index = get_index_of_furthest_node(distance_matrix, used_nodes, first_node_index);
    indexes_of_second_cycle.push_back(furthest_node_index);
    used_nodes[furthest_node_index] = true;

    int first_node_pair_index = get_index_of_closest_node(distance_matrix, used_nodes, first_node_index);
    indexes_of_first_cycle.push_back(first_node_pair_index);
    used_nodes[first_node_pair_index] = true;

    int furthest_node_pair_index = get_index_of_closest_node(distance_matrix, used_nodes, furthest_node_index);
    indexes_of_second_cycle.push_back(furthest_node_pair_index);
    used_nodes[furthest_node_pair_index] = true;

    int i = 4;
    // są już dodane 2 elementy do first cycle i 2 do second cycle
    int best_insert_index;
    int best_node_index;
    int max_regret;

    while (true)
    {
        if (i >= distance_matrix.size())
            break;

        best_insert_index = -1;
        best_node_index = -1;
        max_regret = -1;
        bool firest;

        for (int j = 0; j < used_nodes.size(); j++)
        {
            if (used_nodes[j] == false)
            {
                int insertion_index1;
                int insertion_index2;
                int regret = INT32_MAX;
                int regret2 = INT32_MAX;
                if (indexes_of_first_cycle.size() < first_cycle_full_size)
                {
                    regret = w1 * lowestCost(distance_matrix, indexes_of_first_cycle, j, insertion_index1);
                }
                if (indexes_of_second_cycle.size() < second_cycle_full_size)
                {
                    regret2 = w2 * lowestCost(distance_matrix, indexes_of_second_cycle, j, insertion_index2);
                }
                int func = abs(regret - regret2);
                if (func > max_regret)
                {
                    max_regret = func;
                    best_node_index = j;
                    if (regret < regret2)
                    {
                        best_insert_index = insertion_index1;
                        firest = true;
                    }
                    else
                    {
                        best_insert_index = insertion_index2;
                        firest = false;
                    }
                }
            }
        }
        if (firest)
        {
            indexes_of_first_cycle.insert(indexes_of_first_cycle.begin() + best_insert_index, best_node_index);
        }
        else
            indexes_of_second_cycle.insert(indexes_of_second_cycle.begin() + best_insert_index, best_node_index);
        used_nodes[best_node_index] = true;
        i++;
    }
}

double highest_two_regret_with_weights(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_cycle, int inserted_node_index, int &insertion_index, double weight1, double weight2)
{
    vector<double> insertion_costs(indexes_of_cycle.size() - 1, -1);

    for (int i = 1; i < indexes_of_cycle.size(); i++)
    {
        int insertion_cost = cost_of_insertion(distance_matrix, i, i - 1, inserted_node_index);
        insertion_costs[i - 1] = insertion_cost;
    }

    int lowest_cost_index = -1;
    double lowest_cost = -1;
    int second_lowest_cost_index = -1;
    double second_lowest_cost = -1;

    for (int i = 0; i < insertion_costs.size(); i++)
    {
        if (insertion_costs[i] < lowest_cost || lowest_cost_index == -1)
        {
            lowest_cost_index = i;
            lowest_cost = insertion_costs[i];
        }
    }

    for (int i = 0; i < insertion_costs.size(); i++)
    {
        if ((insertion_costs[i] < second_lowest_cost || second_lowest_cost_index == -1) && i != lowest_cost_index)
        {
            second_lowest_cost_index = i;
            second_lowest_cost = insertion_costs[i];
        }
    }

    insertion_index = lowest_cost_index + 1;
    return weight2 * second_lowest_cost - weight1 * lowest_cost;
}

void two_regret_heuristics_with_weights(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, double weight1, double weight2)
{
    vector<bool> used_nodes(distance_matrix.size(), false);
    int first_cycle_full_size = distance_matrix.size() - distance_matrix.size() / 2;
    int second_cycle_full_size = distance_matrix.size() / 2;

    int first_node_index = get_random_number(0, distance_matrix.size() - 1);
    indexes_of_first_cycle.push_back(first_node_index);
    used_nodes[first_node_index] = true;

    int furthest_node_index = get_index_of_furthest_node(distance_matrix, used_nodes, first_node_index);
    indexes_of_second_cycle.push_back(furthest_node_index);
    used_nodes[furthest_node_index] = true;

    int first_node_pair_index = get_index_of_closest_node(distance_matrix, used_nodes, first_node_index);
    indexes_of_first_cycle.push_back(first_node_pair_index);
    used_nodes[first_node_pair_index] = true;

    int furthest_node_pair_index = get_index_of_closest_node(distance_matrix, used_nodes, furthest_node_index);
    indexes_of_second_cycle.push_back(furthest_node_pair_index);
    used_nodes[furthest_node_pair_index] = true;

    int first_node_third_index = get_index_closest_to_pair(distance_matrix, used_nodes, first_node_index, first_node_pair_index);
    indexes_of_first_cycle.insert(indexes_of_first_cycle.begin() + 1, first_node_third_index);
    used_nodes[first_node_third_index] = true;

    int furthest_node_third_index = get_index_closest_to_pair(distance_matrix, used_nodes, furthest_node_index, furthest_node_pair_index);
    indexes_of_second_cycle.insert(indexes_of_second_cycle.begin() + 1, furthest_node_third_index);
    used_nodes[furthest_node_third_index] = true;

    int i = 3;

    int best_insert_index;
    int best_node_index;
    int max_regret;

    while (true)
    {
        if (i > first_cycle_full_size - 1)
            break;

        best_insert_index = -1;
        best_node_index = -1;
        max_regret = -1;

        for (int j = 0; j < used_nodes.size(); j++)
        {
            if (used_nodes[j] == false)
            {
                int insertion_index;
                int regret = highest_two_regret_with_weights(distance_matrix, indexes_of_first_cycle, j, insertion_index, weight1, weight2);

                if (regret > max_regret || best_node_index == -1)
                {
                    max_regret = regret;
                    best_node_index = j;
                    best_insert_index = insertion_index;
                }
            }
        }

        indexes_of_first_cycle.insert(indexes_of_first_cycle.begin() + best_insert_index, best_node_index);
        used_nodes[best_node_index] = true;

        if (i > second_cycle_full_size - 1)
            break;

        best_insert_index = -1;
        best_node_index = -1;
        max_regret = -1;

        for (int j = 0; j < used_nodes.size(); j++)
        {
            if (used_nodes[j] == false)
            {
                int insertion_index;
                int regret = highest_two_regret_with_weights(distance_matrix, indexes_of_second_cycle, j, insertion_index, weight1, weight2);

                if (regret > max_regret || best_node_index == -1)
                {
                    max_regret = regret;
                    best_node_index = j;
                    best_insert_index = insertion_index;
                }
            }
        }

        indexes_of_second_cycle.insert(indexes_of_second_cycle.begin() + best_insert_index, best_node_index);
        used_nodes[best_node_index] = true;

        i++;
    }
}

float resultFromCycles(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle)
{
    int cs = indexes_of_first_cycle.size();
    int ss = indexes_of_second_cycle.size();
    // cout << cs << " " << ss << endl;
    int total = 0;
    for (int i = 0; i < cs; i++)
    {
        total += distance_matrix[indexes_of_first_cycle[i]][indexes_of_first_cycle[(i + 1) % cs]];
    }
    for (int i = 0; i < ss; i++)
    {
        total += distance_matrix[indexes_of_second_cycle[i]][indexes_of_second_cycle[(i + 1) % ss]];
    }
    return total;
}

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

void saveResults(const string &filename, const vector<int> &first_cycle, const vector<int> &second_cycle, vector<Node> nodes)
{
    ofstream out(filename);
    for (int i = 0; i < first_cycle.size(); i++)
        out << nodes[first_cycle[i]].id << " ";
    out << "\n";
    for (int i = 0; i < second_cycle.size(); i++)
        out << nodes[second_cycle[i]].id << " ";
    out << "\n";
}
void createInitialResult()
{

    vector<Node> nodes = read_coordinates_from_file("kroA200.tsp");
    vector<vector<double>> distance_matrix = create_distance_matrix(nodes);
    vector<Node> nodes2 = read_coordinates_from_file("kroB200.tsp");
    vector<vector<double>> distance_matrix2 = create_distance_matrix(nodes2);

    vector<int> bestFirst;
    vector<int> bestSec;

    // randomRes(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
    // nearest_neighbour(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
    // greedy_cycle(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
    // two_regret_heuristics(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
    // two_regret_heuristics_with_weights(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle, weight1, weight2);

    int sum = 0;
    int mini = INT32_MAX;
    int maxi = 0;
    vector<long long> durations;
    for (int i = 0; i < nrOfTrials; ++i)
    {
        auto start = chrono::high_resolution_clock::now();
        vector<int> indexes_of_first_cycle;
        vector<int> indexes_of_second_cycle;
        // Zmieniasz ponizsze na randomRes / two_regret_heuristics
        two_regret_heuristics(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
        // changeWierzholek(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle, false);
        changeEdge(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle, false);
        // randomChange(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle, 385400);
        int res = resultFromCycles(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
        sum += res;
        if (res < mini)
        {
            mini = res;
            bestFirst = indexes_of_first_cycle;
            bestSec = indexes_of_second_cycle;
        }
        maxi = max(maxi, res);
        durations.push_back(chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() / 1000);
        if (i % 10 == 0)
            cout
                << "Progres:" << i << endl;
    }
    cout << "Mean" << sum / (nrOfTrials + 0.0) << endl;
    cout << "min" << mini << endl;
    cout << "max" << maxi << endl;
    saveResults("trw.txt", bestFirst, bestSec, nodes);

    cout << "MinTime" << *min_element(durations.begin(), durations.end()) << endl;
    cout << "MaxTime" << *max_element(durations.begin(), durations.end()) << endl;
    cout << "MeanTime" << accumulate(durations.begin(), durations.end(), 0LL) / (nrOfTrials + 0.0) << endl;

    durations.clear();
    vector<int> bestFirst2;
    vector<int> bestSec2;
    sum = 0;
    mini = INT32_MAX;
    maxi = 0;
    for (int i = 0; i < nrOfTrials; ++i)
    {
        auto start = chrono::high_resolution_clock::now();
        vector<int> indexes_of_first_cycle;
        vector<int> indexes_of_second_cycle;
        // Zmieniasz ponizsze na randomRes / two_regret_heuristics
        two_regret_heuristics(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle);
        // changeWierzholek(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle, false);
        changeEdge(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle, false);
        // randomChange(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle, 385400);
        int res = resultFromCycles(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle);
        sum += res;
        if (res < mini)
        {
            mini = res;
            bestFirst2 = indexes_of_first_cycle;
            bestSec2 = indexes_of_second_cycle;
        }
        maxi = max(maxi, res);
        durations.push_back(chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() / 1000);
        if (i % 10 == 0)
            cout << "Progres:" << i << endl;
    }
    cout << "Mean" << sum / (nrOfTrials + 0.0) << endl;
    cout << "min" << mini << endl;
    cout << "max" << maxi << endl;
    saveResults("trw2.txt", bestFirst2, bestSec2, nodes2);

    cout << "MinTime" << *min_element(durations.begin(), durations.end()) << endl;
    cout << "MaxTime" << *max_element(durations.begin(), durations.end()) << endl;
    cout << "MeanTime" << accumulate(durations.begin(), durations.end(), 0LL) / (nrOfTrials + 0.0) << endl;
}
int main()
{
    // srand(0);
    createInitialResult();
    return 0;
}
