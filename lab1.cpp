#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <cstdlib>
#include <algorithm>
#include <iterator>
using namespace std;

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
/// @brief Generates random result
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be empty at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be empty at start
void randomRes(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle)
{
    int n = distance_matrix.size();
    vector<int> indexes;
    for (int i = 0; i < n; i++)
    {
        indexes.push_back(i);
    }
    random_shuffle(indexes.begin(), indexes.end());
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

/// @brief Generates greedy result
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be empty at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be empty at start
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

/// @brief Generates 2-regret result
/// @param distance_matrix n*n array of distances between points
/// @param indexes_of_first_cycle Result first cycle -> needs to be empty at start
/// @param indexes_of_second_cycle Result second cycle -> needs to be empty at start
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