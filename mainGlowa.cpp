#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>

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
    return sqrt((n2.x - n1.x) * (n2.x - n1.x) + (n2.y - n1.y) * (n2.y - n1.y));
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
            int distance = distance_matrix[first_index][node_index] + distance_matrix[second_index][node_index];

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
            int distance = distance_matrix[first_index][node_index] + distance_matrix[second_index][node_index];

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

double highest_two_regret(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_cycle, int inserted_node_index, int &insertion_index)
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
    return second_lowest_cost - lowest_cost;
}

void two_regret_heuristics(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle)
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
                int regret = highest_two_regret(distance_matrix, indexes_of_first_cycle, j, insertion_index);

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
                int regret = highest_two_regret(distance_matrix, indexes_of_second_cycle, j, insertion_index);

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
int main()
{
    srand(0);

    vector<Node> nodes = read_coordinates_from_file("kroA200.tsp");
    vector<vector<double>> distance_matrix = create_distance_matrix(nodes);

    vector<int> bestFirst;
    vector<int> bestSec;

    // nearest_neighbour(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
    // greedy_cycle(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
    // two_regret_heuristics(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);

    double weight1 = -1;
    double weight2 = 1;
    // two_regret_heuristics_with_weights(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle, weight1, weight2);
    int sum = 0;
    int mini = INT32_MAX;
    int maxi = 0;
    for (int i = 0; i < 50; ++i)
    {
        vector<int> indexes_of_first_cycle;
        vector<int> indexes_of_second_cycle;
        two_regret_heuristics(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
        int res = resultFromCycles(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
        sum += res;
        if (res < mini)
        {
            mini = res;
            bestFirst = indexes_of_first_cycle;
            bestSec = indexes_of_second_cycle;
        }
        maxi = max(maxi, res);
    }
    cout << "Mean" << sum / 50.0 << endl;
    cout << "min" << mini << endl;
    cout << "max" << maxi << endl;
    saveResults("trw.txt", bestFirst, bestSec, nodes);
    return 0;
}
