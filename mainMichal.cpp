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

#include "lab1.h"
#include "lab2.h"
#include "lab3.h"
#include "lab4.h"

#define nrOfTrials 2
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
        // randomRes(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
        // two_regret_heuristics(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
        // changeWierzholek(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle, false);
        // changeEdge(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle, true);
        // changeEdgeMemory(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
        // changeEdgeCandidates(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle, true, 10);
        // randomChange(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle, 385400);
        // MSLS(distance_matrix, indexes_of_first_cycle, indexes_of_second_cycle);
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
        if (i % nrOfTrials / 10 == 0)
            cout
                << "Progres:" << i << endl;
    }
    cout << "Mean " << sum / (nrOfTrials + 0.0) << endl;
    cout << "min " << mini << endl;
    cout << "max " << maxi << endl;
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
        // randomRes(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle);
        // two_regret_heuristics(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle);
        // changeWierzholek(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle, false);
        // changeEdge(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle, true);
        // changeEdgeMemory(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle);
        // changeEdgeCandidates(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle, true, 10);
        // randomChange(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle, 385400);
        // MSLS(distance_matrix2, indexes_of_first_cycle, indexes_of_second_cycle);
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
        if (i % nrOfTrials / 10 == 0)
            cout << "Progres:" << i << endl;
    }
    cout << "Mean " << sum / (nrOfTrials + 0.0) << endl;
    cout << "min " << mini << endl;
    cout << "max " << maxi << endl;
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
