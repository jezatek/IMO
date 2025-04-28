#ifndef LAB1_H
#define LAB1_H

#include <vector>
using namespace std;

int get_random_number(int min, int max);
void randomRes(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle);
void two_regret_heuristics(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, float w1 = 1, float w2 = 1);
void two_regret_heuristics_with_weights(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, double weight1, double weight2);
#endif