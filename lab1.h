#ifndef LAB1_H
#define LAB1_H

int get_random_number(int min, int max);
void two_regret_heuristics(std::vector<std::vector<double>> &distance_matrix, std::vector<int> &indexes_of_first_cycle, std::vector<int> &indexes_of_second_cycle, float w1 = 1, float w2 = 1);
void two_regret_heuristics_with_weights(std::vector<std::vector<double>> &distance_matrix, std::vector<int> &indexes_of_first_cycle, std::vector<int> &indexes_of_second_cycle, double weight1, double weight2);
#endif