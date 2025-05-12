#ifndef LAB4_H
#define LAB4_H
#include <vector>

float MSLS(std::vector<std::vector<double>> &distance_matrix, std::vector<int> &indexes_of_first_cycle, std::vector<int> &indexes_of_second_cycle);
float ILS(std::vector<std::vector<double>> &distance_matrix, std::vector<int> &indexes_of_first_cycle, std::vector<int> &indexes_of_second_cycle, long long time);
#endif
