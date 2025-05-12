#ifndef LAB3_H
#define LAB3_H
#include <vector>
void changeEdgeMemory(std::vector<std::vector<double>> &distance_matrix, std::vector<int> &indexes_of_first_cycle, std::vector<int> &indexes_of_second_cycle);
void changeEdgeCandidates(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, bool steepest, int numNN = 10);
#endif