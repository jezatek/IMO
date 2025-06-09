#ifndef LAB6_H
#define LAB6_H
#include <vector>
float paraLNS(std::vector<std::vector<double>> &distance_matrix, std::vector<int> &indexes_of_first_cycle, std::vector<int> &indexes_of_second_cycle, long long time, int &iterations, bool LS = true);
void generateStats(vector<vector<double>> &dist, const string goodCycles, const string outputFile);
#endif
