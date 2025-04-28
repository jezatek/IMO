#ifndef LAB2_H
#define LAB2_H
#include <vector>
using namespace std;
vector<int>::iterator nextIter(vector<int>::iterator iter, vector<int> *tab);
vector<int>::iterator prevIter(vector<int>::iterator iter, vector<int> *tab);
void changeWierzholek(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, bool steepest);
void changeEdge(vector<vector<double>> &distance_matrix, vector<int> &indexes_of_first_cycle, vector<int> &indexes_of_second_cycle, bool steepest);
#endif