//
// Created by mrfarinq on 16.06.17.
//

#ifndef SDIZO_3_TRAVELLINGSALESMANPROBLEM_H
#define SDIZO_3_TRAVELLINGSALESMANPROBLEM_H

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <climits>
#include <queue>
#include <algorithm>
#include <map>

using namespace std;
struct Node {
    int cost;                     // path cost till this node
    std::pair<int, int> path;        // path segment (example : 0 1)
    bool bar = false;           // included or excluded path segment
    long parentNodeKey = -1;     // parent id node (in tree)
};

class TravellingSalesmanProblem {

private:
    int amountOfCities;
    int *result;
    int **arrayOfMatrixOfCities;
    int min;
    int reference = std::numeric_limits<int>::max();
    int infinity = 999999999;    // smallest cost found
    int *optimalWay_Solution;
    int length;
    bool setGreedyAlgorithm;
    std::string fileName;
    std::string graphType;


public:
    void DeleteTravellingSalesman();

    void GenerateRandomCities(int amountOfCities = 0, int maxDistanceBetweenCity = 99);

    void PrintCitiesForTheTravellingSalesman();

    void GreedyAlgorithm();

    void PrintSolution();

    std::deque<Node> tree;
    std::vector<int> lastTour;                           // last found tour
    long long int minCost;

    ~TravellingSalesmanProblem();

    TravellingSalesmanProblem();

    void generateFile();

    void readDataFromFile(std::string data);

    void showResult();

    void showCollection();

    void greedyAlgorithm();

    void removeColumn(vector<vector<int>> &matrix, int column, int amountOfCitiesInActualSubset);

    void removeRow(vector<vector<int>> &matrix, int row, int amountOfCitiesInActualSubset);

    void BranchAndBoundAlgorithm();

    int calculateRegret(vector<vector<int>> &m, std::pair<int, int> &path, std::pair<int, int> &pos,
                        int amountOfCitiesInActualSubset);

    int getMinimumRow(vector<vector<int>> &cities, int row, int skippedColumn, int amountOfCitiesInActualSubset);

    int getMinimumColumn(vector<vector<int>> &cities, int row, int skippedColumn, int amountOfCitiesInActualSubset);

    int reduceRow(vector<vector<int>> &cites, int row, int amountOfCitiesInActualSubset);

    int reduceColumn(vector<vector<int>> &cities, int col, int amountOfCitiesInActualSubset);

    void LoadArrayOfMatrixOfCities(long long int **_cities, int _amountOfCities,
                                   std::string _fileName, std::string _graphType);

    int reduceMatrix(vector<vector<int>> &cities, int amountOfCitiesInActualSubset);

    void removeSubTour(vector<vector<int>> &activeRoute, int index, std::pair<int, int> &path,
                       int amountOfCitiesInActualSubset);

    void addLastPath(vector<vector<int>> &m);

    std::vector<int> orderPath(int index, int begin);

    void bruteForce();

    void permutation(int *tab, int i);

    void checkTourCost(vector<vector<int>> &cities);

};


#endif //SDIZO_3_TRAVELLINGSALESMANPROBLEM_H
