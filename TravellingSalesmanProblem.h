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

struct Node {
    int cost;                     // path cost till this node
    std::pair<int, int> path;        // path segment (example : 0 1)
    bool bar = false;           // included or excluded path segment
    long parentNodeKey = -1;     // parent id node (in tree)
};

class TravellingSalesmanProblem {

private:
    int amountOfCities;
    int **arrayOfMatrixOfCities;
    int reference = std::numeric_limits<int>::max();
    int *optimalWay_Solution;
    int length;
    bool setGreedyAlgorithm;
    std::string fileName;
    std::string graphType;
    bool randomGeneratorData;
    std::deque<Node> tree;
    std::vector<int> lastTour;

    std::vector<int> orderPath(int index, int begin);

public:
    TravellingSalesmanProblem();

    ~TravellingSalesmanProblem();

    void DeleteTravellingSalesman();

    void LoadArrayOfMatrixOfCities(long long int **_cities, int _amountOfCities,
                                   std::string _fileName, std::string _graphType);

    void GenerateRandomCities(int amountOfCities = 0, int maxDistanceBetweenCity = 99);

    void PrintCitiesForTheTravellingSalesman(bool printInputMatrix);

    void GreedyAlgorithm();

    void PrintSolution();

    void removeColumn(std::vector<std::vector<int>> &matrix, int column);

    void removeRow(std::vector<std::vector<int>> &matrix, int row);

    void BranchAndBoundAlgorithm();

    int calculateRegret(std::vector<std::vector<int>> &m, std::pair<int, int> &path, std::pair<int, int> &pos,
                        int amountOfCitiesInActualSubset);

    int
    getMinimumRow(std::vector<std::vector<int>> &cities, int row, int skippedColumn, int amountOfCitiesInActualSubset);

    int getMinimumColumn(std::vector<std::vector<int>> &cities, int row, int skippedColumn,
                         int amountOfCitiesInActualSubset);

    int reduceRow(std::vector<std::vector<int>> &cites, int row, int amountOfCitiesInActualSubset);

    int reduceColumn(std::vector<std::vector<int>> &cities, int col, int amountOfCitiesInActualSubset);

    int reduceMatrix(std::vector<std::vector<int>> &cities, int amountOfCitiesInActualSubset);

    void removeSubTour(std::vector<std::vector<int>> &activeRoute, int index, std::pair<int, int> &path,
                       int amountOfCitiesInActualSubset);

    void addLastPath(std::vector<std::vector<int>> &m);

    void checkTourCost(std::vector<std::vector<int>> &cities);
};


#endif //SDIZO_3_TRAVELLINGSALESMANPROBLEM_H
