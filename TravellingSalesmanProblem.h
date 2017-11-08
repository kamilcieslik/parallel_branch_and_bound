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
    int lowerBound;                     // path cost till this node
    std::pair<int, int> path;        // path segment (example : 0 1)
    bool bar = false;           // included or excluded path segment
    long parent = -1;     // parent id node (in tree)
};

class TravellingSalesmanProblem {

private:
    int amountOfCities;
    int **arrayOfMatrixOfCities;
    int *optimalWay_GreedyAlgorithmSolution;
    int length;
    bool setGreedyAlgorithm;

    int upperBound;
    std::string fileName;
    std::string graphType;
    bool randomGeneratorData;
    std::deque<Node> treeOfSubsets;
    std::vector<int> optimalWay_BranchAndBoundSolution;

    int StandarizationOfMatrix(std::vector<std::vector<int>> &cities, int amountOfCitiesInActualSubset);

    int GetMinimumRow(std::vector<std::vector<int>> &cities, int row, int skippedColumn,
                      int amountOfCitiesInActualSubset);

    int GetMinimumColumn(std::vector<std::vector<int>> &cities, int row, int skippedColumn,
                         int amountOfCitiesInActualSubset);

    int SubtractMinimalValuesFromTheRows(std::vector<std::vector<int>> &cites, int row,
                                         int amountOfCitiesInActualSubset);

    int SubtractMinimalValuesFromTheColumns(std::vector<std::vector<int>> &cities, int col,
                                            int amountOfCitiesInActualSubset);

    int CalculateCostOfResignation(std::vector<std::vector<int>> &m, std::pair<int, int> &path,
                                   std::pair<int, int> &pos,
                                   int amountOfCitiesInActualSubset);

    void MatrixShortening(std::vector<std::vector<int>> &matrix, int row, int col);

    void EliminationOfSubtour(std::vector<std::vector<int>> &activeRoute, int index, std::pair<int, int> &path,
                              int amountOfCitiesInActualSubset);

    void AddLastPath(std::vector<std::vector<int>> &m);

    void SetOptimalWay(std::vector<std::vector<int>> &m, int index, int begin);

public:
    TravellingSalesmanProblem();

    ~TravellingSalesmanProblem();

    void ReadCitiesFromFile(std::string path);

    void DeleteTravellingSalesman();

    void LoadArrayOfMatrixOfCities(long long int **_cities, int _amountOfCities,
                                   std::string _fileName, std::string _graphType);

    void GenerateRandomCities(int amountOfCities = 0, int maxDistanceBetweenCity = 99);

    void PrintCitiesForTheTravellingSalesman(bool printInputMatrix);

    void GreedyAlgorithm();

    void BranchAndBoundAlgorithm();

    void PrintSolution();

    int GetTourLength(std::string whichAlgorithm);
};


#endif //SDIZO_3_TRAVELLINGSALESMANPROBLEM_H
