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

struct Subset {
    std::pair<int, int> route;
    int lowerBound;
    long parent;
    bool isK1;

    Subset();
};

class TravellingSalesmanProblem {

private:
    int amountOfCities;
    int **arrayOfMatrixOfCities;
    int *allPossiblePermutations;
    int *optimalWay_BruteForceAlgorithmSolution;
    int length;
    bool setBruteForceAlgorithm;

    int upperBound;
    std::string fileName;
    std::string graphType;
    bool randomGeneratorData;
    std::vector<Subset> treeOfSubsets;
    std::vector<int> optimalWay_BranchAndBoundSolution;

    void PrepareMatrix(std::vector<std::vector<int>> &matrix);

    int StandarizationOfMatrix(std::vector<std::vector<int>> &activeMatrix);

    int GetMinimumRow(std::vector<std::vector<int>> &activeMatrix, int row, int skippedColumn = -1);

    int GetMinimumColumn(std::vector<std::vector<int>> &activeMatrix, int row, int skippedColumn = -1);

    int SubtractMinimalValuesFromTheRows(std::vector<std::vector<int>> &activeMatrix, int row);

    int SubtractMinimalValuesFromTheColumns(std::vector<std::vector<int>> &cities, int col);

    int CalculateCostOfResignation(std::vector<std::vector<int>> &activeMatrix, std::pair<int, int> &route,
                                   std::pair<int, int> &positionOfMatrixCell);

    void MatrixShortening(std::vector<std::vector<int>> &activeMatrix, int row, int col);

    void EliminationOfSubtour(std::vector<std::vector<int>> &activeMatrix, int index, std::pair<int, int> &route);

    void SetOptimalWay(std::vector<std::vector<int>> &activeMatrix, int index);

    void CalculateTheMostOptimalPermutation(int recursive_param);

public:
    TravellingSalesmanProblem();

    ~TravellingSalesmanProblem();

    void DeleteTravellingSalesman();

    void LoadArrayOfMatrixOfCities(long long int **_cities, int _amountOfCities,
                                   std::string _fileName, std::string _graphType);

    void ReadCitiesFromNormalFile(std::string path);

    void GenerateRandomCities(int amountOfCities = 0, int maxDistanceBetweenCity = 99);

    void PrintCitiesForTheTravellingSalesman(bool printInputMatrix);

    void BruteForceAlgorithm();

    void BranchAndBoundAlgorithm();

    void PrintSolution();

    int GetTourLength(std::string whichAlgorithm);
};


#endif //SDIZO_3_TRAVELLINGSALESMANPROBLEM_H
