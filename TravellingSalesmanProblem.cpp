//
// Created by mrfarinq on 16.06.17.
//

#include "TravellingSalesmanProblem.h"
#include <stack>

Subset::Subset() : isK1(true), parent(INT_MIN) {
}

TravellingSalesmanProblem::TravellingSalesmanProblem() : amountOfCities(0), arrayOfMatrixOfCities(nullptr),
                                                         optimalWay_BruteForceAlgorithmSolution(nullptr),
                                                         upperBound(INT_MAX), allPossiblePermutations(nullptr) {
}

TravellingSalesmanProblem::~TravellingSalesmanProblem() {
    DeleteTravellingSalesman();
}

void TravellingSalesmanProblem::DeleteTravellingSalesman() {
    for (auto i = 0; i < amountOfCities; i++) {
        delete[] arrayOfMatrixOfCities[i];
    }
    delete[] arrayOfMatrixOfCities;
    arrayOfMatrixOfCities = nullptr;

    if (optimalWay_BruteForceAlgorithmSolution != nullptr) {
        delete[] optimalWay_BruteForceAlgorithmSolution;
        optimalWay_BruteForceAlgorithmSolution = nullptr;
    }

    if (allPossiblePermutations != nullptr) {
        delete[] allPossiblePermutations;
        allPossiblePermutations = nullptr;
    }
}

void TravellingSalesmanProblem::LoadArrayOfMatrixOfCities(long long int **_cities, int _amountOfCities,
                                                          std::string _fileName, std::string _graphType) {
    randomGeneratorData = false;
    if (arrayOfMatrixOfCities != nullptr) {
        for (int i = 0; i < amountOfCities; i++)
            delete[] arrayOfMatrixOfCities[i];
        delete[] arrayOfMatrixOfCities;
    }

    amountOfCities = _amountOfCities;

    arrayOfMatrixOfCities = new int *[amountOfCities];
    for (int i = 0; i < amountOfCities; i++)
        arrayOfMatrixOfCities[i] = new int[amountOfCities];

    for (int i = 0; i < amountOfCities; i++) {
        for (int j = 0; j < amountOfCities; j++) {
            arrayOfMatrixOfCities[i][j] = (int) _cities[i][j];
        }
    }

    fileName = _fileName;
    graphType = _graphType;

    std::cout << "Wczytywanie zakończone pomyślnie." << std::endl;
}

void TravellingSalesmanProblem::ReadCitiesFromNormalFile(std::string path) {
    if (arrayOfMatrixOfCities != nullptr)
        DeleteTravellingSalesman();

    std::fstream file(path, std::ios::in);
    if (file.is_open()) {
        file >> amountOfCities;

        arrayOfMatrixOfCities = new int *[amountOfCities];
        for (auto i = 0; i < amountOfCities; i++) {
            arrayOfMatrixOfCities[i] = new int[amountOfCities];
        }
        int *securityMatrixForReadingPerLine = new int[amountOfCities];

        for (auto i = 0; i < amountOfCities; i++) {
            for (auto j = 0; j < amountOfCities; j++) {
                if (file.fail()) throw std::logic_error("Błąd odczytu danych w pliku.");
                file >> securityMatrixForReadingPerLine[j];
            }

            for (auto j = 0; j < amountOfCities; j++) {
                arrayOfMatrixOfCities[i][j] = securityMatrixForReadingPerLine[j];
            }
        }
        delete[] securityMatrixForReadingPerLine;
        file.close();
        std::cout << "Wczytywanie zakończone pomyślnie." << std::endl;
    } else {
        std::cout << "Błąd otwarcia pliku." << std::endl;
    }
}

void TravellingSalesmanProblem::GenerateRandomCities(int amountOfCities, int maxDistanceBetweenCity) {
    randomGeneratorData = true;
    if (arrayOfMatrixOfCities != nullptr)
        DeleteTravellingSalesman();

    if (amountOfCities == 0) {
        std::cout << "Podaj ilość miast: ";
        std::cin >> this->amountOfCities;
        if (this->amountOfCities < 1) {
            throw std::invalid_argument("Liczba miast nie może być mniejsza od 1.");
        }

        arrayOfMatrixOfCities = new int *[this->amountOfCities];
        for (auto i = 0; i < this->amountOfCities; i++) {
            arrayOfMatrixOfCities[i] = new int[this->amountOfCities];
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dist_distancesBetweenCities(1, maxDistanceBetweenCity);

        for (auto i = 0; i < this->amountOfCities; i++) {
            for (auto j = 0; j < this->amountOfCities; j++) {
                arrayOfMatrixOfCities[i][j] = dist_distancesBetweenCities(gen);
            }
        }
    } else {
        this->amountOfCities = amountOfCities;

        arrayOfMatrixOfCities = new int *[this->amountOfCities];
        for (auto i = 0; i < this->amountOfCities; i++) {
            arrayOfMatrixOfCities[i] = new int[this->amountOfCities];
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dist_distancesBetweenCities(1, maxDistanceBetweenCity);

        for (auto i = 0; i < this->amountOfCities; i++) {
            for (auto j = 0; j < this->amountOfCities; j++) {
                arrayOfMatrixOfCities[i][j] = dist_distancesBetweenCities(gen);
            }
        }
    }
}

void TravellingSalesmanProblem::PrintCitiesForTheTravellingSalesman(bool printInputMatrix) {
    if (arrayOfMatrixOfCities == nullptr)
        throw std::logic_error("Brak miast do wyświetlenia.");

    std::cout << "\e[1mProblem\e[0m" << std::endl;
    std::cout << "-------------------" << std::endl;
    if (printInputMatrix) {
        std::cout << std::endl;
        std::cout << "\t";
        for (auto i = 0; i < amountOfCities; i++) {
            std::cout << i << ".\t";
        }
        std::cout << "\v" << std::endl;
        for (auto i = 0; i < amountOfCities; i++) {
            for (auto j = 0; j < amountOfCities; j++) {
                if (j == 0) {
                    if (arrayOfMatrixOfCities[i][j] < 0) {
                        if (i == j)
                            std::cout << i << ".\t\b" << "\e[1mINF\e[0m";
                        else
                            std::cout << i << ".\t\b" << arrayOfMatrixOfCities[i][j];
                    } else {
                        if (i == j)
                            std::cout << i << ".\t" << "\e[1mINF\e[0m";
                        else
                            std::cout << i << ".\t" << arrayOfMatrixOfCities[i][j];
                    }
                } else {
                    if (arrayOfMatrixOfCities[i][j] < 0) {
                        if (i == j)
                            std::cout << "\t\b" << "\e[1mINF\e[0m";
                        else
                            std::cout << "\t\b" << arrayOfMatrixOfCities[i][j];
                    } else {
                        if (i == j)
                            std::cout << "\t" << "\e[1mINF\e[0m";
                        else
                            std::cout << "\t" << arrayOfMatrixOfCities[i][j];
                    }
                }
            }
            std::cout << "\v" << std::endl;
        }
    }
    if (!randomGeneratorData) {
        std::cout << "Name of TSPLIB file:\t" << fileName << std::endl;
        std::cout << "Graph type:\t\t" << graphType << std::endl;

    } else {
        std::cout << "Dane wejściowe zostały wygenerowane losowo." << std::endl;
        std::cout << "Graph type:\t\tATSP" << std::endl;
    }
    std::cout << "Number of cities:\t" << amountOfCities << std::endl;
}

// --------------------------------------------------------------------------------
// Funkcja rekurencyjna przeszukująca permutacje na potrzeby algorytmu zupełnego.
// --------------------------------------------------------------------------------
void TravellingSalesmanProblem::CalculateTheMostOptimalPermutation(int recursive_param) {
    if (recursive_param == amountOfCities - 1) {
        int lengthInThisPermutation = 0;
        for (auto i = 0; i < amountOfCities - 1; i++) {
            lengthInThisPermutation += arrayOfMatrixOfCities[allPossiblePermutations[i]][allPossiblePermutations[i +
                                                                                                                 1]];
        }
        lengthInThisPermutation += arrayOfMatrixOfCities[allPossiblePermutations[amountOfCities -
                                                                                 1]][allPossiblePermutations[0]];
        if (lengthInThisPermutation < length) {
            length = lengthInThisPermutation;
            for (auto i = 0; i < amountOfCities; i++) {
                optimalWay_BruteForceAlgorithmSolution[i] = allPossiblePermutations[i];
            }
            optimalWay_BruteForceAlgorithmSolution[amountOfCities] = allPossiblePermutations[0];
        }
    } else {
        for (auto i = recursive_param; i < amountOfCities; i++) {
            std::swap(allPossiblePermutations[recursive_param], allPossiblePermutations[i]);
            CalculateTheMostOptimalPermutation(recursive_param + 1);
            std::swap(allPossiblePermutations[recursive_param], allPossiblePermutations[i]);
        }
    }
}

// -------------------------------------------------------------------
// Algorytm zupełny dla problemu komiwojażera.
// -------------------------------------------------------------------
void TravellingSalesmanProblem::BruteForceAlgorithm() {
    if (arrayOfMatrixOfCities == nullptr)
        throw std::logic_error("Brak miast do przeprowadzenia algorytmu problemu komiwojażera.");

    if (optimalWay_BruteForceAlgorithmSolution != nullptr)
        delete[] optimalWay_BruteForceAlgorithmSolution;

    setBruteForceAlgorithm = true;
    optimalWay_BruteForceAlgorithmSolution = new int[amountOfCities + 1];

    allPossiblePermutations = new int[amountOfCities];
    for (int i = 0; i < amountOfCities; i++) {
        allPossiblePermutations[i] = i;
    }

    length = INT_MAX;
    CalculateTheMostOptimalPermutation(0);

    delete[] allPossiblePermutations;
    allPossiblePermutations = nullptr;
}

// -------------------------------------------------------------------
// Algorytm podziału i ograniczeń dla problemu komiwojażera.
// -------------------------------------------------------------------
void TravellingSalesmanProblem::BranchAndBoundAlgorithm() {
    if (arrayOfMatrixOfCities == nullptr)
        throw std::logic_error("Brak miast do przeprowadzenia algorytmu problemu komiwojażera.");

    setBruteForceAlgorithm = false;

    upperBound = INT_MAX;
    treeOfSubsets.clear();
    optimalWay_BranchAndBoundSolution.clear();

    std::vector<std::vector<int>> InputMatrixOfCities(static_cast<unsigned long>(amountOfCities + 1),
                                                      std::vector<int>(static_cast<unsigned long>(amountOfCities + 1)));
    Subset subsetK2;
    subsetK2.isK1 = false;
    Subset subsetK1;

    PrepareMatrix(InputMatrixOfCities);

    std::pair<int, int> positionOfMatrixCell;
    std::pair<int, std::vector<std::vector<int>>> matrix;
    std::stack<std::pair<int, std::vector<std::vector<int>>>> stackOfMatrices;

    int id = 0;
    matrix.first = id;
    matrix.second = InputMatrixOfCities;
    stackOfMatrices.push(matrix);

    while (!stackOfMatrices.empty()) {
        id = stackOfMatrices.top().first;
        std::vector<std::vector<int>> actualMatrix(stackOfMatrices.top().second);
        stackOfMatrices.pop();

        subsetK1.lowerBound = StandarizationOfMatrix(actualMatrix);
        if (id == 0) {
            treeOfSubsets.push_back(subsetK1);
        }

        while (actualMatrix.size() > 3 and treeOfSubsets[id].lowerBound < upperBound) {
            subsetK2.lowerBound = treeOfSubsets[id].lowerBound +
                                  CalculateCostOfResignation(actualMatrix, subsetK1.route, positionOfMatrixCell);
            subsetK2.parent = id;
            treeOfSubsets.push_back(subsetK2);

            if (subsetK2.lowerBound < upperBound) {
                matrix.first = (int) (treeOfSubsets.size() - 1);
                matrix.second = actualMatrix;
                matrix.second[positionOfMatrixCell.first][positionOfMatrixCell.second] = INT_MAX;
                stackOfMatrices.push(matrix);
            }

            MatrixShortening(actualMatrix, positionOfMatrixCell.first, positionOfMatrixCell.second);
            EliminationOfSubtour(actualMatrix, (int) (treeOfSubsets.size() - 1), subsetK1.route);

            subsetK1.lowerBound = treeOfSubsets[id].lowerBound + StandarizationOfMatrix(actualMatrix);
            subsetK1.parent = id;
            treeOfSubsets.push_back(subsetK1);

            id = (int) (treeOfSubsets.size() - 1);
        }

        if (actualMatrix.size() == 3) {
            if (subsetK1.lowerBound < upperBound) {
                upperBound = subsetK1.lowerBound;
                SetOptimalWay(actualMatrix, (int) (treeOfSubsets.size() + 1));
            }
        }
    }
}

void TravellingSalesmanProblem::PrepareMatrix(std::vector<std::vector<int>> &matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        matrix[i][0] = i;
        matrix[0][i] = i;
    }

    for (auto i = 0; i < matrix.size() - 1; i++) {
        for (auto j = 0; j < matrix.size() - 1; j++) {
            matrix[i + 1][j + 1] = arrayOfMatrixOfCities[i][j];
        }
        matrix[i + 1][i + 1] = INT_MAX;
    }
}

void TravellingSalesmanProblem::EliminationOfSubtour(std::vector<std::vector<int>> &activeMatrix, int index,
                                                     std::pair<int, int> &route) {
    std::vector<std::pair<int, int> > _routes;
    // Research of all the included path
    while (index != 0) { // Iterate until we are not arrived at the root
        if (treeOfSubsets[index].isK1) {
            _routes.push_back(treeOfSubsets[index].route);
        }
        index = (int) treeOfSubsets[index].parent;
    }

    // Research of the longest subtour
    std::deque<int> subtour = {route.first, route.second};
    bool found = true;
    while (found) {
        found = false;
        for (const std::pair<int, int> &_route : _routes) {
            // Check that "segment" go ahead in a subtour
            if (_route.second == subtour.front()) {
                subtour.push_front(_route.second);
                subtour.push_front(_route.first);
                found = true;
                break;
            }
                // Check that "segment" go behind in a subtour
            else if (_route.first == subtour.back()) {
                subtour.push_back(_route.first);
                subtour.push_back(_route.second);
                found = true;
                break;
            }
        }
    }

    std::pair<int, int> positionOfMatrixCell;
    int founds = 0;
    // Research of the segment to delete in the matrix
    for (int i = 1; i < activeMatrix.size(); i++) {
        if (activeMatrix[i][0] == subtour.back()) {
            positionOfMatrixCell.first = i;
            founds++;
        }
        if (activeMatrix[0][i] == subtour.front()) {
            positionOfMatrixCell.second = i;
            founds++;
        }
    }

    // If the segment to delete has been found, then delete it by giving him an infinite cost
    if (founds == 2) {
        //m.setValue(pos.first, pos.second, this->INT_MAX);
        activeMatrix[positionOfMatrixCell.first][positionOfMatrixCell.second] = INT_MAX;
    }
}

int
TravellingSalesmanProblem::CalculateCostOfResignation(std::vector<std::vector<int>> &activeMatrix,
                                                      std::pair<int, int> &route,
                                                      std::pair<int, int> &positionOfMatrixCell) {
    int max = INT_MIN;
    for (int i = 1; i < activeMatrix.size(); i++) {
        for (int j = 1; j < activeMatrix.size(); j++) {
            if (activeMatrix[i][j] == 0) {
                int val = GetMinimumRow(activeMatrix, i, j) +
                          GetMinimumColumn(activeMatrix, j, i);
                if (max < val || max < 0) {
                    max = val;
                    positionOfMatrixCell.first = i;
                    positionOfMatrixCell.second = j;

                    route.first = activeMatrix[i][0];
                    route.second = activeMatrix[0][j];
                }
            }
        }
    }

    return max;
}

int TravellingSalesmanProblem::GetMinimumRow(std::vector<std::vector<int>> &activeMatrix, int row, int skippedColumn) {
    int min = INT_MAX;
    for (int i = 1; i < activeMatrix.size(); i++) {
        int currentValue = activeMatrix[row][i];
        if (currentValue != INT_MAX && i != skippedColumn) {
            min = (min < currentValue ? min : currentValue);
        }
    }

    return min;
}

int TravellingSalesmanProblem::GetMinimumColumn(std::vector<std::vector<int>> &activeMatrix, int column,
                                                int skippedColumn) {
    int min = INT_MAX;
    for (int i = 1; i < activeMatrix.size(); i++) {
        int currentValue = activeMatrix[i][column];
        if (currentValue != INT_MAX && i != skippedColumn) {
            min = (min < currentValue ? min : currentValue);
        }
    }

    return min;
}

int
TravellingSalesmanProblem::SubtractMinimalValuesFromTheRows(std::vector<std::vector<int>> &activeMatrix, int row) {
    int min = GetMinimumRow(activeMatrix, row);
    for (int i = 1; i < activeMatrix.size(); i++) {
        if (activeMatrix[row][i] != INT_MAX) {
            activeMatrix[row][i] = activeMatrix[row][i] - min;
        }
    }

    return min;
}

int
TravellingSalesmanProblem::SubtractMinimalValuesFromTheColumns(std::vector<std::vector<int>> &activeMatrix, int col) {
    int min = GetMinimumColumn(activeMatrix, col);
    for (int i = 1; i < activeMatrix.size(); i++) {
        if (activeMatrix[i][col] != INT_MAX) {
            activeMatrix[i][col] = activeMatrix[i][col] - min;
        }
    }

    return min;
}

int TravellingSalesmanProblem::StandarizationOfMatrix(std::vector<std::vector<int>> &activeMatrix) {
    int minRowTotal = 0;
    for (int i = 1; i < activeMatrix.size(); i++) {
        minRowTotal += SubtractMinimalValuesFromTheRows(activeMatrix, i);
    }

    int minColTotal = 0;
    for (int i = 1; i < activeMatrix.size(); i++) {
        minColTotal += SubtractMinimalValuesFromTheColumns(activeMatrix, i);
    }

    return minRowTotal + minColTotal;
}

void TravellingSalesmanProblem::SetOptimalWay(std::vector<std::vector<int>> &activeMatrix, int index) {
    Subset K1;

    for (int i = 1; i < 3; i++) {
        for (int j = 1; j < 3; j++) {
            if (activeMatrix[i][j] == 0) {
                K1.lowerBound = treeOfSubsets.back().lowerBound;
                K1.parent = treeOfSubsets.size() - 1;
                K1.route.first = activeMatrix[i][0];
                K1.route.second = activeMatrix[0][j];
                treeOfSubsets.push_back(K1);
            }
        }
    }

    std::vector<std::pair<int, int> > route;
    while (index != 0) {
        if (treeOfSubsets[index].isK1) {
            route.push_back(treeOfSubsets[index].route);
        }
        index = (int) treeOfSubsets[index].parent;
    }

    std::vector<int> optimalWay;
    int pathSize = (int) route.size();
    for (int i = 0; i < route.size(); i++) {
        if (route[i].first == 1) {
            optimalWay.push_back(route[i].first);
            optimalWay.push_back(route[i].second);
            route.erase(route.begin() + i);
        }
    }

    while (optimalWay.size() != pathSize) {
        for (int i = 0; i < route.size(); i++) {
            if (optimalWay.back() == route[i].first) {
                optimalWay.push_back(route[i].second);
                route.erase(route.begin() + i);
            }
        }
    }

    optimalWay_BranchAndBoundSolution = optimalWay;
}

void TravellingSalesmanProblem::MatrixShortening(std::vector<std::vector<int>> &activeMatrix, int row, int col) {
    auto it_row = activeMatrix.begin() + row;
    activeMatrix.erase(it_row);

    for (int i = 0; i < activeMatrix.size(); i++) {
        auto it_col = activeMatrix[i].begin() + col;
        activeMatrix[i].erase(it_col);
    }
}

void TravellingSalesmanProblem::PrintSolution() {
    std::cout << "\e[1mSolution\e[0m" << std::endl;
    if (setBruteForceAlgorithm) {
        std::cout << "\e[1mFull Search Algorithm\e[0m" << std::endl;
    } else {
        std::cout << "\e[1mBranch and Bound Algorithm\e[0m" << std::endl;
    }

    std::cout << "-------------------" << std::endl;

    if (setBruteForceAlgorithm) {
        std::cout << "Length\t= " << length << std::endl;
        std::cout << "Path\t= ";
        for (auto i = 0; i < amountOfCities; i++) {
            std::cout << optimalWay_BruteForceAlgorithmSolution[i] << " - ";
        }
        std::cout << "0" << std::endl;
    } else {
        std::cout << "Length\t= " << this->upperBound << std::endl;
        std::cout << "Path\t= ";
        for (auto i = 0; i < amountOfCities; i++) {
            std::cout << optimalWay_BranchAndBoundSolution[i] - 1 << " - ";
        }
        std::cout << "0" << std::endl;
    }
}

int TravellingSalesmanProblem::GetTourLength(std::string whichAlgorithm) {
    if (whichAlgorithm == "bruteforce")
        return length;
    else if ("branchandbound")
        return upperBound;
    return 0;
}
