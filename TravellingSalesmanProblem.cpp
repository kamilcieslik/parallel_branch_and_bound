//
// Created by mrfarinq on 16.06.17.
//

#include "TravellingSalesmanProblem.h"
#include <stack>

Subset::Subset() : isK1(true), parent(INT_MIN) {
}

TravellingSalesmanProblem::TravellingSalesmanProblem() : amountOfCities(0), arrayOfMatrixOfCities(nullptr),
                                                         optimalWay_GreedyAlgorithmSolution(nullptr),
                                                         upperBound(INT_MAX) {
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

    if (optimalWay_GreedyAlgorithmSolution != nullptr) {
        delete[] optimalWay_GreedyAlgorithmSolution;
        optimalWay_GreedyAlgorithmSolution = nullptr;
    }
}

void TravellingSalesmanProblem::ReadCitiesFromFile(std::string path) {
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
    } else {
        std::cout << "Błąd otwarcia pliku.\n";
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

// -------------------------------------------------------------------
// Algorytm zachłanny dla problemu komiwojażera.
// -------------------------------------------------------------------
void TravellingSalesmanProblem::GreedyAlgorithm() {
    if (arrayOfMatrixOfCities == nullptr)
        throw std::logic_error("Brak miast do przeprowadzenia algorytmu problemu komiwojażera.");

    if (optimalWay_GreedyAlgorithmSolution != nullptr)
        delete[] optimalWay_GreedyAlgorithmSolution;

    setGreedyAlgorithm = true;
    optimalWay_GreedyAlgorithmSolution = new int[amountOfCities];

    bool *visitedCities = new bool[amountOfCities];
    for (int i = 0; i < amountOfCities; i++) {
        visitedCities[i] = false;
    }

    length = 0;
    int currentMinLength;

    int nextCity = 0;
    int city = nextCity;
    visitedCities[city] = true;

    optimalWay_GreedyAlgorithmSolution[0] = nextCity;

    for (auto j = 0; j < amountOfCities - 1; j++) {
        city = nextCity;
        currentMinLength = INT_MAX;
        for (auto i = 0; i < amountOfCities; i++) {
            if (arrayOfMatrixOfCities[city][i] < currentMinLength && !visitedCities[i]) {
                currentMinLength = arrayOfMatrixOfCities[city][i];
                nextCity = i;
            }
        }
        visitedCities[nextCity] = true;
        optimalWay_GreedyAlgorithmSolution[j] = nextCity;
        length += arrayOfMatrixOfCities[city][nextCity];
    }
    optimalWay_GreedyAlgorithmSolution[amountOfCities - 1] = 0;
    length += arrayOfMatrixOfCities[optimalWay_GreedyAlgorithmSolution[amountOfCities - 2]][0];

    delete[] visitedCities;
}

// -------------------------------------------------------------------
// Algorytm podziału i ograniczeń dla problemu komiwojażera.
// -------------------------------------------------------------------
void TravellingSalesmanProblem::BranchAndBoundAlgorithm() {
    if (arrayOfMatrixOfCities == nullptr)
        throw std::logic_error("Brak miast do przeprowadzenia algorytmu problemu komiwojażera.");

    setGreedyAlgorithm = false;

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


void TravellingSalesmanProblem::PrintSolution() {
    std::cout << "\e[1mSolution\e[0m" << std::endl;
    if (setGreedyAlgorithm) {
        std::cout << "\e[1mGreedy Algorithm\e[0m" << std::endl;
    } else {
        std::cout << "\e[1mBranch and Bound Algorithm\e[0m" << std::endl;
    }

    std::cout << "-------------------" << std::endl;

    if (setGreedyAlgorithm) {
        std::cout << "Length\t= " << length << std::endl;
        std::cout << "Path\t= ";
        std::cout << "0 - ";
        for (auto i = 0; i < amountOfCities; i++) {
            if (i == amountOfCities - 1) {
                std::cout << optimalWay_GreedyAlgorithmSolution[i] << std::endl;
            } else {
                std::cout << optimalWay_GreedyAlgorithmSolution[i] << " - ";
            }
        }
    } else {
        std::cout << "Length\t= " << this->upperBound << std::endl;
        std::cout << "Path\t= ";
        for (auto i:this->optimalWay_BranchAndBoundSolution) {
            std::cout << i - 1 << " - ";
        }
        std::cout << "0" << std::endl;
    }
}

void TravellingSalesmanProblem::PrepareMatrix(std::vector<std::vector<int>> &m) {
    for (int i = 0; i < m.size(); i++) {
        m[i][0] = i;
        m[0][i] = i;
    }

    for (auto i = 0; i < m.size() - 1; i++) {
        for (auto j = 0; j < m.size() - 1; j++) {
            m[i + 1][j + 1] = arrayOfMatrixOfCities[i][j];
        }
        m[i + 1][i + 1] = INT_MAX;
    }
}

void TravellingSalesmanProblem::EliminationOfSubtour(std::vector<std::vector<int>> &activeRoute, int index,
                                                     std::pair<int, int> &path) {
    std::vector<std::pair<int, int> > paths;
    // Research of all the included path
    while (index != 0) { // Iterate until we are not arrived at the root
        if (treeOfSubsets[index].isK1) {
            paths.push_back(treeOfSubsets[index].route);
        }
        index = (int) treeOfSubsets[index].parent;
    }

    // Research of the longest subtour
    std::deque<int> subtour = {path.first, path.second};
    bool found = true;
    while (found) {
        found = false;
        for (const std::pair<int, int> &segment : paths) {
            // Check that "segment" go ahead in a subtour
            if (segment.second == subtour.front()) {
                subtour.push_front(segment.second);
                subtour.push_front(segment.first);
                found = true;
                break;
            }
                // Check that "segment" go behind in a subtour
            else if (segment.first == subtour.back()) {
                subtour.push_back(segment.first);
                subtour.push_back(segment.second);
                found = true;
                break;
            }
        }
    }

    std::pair<int, int> pos;
    int founds = 0;
    // Research of the segment to delete in the matrix
    for (int i = 1; i < activeRoute.size(); i++) {
        if (activeRoute[i][0] == subtour.back()) {
            pos.first = i;
            founds++;
        }
        if (activeRoute[0][i] == subtour.front()) {
            pos.second = i;
            founds++;
        }
    }

    // If the segment to delete has been found, then delete it by giving him an infinite cost
    if (founds == 2) {
        //m.setValue(pos.first, pos.second, this->INT_MAX);
        activeRoute[pos.first][pos.second] = INT_MAX;
    }
}

int
TravellingSalesmanProblem::CalculateCostOfResignation(std::vector<std::vector<int>> &m, std::pair<int, int> &path,
                                                      std::pair<int, int> &pos) {
    int max = INT_MIN;
    for (int i = 1; i < m.size(); i++) {
        for (int j = 1; j < m.size(); j++) {
            if (m[i][j] == 0) {
                int val = GetMinimumRow(m, i, j, m.size()) +
                          GetMinimumColumn(m, j, i, m.size());
                if (max < val || max < 0) {
                    max = val;
                    pos.first = i;
                    pos.second = j;

                    path.first = m[i][0];
                    path.second = m[0][j];
                }
            }
        }
    }

    return max;
}

int TravellingSalesmanProblem::GetMinimumRow(std::vector<std::vector<int>> &cities, int row, int skippedColumn,
                                             int amountOfCitiesInActualSubset) {
    int min = INT_MAX;
    for (int i = 1; i < amountOfCitiesInActualSubset; i++) {
        int currentValue = cities[row][i];
        if (currentValue != INT_MAX && i != skippedColumn) {
            min = (min < currentValue ? min : currentValue);
        }
    }
    return min;
}

int TravellingSalesmanProblem::GetMinimumColumn(std::vector<std::vector<int>> &cities, int column, int skippedColumn,
                                                int amountOfCitiesInActualSubset) {
    int min = INT_MAX;
    for (int i = 1; i < amountOfCitiesInActualSubset; i++) {
        int currentValue = cities[i][column];
        if (currentValue != INT_MAX && i != skippedColumn) {
            min = (min < currentValue ? min : currentValue);
        }
    }
    return min;
}

int
TravellingSalesmanProblem::SubtractMinimalValuesFromTheRows(std::vector<std::vector<int>> &cities, int row,
                                                            int amountOfCitiesInActualSubset) {
    int min = GetMinimumRow(cities, row, -1, amountOfCitiesInActualSubset);
    for (int i = 1; i < amountOfCitiesInActualSubset; i++) {
        if (cities[row][i] != INT_MAX) {
            cities[row][i] = cities[row][i] - min;
        }
    }
    return min;
}

int TravellingSalesmanProblem::SubtractMinimalValuesFromTheColumns(std::vector<std::vector<int>> &cities, int col,
                                                                   int amountOfCitiesInActualSubset) {
    int min = GetMinimumColumn(cities, col, -1, amountOfCitiesInActualSubset);
    for (int i = 1; i < amountOfCitiesInActualSubset; i++) {
        if (cities[i][col] != INT_MAX) {
            cities[i][col] = cities[i][col] - min;
        }
    }
    return min;
}

int TravellingSalesmanProblem::StandarizationOfMatrix(std::vector<std::vector<int>> &cities) {
    int minRowTotal = 0;
    for (int i = 1; i < cities.size(); i++) {
        minRowTotal += SubtractMinimalValuesFromTheRows(cities, i, cities.size());
    }

    int minColTotal = 0;
    for (int i = 1; i < cities.size(); i++) {
        minColTotal += SubtractMinimalValuesFromTheColumns(cities, i, cities.size());
    }

    return minRowTotal + minColTotal;
}

void TravellingSalesmanProblem::SetOptimalWay(std::vector<std::vector<int>> &m, int index) {
    Subset K1;

    for (int i = 1; i < 3; i++) {
        for (int j = 1; j < 3; j++) {
            if (m[i][j] == 0) {
                K1.lowerBound = treeOfSubsets.back().lowerBound;
                K1.parent = treeOfSubsets.size() - 1;
                K1.route.first = m[i][0];
                K1.route.second = m[0][j];
                treeOfSubsets.push_back(K1);
            }
        }
    }

    std::vector<std::pair<int, int> > route;
    // Retrieval of the path stored in a branch's tree
    while (index != 0) {    // Iterate until we are not arrived at the root
        if (treeOfSubsets[index].isK1) {     // If it is a node without regret cost
            route.push_back(treeOfSubsets[index].route);   // then we add this segment to the path
        }
        index = (int) treeOfSubsets[index].parent;
    }

    std::vector<int> tour;
    // Research of the path segment containing begin
    int pathSize = (int) route.size();
    for (int i = 0; i < pathSize; i++) {
        if (route[i].first == 1) {
            tour.push_back(route[i].first);
            tour.push_back(route[i].second);
            route.erase(route.begin() + i);
        }
    }

    // Ordering of the rest of the tour
    while (tour.size() != pathSize) {
        for (int i = 0; i < route.size(); i++) {
            if (tour.back() == route[i].first) {
                tour.push_back(route[i].second);
                route.erase(route.begin() + i);
            }
        }
    }

    optimalWay_BranchAndBoundSolution = tour;
}

void TravellingSalesmanProblem::MatrixShortening(std::vector<std::vector<int>> &data, int row, int col) {
    auto it_row = data.begin() + row;
    data.erase(it_row);

    for (int i = 0; i < data.size(); i++) {
        auto it_col = data[i].begin() + col;
        data[i].erase(it_col);
    }
}

int TravellingSalesmanProblem::GetTourLength(std::string whichAlgorithm) {
    if (whichAlgorithm == "greedy")
        return length;
    else if ("branchandbound")
        return upperBound;
    return 0;
}


