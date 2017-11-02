//
// Created by mrfarinq on 16.06.17.
//

#include "TravellingSalesmanProblem.h"
#include <stack>


TravellingSalesmanProblem::TravellingSalesmanProblem() : amountOfCities(0), arrayOfMatrixOfCities(nullptr),
                                                         optimalWay_Solution(nullptr) {
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

    if (optimalWay_Solution != nullptr) {
        delete[] optimalWay_Solution;
        optimalWay_Solution = nullptr;
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

    if (optimalWay_Solution != nullptr)
        delete[] optimalWay_Solution;

    setGreedyAlgorithm = true;
    optimalWay_Solution = new int[amountOfCities];

    bool *visitedCities = new bool[amountOfCities];
    for (int i = 0; i < amountOfCities; i++) {
        visitedCities[i] = false;
    }

    length = 0;
    int currentMinLength;

    int nextCity = 0;
    int city = nextCity;
    visitedCities[city] = true;

    optimalWay_Solution[0] = nextCity;

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
        optimalWay_Solution[j] = nextCity;
        length += arrayOfMatrixOfCities[city][nextCity];
    }
    optimalWay_Solution[amountOfCities - 1] = 0;
    length += arrayOfMatrixOfCities[optimalWay_Solution[amountOfCities - 2]][0];

    delete[] visitedCities;
}


void TravellingSalesmanProblem::BranchAndBoundAlgorithm() {
    if (arrayOfMatrixOfCities == nullptr)
        throw std::logic_error("Brak miast do przeprowadzenia algorytmu problemu komiwojażera.");

    setGreedyAlgorithm = false;

    std::vector<std::vector<int>> data(static_cast<unsigned long>(amountOfCities + 1),
                                       std::vector<int>(static_cast<unsigned long>(amountOfCities + 1)));

    for (int i = 0; i < amountOfCities + 1; i++) {
        data[i][0] = i;
        data[0][i] = i;
    }

    for (auto i = 0; i < amountOfCities; i++) {
        for (auto j = 0; j < amountOfCities; j++) {
            data[i + 1][j + 1] = arrayOfMatrixOfCities[i][j];
        }
        data[i + 1][i + 1] = INT_MAX;
    }

    Node normalNode;     // node with regret
    Node regretNode;     // node without regret
    regretNode.bar = true;
    std::pair<int, int> pos;     // var to store the position of a cell in the matrix
    std::stack<std::pair<int, std::vector<std::vector<int>>> > matrices;    // stack containing the necessary matrix to pursue other branch of the tree
    std::pair<int, std::vector<std::vector<int>>> matrix;   // matrix associated to a node
    // Init of the stack with the initial distances matrix
    matrix.first = 0;
    matrix.second = data;
    matrices.push(matrix);

    while (!matrices.empty()) {
        // Iterate till the stack is empty
        int id = matrices.top().first;

        std::vector<std::vector<int>> m(matrices.top().second);
        matrices.pop();

        // Reduction of the matrix and computation of the minimum sum (raw + col)
        normalNode.cost = reduceMatrix(m, (int) m.size());
        if (id == 0) {      // root tree case
            tree.push_back(normalNode);
        }

        /* Until it ends up with a 2x2 matrix (3x3 du to the indexes storage)
         * and until the current node is lower than the reference value */

        while (m.size() > 3 and tree[id].cost < this->reference) {
            // Compute the node with regret
            regretNode.cost = tree[id].cost + calculateRegret(m, normalNode.path, pos, (int) m.size());

            regretNode.parentNodeKey = id;
            tree.push_back(regretNode);

            // Storing of the matrix
            if (regretNode.cost < this->reference) {
                matrix.first = (int) (tree.size() - 1);
                matrix.second = m;
                matrix.second[pos.first][pos.second] = INT_MAX;   // Suppression case i, j pour une potentielle recherche ulterieur
                matrices.push(matrix);
            }

            // Deletion raw col
            removeRow(m, pos.first);
            removeColumn(m, pos.second);

            // Subtour deletion
            removeSubTour(m, (int) (tree.size() - 1), normalNode.path, (int) m.size());

            // Compute the node without regret
            normalNode.cost = tree[id].cost + reduceMatrix(m, (int) m.size());
            normalNode.parentNodeKey = id;
            tree.push_back(normalNode);

            id = (int) (tree.size() - 1);

            //cout<<tree[id].cost<<" "<<this->reference<<endl;
        }
        // Update of the best tour and the reference value
        if (m.size() == 3) {
            if (normalNode.cost < this->reference) {
                addLastPath(m);
                this->reference = normalNode.cost;
                this->lastTour = orderPath((int) (tree.size() - 1), 1);

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
                std::cout << optimalWay_Solution[i] << std::endl;
            } else {
                std::cout << optimalWay_Solution[i] << " - ";
            }
        }
    } else {
        std::cout << "Length\t= " << this->reference << std::endl;
        std::cout << "Path\t= ";
        for (auto i:this->lastTour) {
            std::cout << i - 1 << " - "; //do zmiany (indeks w dół)
        }
        std::cout << "0" << std::endl;
    }
}

// Add the two last segments of the tour when the matrix is 2x2
void TravellingSalesmanProblem::addLastPath(std::vector<std::vector<int>> &m) {
    Node normalNode;

    for (int i = 1; i < 3; i++) {
        for (int j = 1; j < 3; j++) {
            if (m[i][j] == 0) {
                normalNode.path.first = m[i][0];
                normalNode.path.second = m[0][j];
                normalNode.cost = tree.back().cost;
                normalNode.parentNodeKey = tree.size() - 1;
                tree.push_back(normalNode);
            }
        }
    }
}

void TravellingSalesmanProblem::removeSubTour(std::vector<std::vector<int>> &activeRoute, int index,
                                              std::pair<int, int> &path,
                                              int amountOfCitiesInActualSubset) {
    int size = amountOfCitiesInActualSubset;
    std::pair<int, int> pos;
    int founds = 0;
    std::vector<std::pair<int, int> > paths;

    // Research of all the included path
    while (index != 0) { // Iterate until we are not arrived at the root
        if (!tree[index].bar) {
            paths.push_back(tree[index].path);
        }
        index = (int) tree[index].parentNodeKey;
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

    // Research of the segment to delete in the matrix
    for (int i = 1; i < size; i++) {
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
TravellingSalesmanProblem::calculateRegret(std::vector<std::vector<int>> &m, std::pair<int, int> &path,
                                           std::pair<int, int> &pos,
                                           int amountOfCitiesInActualSubset) {
    int size = amountOfCitiesInActualSubset;
    int max = -1;
    for (int i = 1; i < size; i++) {
        for (int j = 1; j < size; j++) {
            if (m[i][j] == 0) {
                int val = getMinimumRow(m, i, j, size) + getMinimumColumn(m, j, i, size);
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

int TravellingSalesmanProblem::getMinimumRow(std::vector<std::vector<int>> &cities, int row, int skippedColumn,
                                             int amountOfCitiesInActualSubset) {
    int nbCol = amountOfCitiesInActualSubset;
    int min = INT_MAX;
    for (int i = 1; i < nbCol; i++) {
        int currentValue = cities[row][i];
        if (currentValue != INT_MAX && i != skippedColumn) {
            min = (min < currentValue ? min : currentValue);
        }
    }
    return min;
}

int TravellingSalesmanProblem::getMinimumColumn(std::vector<std::vector<int>> &cities, int column, int skippedColumn,
                                                int amountOfCitiesInActualSubset) {
    int nbRow = amountOfCitiesInActualSubset;
    int min = INT_MAX;
    for (int i = 1; i < nbRow; i++) {
        int currentValue = cities[i][column];
        if (currentValue != INT_MAX && i != skippedColumn) {
            min = (min < currentValue ? min : currentValue);
        }
    }
    return min;
}

int
TravellingSalesmanProblem::reduceRow(std::vector<std::vector<int>> &cities, int row, int amountOfCitiesInActualSubset) {
    int nbCol = amountOfCitiesInActualSubset;
    int min = getMinimumRow(cities, row, -1, nbCol);
    for (int i = 1; i < nbCol; i++) {
        if (cities[row][i] != INT_MAX) {
            cities[row][i] = cities[row][i] - min;
        }
    }
    return min;
}

int TravellingSalesmanProblem::reduceColumn(std::vector<std::vector<int>> &cities, int col,
                                            int amountOfCitiesInActualSubset) {
    int nbRow = amountOfCitiesInActualSubset;

    int min = getMinimumColumn(cities, col, -1, amountOfCitiesInActualSubset);
    for (int i = 1; i < nbRow; i++) {
        if (cities[i][col] != INT_MAX) {
            cities[i][col] = cities[i][col] - min;
        }
    }
    return min;
}

int TravellingSalesmanProblem::reduceMatrix(std::vector<std::vector<int>> &cities, int amountOfCitiesInActualSubset) {
    int nbRow = amountOfCitiesInActualSubset;
    int minRowTotal = 0;
    for (int i = 1; i < nbRow; i++) {
        minRowTotal += reduceRow(cities, i, nbRow);
    }

    int nbCol = amountOfCitiesInActualSubset;
    int minColTotal = 0;
    for (int i = 1; i < nbCol; i++) {
        minColTotal += reduceColumn(cities, i, nbCol);
    }

    return minRowTotal + minColTotal;
}

std::vector<int> TravellingSalesmanProblem::orderPath(int index, int begin) {
    std::vector<std::pair<int, int> > path;
    std::vector<int> tour;

    // Retrieval of the path stored in a branch's tree
    while (index != 0) {    // Iterate until we are not arrived at the root
        if (!tree[index].bar) {     // If it is a node without regret cost
            path.push_back(tree[index].path);   // then we add this segment to the path
        }
        index = (int) tree[index].parentNodeKey;
    }

    // Research of the path segment containing begin
    int pathSize = (int) path.size();
    for (int i = 0; i < pathSize; i++) {
        if (path[i].first == begin) {
            tour.push_back(path[i].first);
            tour.push_back(path[i].second);
            path.erase(path.begin() + i);
        }
    }

    // Ordering of the rest of the tour
    while (tour.size() != pathSize) {
        for (int i = 0; i < path.size(); i++) {
            if (tour.back() == path[i].first) {
                tour.push_back(path[i].second);
                path.erase(path.begin() + i);
            }
        }
    }

    return tour;
}

void TravellingSalesmanProblem::checkTourCost(std::vector<std::vector<int>> &cities) {
    int cost = 0;
    int size = (int) this->lastTour.size();
    for (int i = 0; i < size - 1; i++) {
        cost += this->arrayOfMatrixOfCities[this->lastTour[i]][this->lastTour[i + 1]];
    }
    cost += this->arrayOfMatrixOfCities[this->lastTour.back()][this->lastTour.front()];
    std::cout << "Cost check " << cost << " ";
}

void TravellingSalesmanProblem::removeRow(std::vector<std::vector<int>> &data, int row) {
    typename std::vector<std::vector<int> >::iterator it = data.begin() + row;
    data.erase(it);
}

void TravellingSalesmanProblem::removeColumn(std::vector<std::vector<int>> &data, int col) {
    for (int row = 0; row < data.size(); row++) {
        typename std::vector<int>::iterator it = data[row].begin() + col;
        data[row].erase(it);
    }
}


