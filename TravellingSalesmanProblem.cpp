//
// Created by mrfarinq on 16.06.17.
//

#include "TravellingSalesmanProblem.h"
#include <numeric>
#include <stack>

using namespace std;

TravellingSalesmanProblem::TravellingSalesmanProblem()
{
    result=nullptr;
    arrayOfMatrixOfCities=nullptr;
    amountOfCities=0;

}

TravellingSalesmanProblem::~TravellingSalesmanProblem()
{


    if(arrayOfMatrixOfCities!=nullptr){
        for (int i = 0; i < amountOfCities; i++)
            delete[] arrayOfMatrixOfCities[i];
        delete[] arrayOfMatrixOfCities;
    }

    if(result!=nullptr)
    {
        delete []result;
    }

    result=nullptr;
    arrayOfMatrixOfCities=nullptr;

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
            arrayOfMatrixOfCities[i][j] = _cities[i][j];
        }
    }

    fileName = _fileName;
    graphType = _graphType;
}


void TravellingSalesmanProblem::generateFile() {

    if(arrayOfMatrixOfCities!=nullptr){
        for (int i = 0; i < amountOfCities; i++)
            delete[] arrayOfMatrixOfCities[i];
        delete[] arrayOfMatrixOfCities;
    }
    arrayOfMatrixOfCities=nullptr;
    srand(time(NULL));

    cout << "podaj liczbe miast: ";
    cin >> amountOfCities;

    vector <int> help;
    vector <vector < int> > vec;

    for (int i = 0; i < amountOfCities; i++) {
        for (int j = i; j < amountOfCities; j++) {
            vec.push_back(vector<int>(amountOfCities, 0));
        }
    }

    for (int i = 0; i < amountOfCities; i++) {
        for (int j = 0; j < i; j++) {
            vec[i][j] = rand() % 50 + 1;
            vec[j][i] = vec[i][j];
        }
    }


    fstream file("file2.txt", ios::out);
    if (file.good()) {
        file << amountOfCities << "\n";
        for (int i = 0; i < amountOfCities; i++) {
            // file << i << " ";
            for (int j = 0; j < amountOfCities; j++) {
                file << vec[i][j] << " ";
            }
            file << "\n";
        }
        file.close();
    }

    vec.clear();
}

void TravellingSalesmanProblem::showCollection()
{
    cout<<"liczba miast: "<<amountOfCities<<endl;
    cout<<"    ";
    for (int i = 0; i < amountOfCities; i++) {
        cout<<std::left<<setw(10)<< i ;

    }
    cout<<endl;
    cout<<"-----------------------------"<<endl;
    for (int i = 0; i < amountOfCities; i++) {
        cout<<setw(5)<<i<<" |";
        for (int j = 0; j < amountOfCities; j++) {
            cout <<std::left<<setw(10)<< arrayOfMatrixOfCities[i][j];
        }
        cout << endl;
    }

}
void TravellingSalesmanProblem::showResult()
{
    int l = 0;
    std::cout << std::endl;

    for (int k = 0; k <= amountOfCities-1; k++)
    {
        std::cout << result[k] << " - ";
    }
    cout<<result[0]<<endl;
    for (int k = 0; k < amountOfCities - 1; k++)
        l += arrayOfMatrixOfCities[result[k]][result[k + 1]];
    l += arrayOfMatrixOfCities[result[amountOfCities - 1]][result[0]];
    std::cout << "Value = " << l << std::endl;
}


void TravellingSalesmanProblem::readDataFromFile(std::string data)
{


    if(arrayOfMatrixOfCities!=nullptr){
        for (int i = 0; i < amountOfCities; i++)
            delete[] arrayOfMatrixOfCities[i];
        delete[] arrayOfMatrixOfCities;
    }
    arrayOfMatrixOfCities=nullptr;

    fstream file(data, ios::in);

    if (file.is_open()) {
        file >> amountOfCities;

        if (file.fail()){
            cout << "blad wczytywania pliku" << endl;
        }
        else
        {
            arrayOfMatrixOfCities = new int*[amountOfCities];
            for (int i = 0; i < amountOfCities; i++)
                arrayOfMatrixOfCities[i] = new  int[amountOfCities];
            for (int i = 0; i < amountOfCities; i++)
            {
                for (int j = 0; j < amountOfCities; j++)
                    file >> arrayOfMatrixOfCities[i][j];
            }

        }


        file.close();
    }
    else
        cout << "blad wczytywania pliku" << endl;


}


void TravellingSalesmanProblem::GenerateRandomCities(int amountOfCities, int maxDistanceBetweenCity) {
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

void TravellingSalesmanProblem::PrintCitiesForTheTravellingSalesman() {
    if (arrayOfMatrixOfCities == nullptr)
        throw std::logic_error("Brak miast do wyświetlenia.");

    std::cout << "\e[1mProblem\e[0m" << std::endl;
    std::cout << "-------------------" << std::endl;
    std::cout << "Name of TSPLIB file:\t" << fileName << std::endl;
    std::cout << "Graph type:\t\t" << graphType << std::endl;
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

    if (optimalWay_Solution != nullptr)
        delete[] optimalWay_Solution;

    setGreedyAlgorithm = false;

    vector<vector<int>> data(static_cast<unsigned long>(amountOfCities+1), vector<int>(static_cast<unsigned long>(amountOfCities+1)));

    for(int i=0;i<amountOfCities+1;i++){
        data[i][0]= i;
        data[0][i]=i;
    }
    for (auto i = 0; i < amountOfCities; i++) {
        for (auto j = 0; j < amountOfCities; j++) {
            data[i+1][j+1] = arrayOfMatrixOfCities[i][j];
        }

        data[i+1][i+1]=infinity;
    }


    for (auto i = 0; i < data.size(); i++) {
        for (auto j = 0; j < data.size(); j++) {
            std::cout<<"        "<<data[i][j];
        }
        std::cout<<endl;
    }
    Node normalNode;     // node with regret
    Node regretNode;     // node without regret
    regretNode.bar = true;
    pair<int, int> pos;     // var to store the position of a cell in the matrix
    stack<pair<int, vector<vector<int>>> > matrices;    // stack containing the necessary matrix to pursue other branch of the tree
    pair<int, vector<vector<int>>> matrix;   // matrix associated to a node
    int num=data.size();
    // Init of the stack with the initial distances matrix
    matrix.first = 0;
    matrix.second = data;
    matrices.push(matrix);
    vector<vector<int>> m;

    while (!matrices.empty()) {
        // Iterate till the stack is empty
        int id = matrices.top().first;

        vector<vector<int>> m(matrices.top().second);
        matrices.pop();

        // Reduction of the matrix and computation of the minimum sum (raw + col)
        normalNode.cost = reduceMatrix(m,m.size());
        if (id == 0) {      // root tree case
            tree.push_back(normalNode);
        }

        /* Until it ends up with a 2x2 matrix (3x3 du to the indexes storage)
         * and until the current node is lower than the reference value */

        while (m.size() > 3 and tree[id].cost < this->reference) {


            if ((tree.size() - 1) % 10000 == 0) {
                cout << "\r" << tree.size() - 1 << " nodes ..." << std::flush;
            }


            // Compute the node with regret
            regretNode.cost = tree[id].cost + calculateRegret(m, normalNode.path, pos,m.size());


            regretNode.parentNodeKey = id;
            tree.push_back(regretNode);

            // Storing of the matrix
            if (regretNode.cost < this->reference) {
                matrix.first = tree.size() - 1;
                matrix.second = m;
                matrix.second[pos.first][pos.second]=infinity;   // Suppression case i, j pour une potentielle recherche ulterieur
                matrices.push(matrix);
            }

            // Deletion raw col
            removeRow(m,pos.first,m.size());
            removeColumn(m,pos.second,m.size());


            // Subtour deletion
            removeSubTour(m, tree.size() - 1, normalNode.path,m.size());

            // Compute the node without regret
            normalNode.cost = tree[id].cost + reduceMatrix(m,m.size());
            normalNode.parentNodeKey = id;
            tree.push_back(normalNode);

            id = tree.size() - 1;


            //cout<<tree[id].cost<<" "<<this->reference<<endl;
        }
        num=m.size();
        // Update of the best tour and the reference value
        if (m.size() == 3) {
            if (normalNode.cost < this->reference) {
                addLastPath(m);
                this->reference = normalNode.cost;
                this->lastTour = orderPath(tree.size() - 1, 1);

#ifdef DEBUG

                cout << "\r";
                checkTourCost(data);
                cout << "Cost " << this->reference;
                cout << " Tour ";
                for (int i = 0; i < this->lastTour.size(); i++) {
                    cout << this->lastTour[i] << " ";
                }
                cout << "Node " << tree.size() - 1;
                cout << endl;
#endif
                PrintSolution();

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
            std::cout << i << " - ";
        }
        std::cout << "1" << std::endl; //do zmiany (indeks w dół)
    }
}

// Add the two last segments of the tour when the matrix is 2x2
void TravellingSalesmanProblem ::addLastPath(vector<vector<int>> & m) {
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

void TravellingSalesmanProblem::removeSubTour(vector<vector<int>> & activeRoute, int index, pair<int, int> &path,int amountOfCitiesInActualSubset) {
    int size = amountOfCitiesInActualSubset;
    pair<int, int> pos;
    int founds = 0;
    vector<pair<int, int> > paths;

    // Research of all the included path
    while (index != 0) { // Iterate until we are not arrived at the root
        if (tree[index].bar == false) {
            paths.push_back(tree[index].path);
        }
        index = tree[index].parentNodeKey;
    }

    // Research of the longest subtour
    deque<int> subtour = {path.first, path.second};
    bool found = true;
    while (found) {
        found = false;
        for (const pair<int, int> &segment : paths) {
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
        //m.setValue(pos.first, pos.second, this->infinity);
        activeRoute[pos.first][pos.second]=infinity;
    }
}

int TravellingSalesmanProblem::calculateRegret(vector<vector<int>> &m, std::pair<int, int> &path, std::pair<int, int> &pos,int amountOfCitiesInActualSubset ) {
    int size = amountOfCitiesInActualSubset;
    int max = -1;
    for (int i = 1; i < size; i++) {
        for (int j = 1; j < size; j++) {
            if (m[i][j] == 0) {
                int val = getMinimumRow(m, i, j,size) + getMinimumColumn(m, j, i,size);
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

int TravellingSalesmanProblem::getMinimumRow(vector<vector<int>> & cities, int row, int skippedColumn,int amountOfCitiesInActualSubset) {
    int nbCol = amountOfCitiesInActualSubset;
    int min = this->infinity;
    for (int i = 1; i < nbCol; i++) {
        int currentValue = cities[row][ i];
        if (currentValue != this->infinity && i != skippedColumn) {
            min = (min < currentValue ? min : currentValue);
        }
    }
    return min;
}

int TravellingSalesmanProblem::getMinimumColumn(vector<vector<int>> &cities, int column, int skippedColumn,int amountOfCitiesInActualSubset) {
    int nbRow = amountOfCitiesInActualSubset;
    int min = this->infinity;
    for (int i = 1; i < nbRow; i++) {
        int currentValue = cities[i][column];
        if (currentValue != this->infinity && i != skippedColumn) {
            min = (min < currentValue ? min : currentValue);
        }
    }
    return min;
}

int TravellingSalesmanProblem::reduceRow(vector<vector<int>> &cities, int row,int amountOfCitiesInActualSubset) {
    int nbCol = amountOfCitiesInActualSubset;
    int min = getMinimumRow(cities, row,-1,nbCol);
    for (int i = 1; i < nbCol; i++) {
        if (cities[row][i] != this->infinity) {
            cities[row][i] =cities[row][i] - min;
        }
    }
    return min;
}

int TravellingSalesmanProblem::reduceColumn(vector<vector<int>> &cities, int col,int amountOfCitiesInActualSubset) {
    int nbRow = amountOfCitiesInActualSubset;


    int min = getMinimumColumn(cities, col,-1,amountOfCitiesInActualSubset);
    for (int i = 1; i < nbRow; i++) {
        if (cities[i][col] != this->infinity) {
            cities[i][col]= cities[i][col] - min;
        }
    }
    return min;
}

int TravellingSalesmanProblem::reduceMatrix(vector<vector<int>> &cities,int amountOfCitiesInActualSubset){
    int nbRow = amountOfCitiesInActualSubset;
    int minRowTotal = 0;
    for (int i = 1; i < nbRow; i++) {
        minRowTotal += reduceRow(cities, i,nbRow);
    }

    int nbCol = amountOfCitiesInActualSubset;
    int minColTotal = 0;
    for (int i = 1; i < nbCol; i++) {
        minColTotal += reduceColumn(cities, i,nbCol);
    }

    return minRowTotal + minColTotal;
}

std::vector<int> TravellingSalesmanProblem::orderPath(int index, int begin) {
    vector<pair<int, int> > path;
    vector<int> tour;

    // Retrieval of the path stored in a branch's tree
    while (index != 0) {    // Iterate until we are not arrived at the root
        if (tree[index].bar == false) {     // If it is a node without regret cost
            path.push_back(tree[index].path);   // then we add this segment to the path
        }
        index = tree[index].parentNodeKey;
    }

    // Research of the path segment containing begin
    int pathSize = path.size();
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

void TravellingSalesmanProblem::checkTourCost(vector<vector<int>> &cities) {
    int cost = 0;
    int size = this->lastTour.size();
    for (int i = 0; i < size - 1; i++) {
        cost += this->arrayOfMatrixOfCities[this->lastTour[i]][this->lastTour[i + 1]];
    }
    cost += this->arrayOfMatrixOfCities[this->lastTour.back()][this->lastTour.front()];
    cout << "Cost check " << cost << " ";
}
void TravellingSalesmanProblem::removeRow(vector<vector<int>> &data, int row,int amountOfCitiesInActualSubset)
{
    typename vector<vector<int> >::iterator it = data.begin() + row;
    data.erase(it);


}

void TravellingSalesmanProblem::removeColumn(vector<vector<int>> &data, int col,int amountOfCitiesInActualSubset)
{
    for (int row = 0; row < data.size(); row++) {
        typename vector<int>::iterator it = data[row].begin() + col;
        data[row].erase(it);
    }

}


