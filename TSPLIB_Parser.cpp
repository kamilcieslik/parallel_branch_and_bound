//
// Created by mrfarinq on 26.10.17.
//

#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include "TSPLIB_Parser.h"

TSPLIB_Parser::TSPLIB_Parser(std::string path) : arrayOfMatrixCities(nullptr) {
    std::ifstream file(path, std::ios::in);
    if (file.is_open()) {
        if (readProblem(file)) {
        } else {
            throw std::invalid_argument("Błąd odczytu pliku.");
        }
    } else
        throw std::logic_error("Błąd odczytu pliku.");
}

TSPLIB_Parser::~TSPLIB_Parser() {
    if (arrayOfMatrixCities != nullptr) {
        for (auto i = 0; i < dimension; i++)
            delete[] arrayOfMatrixCities[i];
        delete[] arrayOfMatrixCities;
    }
    arrayOfMatrixCities = nullptr;
}

bool TSPLIB_Parser::readProblem(std::ifstream &inputFile) {
    std::string line;
    bool isMatrixType = false;
    bool isCoordinatesType = false;
    while (inputFile) {
        getline(inputFile, line);
        if (line == "EOF" || line == "DISPLAY_DATA_SECTION") {
            break;
        }
        if (line.find(delimiter) != std::string::npos) {
            std::string parameter = line.substr(0, line.find(delimiter));
            std::string value = line.substr(line.find(delimiter) + 1, line.npos);

            if (!checkParameter(trim(parameter), trim(value))) {
                return false;
            }
        }
        if (isMatrixType) {
            std::stringstream stream(line);
            int n;
            while (stream >> n) {
                this->numbers.push_back((long long int &&) n);
            }
        }
        if (isCoordinatesType) {
            std::stringstream stream(line);
            int n;
            stream >> n;
            while (stream >> n) {
                this->numbers.push_back((long long int &&) n);
            }
        }
        if (line == "EDGE_WEIGHT_SECTION") {
            isMatrixType = true;
        }
        if (line == "NODE_COORD_SECTION") {
            isCoordinatesType = true;
        }
    }

    GenerateMatrix();
    return true;
}

bool TSPLIB_Parser::checkParameter(std::string keyword, std::string value) {
    if (keyword == "NAME") {
        this->fileName = value;
    } else if (keyword == "TYPE") {
        if ((value == "TSP") || (value == "ATSP")) {
            this->type = value;
        } else {
            std::cout << "Parametr - " << keyword << " - " << value << " jest niewspierany." << std::endl;
            return 0;
        }
    } else if (keyword == "COMMENT") {}
        //this->comment = value;
    else if (keyword == "DIMENSION")
        this->dimension = stoi(value);
    else if (keyword == "EDGE_WEIGHT_TYPE") {
        if (value == "EXPLICIT")
            this->edgeWeightType = value;
        else if (value == "EUC_2D")
            this->edgeWeightType = value;
        else if (value == "ATT")
            this->edgeWeightFormat = value;
        else {
            std::cout << "Parametr " << keyword << " - " << value << " jest niewspierany." << std::endl;
            return 0;
        }
    } else if (keyword == "EDGE_WEIGHT_FORMAT") {
        if ((value == "FULL_MATRIX") ||
            (value == "UPPER_ROW") ||
            (value == "LOWER_ROW") ||
            (value == "UPPER_DIAG_ROW") ||
            (value == "LOWER_DIAG_ROW") ||
            (value == "UPPER_COL") ||
            (value == "LOWER_COL") ||
            (value == "UPPER_DIAG_COL") ||
            (value == "LOWER_DIAG_COL"))
            this->edgeWeightFormat = value;
        else {
            std::cout << "Parametr " << keyword << " - " << value << " jest niewspierany." << std::endl;
            return 0;
        }
    } else if (keyword == "DISPLAY_DATA_TYPE") {

    } else {
        std::cout << "Nieznany parametr " << keyword << "." << std::endl;
        return false;
    }
    return true;
}

std::string TSPLIB_Parser::trim(std::string s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

bool TSPLIB_Parser::GenerateMatrix() {
    if (arrayOfMatrixCities != nullptr) {
        for (auto i = 0; i < dimension; i++)
            delete[] arrayOfMatrixCities[i];
        delete[] arrayOfMatrixCities;
    }
    arrayOfMatrixCities = new long long int *[dimension];
    for (auto i = 0; i < dimension; i++)
        arrayOfMatrixCities[i] = new long long int[dimension];

    if (edgeWeightType == "EUC_2D")
        EuclidesMatrix();

    if (edgeWeightFormat == "ATT")
        PseudoEuclidesMatrix();

    if (edgeWeightFormat == "FULL_MATRIX")
        FullMatrix();

    if ((edgeWeightFormat == "UPPER_ROW") || (edgeWeightFormat == "LOWER_COL"))
        UpperRowMatrix();

    if ((edgeWeightFormat == "LOWER_ROW") || (edgeWeightFormat == "UPPER_COL"))
        LowerRowMatrix();

    if ((edgeWeightFormat == "UPPER_DIAG_ROW") || (edgeWeightFormat == "LOWER_DIAG_COL"))
        UpperDiagRowMatrix();

    if ((edgeWeightFormat == "LOWER_DIAG_ROW") || (edgeWeightFormat == "UPPER_DIAG_COL"))
        LowerDiagRowMatrix();

    return true;
}

void TSPLIB_Parser::FullMatrix() {
    for (auto i = 0; i < this->dimension; i++) {
        for (auto j = 0; j < this->dimension; j++) {
            if (i != j) {
                arrayOfMatrixCities[i][j] = this->numbers[i * this->dimension + j];
            }
        }
    }
}

void TSPLIB_Parser::UpperRowMatrix() {
    int counter = 0;
    for (auto i = 0; i < this->dimension - 1; i++) {
        for (auto j = i + 1; j < this->dimension; j++) {
            arrayOfMatrixCities[i][j] = this->numbers[counter];
            arrayOfMatrixCities[j][i] = this->numbers[counter];
            counter++;
        }
    }
}

void TSPLIB_Parser::LowerRowMatrix() {
    int counter = 0;
    for (auto i = 1; i < this->dimension; i++) {
        for (auto j = 0; j < i; j++) {
            arrayOfMatrixCities[i][j] = this->numbers[counter];
            arrayOfMatrixCities[j][i] = this->numbers[counter];
            counter++;
        }
    }
}

void TSPLIB_Parser::UpperDiagRowMatrix() {
    int counter = 0;
    for (auto i = 0; i < this->dimension; i++) {
        for (auto j = i; j < this->dimension; j++) {
            if (i != j) {
                arrayOfMatrixCities[i][j] = this->numbers[counter];
                arrayOfMatrixCities[j][i] = this->numbers[counter];
            }
            counter++;
        }
    }
}

void TSPLIB_Parser::LowerDiagRowMatrix() {
    int counter = 0;
    for (auto i = 0; i < this->dimension; i++) {
        for (auto j = 0; j < i + 1; j++) {
            if (j != i) {
                arrayOfMatrixCities[i][j] = this->numbers[counter];
                arrayOfMatrixCities[j][i] = this->numbers[counter];
            }
            counter++;
        }
    }
}

void TSPLIB_Parser::EuclidesMatrix() {
    int i = 0;
    int j = 0;
    for (auto k = 0; k < (2 * dimension); k = k + 2) {
        for (auto l = 0; l < (2 * dimension); l = l + 2) {

            arrayOfMatrixCities[i][j] = (long long int) round(
                    sqrt((numbers[k] - numbers[l]) * (numbers[k] - numbers[l]) +
                         (numbers[k + 1] - numbers[l + 1]) * (numbers[k + 1] - numbers[l + 1])));
            j++;
        }
        j = 0;
        i++;
    }
}

void TSPLIB_Parser::PseudoEuclidesMatrix() {
    double rkl;
    long long int tkl;
    int i = 0;
    int j = 0;
    for (auto k = 0; k < (2 * dimension); k = k + 2) {
        for (auto l = 0; l < (2 * dimension); l = l + 2) {

            rkl = sqrt(((numbers[k] - numbers[l]) * (numbers[k] - numbers[l]) +
                        (numbers[k + 1] - numbers[l + 1]) * (numbers[k + 1] - numbers[l + 1])) / 10.0);
            tkl = (long long int) round(rkl);
            if (tkl < rkl) arrayOfMatrixCities[i][j] = tkl + 1;
            else arrayOfMatrixCities[i][j] = tkl;

            j++;
        }
        j = 0;
        i++;
    }
}


int TSPLIB_Parser::GetDimension() {
    return dimension;
}

std::string TSPLIB_Parser::GetFileName() {
    return fileName;
}

std::string TSPLIB_Parser::GetGraphType() {
    return type;
}

long long int **TSPLIB_Parser::GetArrayOfMatrixCities() {
    return arrayOfMatrixCities;
}
