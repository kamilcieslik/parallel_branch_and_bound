//
// Created by mrfarinq on 26.10.17.
//

#ifndef COURSEPROJECT_SDZ3_TSPLIB_PARSER_H
#define COURSEPROJECT_SDZ3_TSPLIB_PARSER_H

#include <vector>
#include <string>

class TSPLIB_Parser {
private:
    std::string fileName;
    std::string type;
    std::string edgeWeightType;
    std::string edgeWeightFormat;
    std::vector<long long int> numbers;
    int dimension;
    long long int **arrayOfMatrixCities;

public:
    long long int** GetArrayOfMatrixCities();

    int GetDimension();

    std::string GetFileName();

    std::string GetGraphType();

    bool checkParameter(std::string, std::string);

    std::string RemovePotentialSpacesBeforeAndAfter(std::string);

    bool readProblem(std::ifstream &);

    bool GenerateMatrix();

    void EuclidesMatrix();

    void PseudoEuclidesMatrix();

    void FullMatrix();

    void UpperRowMatrix();

    void LowerRowMatrix();

    void UpperDiagRowMatrix();

    void LowerDiagRowMatrix();

    TSPLIB_Parser(std::string path);

    ~TSPLIB_Parser();
};

#endif //COURSEPROJECT_SDZ3_TSPLIB_PARSER_H
