//
// Created by mrfarinq on 22.05.17.
//

#include "AlgorithmTest.h"
#include "TimeMeasurement.h"
#include "TravellingSalesmanProblem.h"
#include "TSPLIB_Parser.h"

AlgorithmTest::AlgorithmTest() {

}

AlgorithmTest::~AlgorithmTest() {

}

void AlgorithmTest::TravellingSalesmanProblem_Test_Greedy(int numberOfRepetitions) {
    TravellingSalesmanProblem s;
    TimeMeasurement t;
    std::vector<double> results;
    int amountOfCities = 33;
    double sum = 0;
    std::string path = "test/atsp/ftv33.atsp";
    std::ofstream file;
    file.open("test_atsp_greedy.txt", std::ios::out);
    file << "Test czasów wykonywania algorytmu zachłannego problemu komiwojażera z dnia - " << t.currentDateTime()
         << ".\nDla każdego zestawu danych wyniki uśrednione z " << numberOfRepetitions
         << " losowych instancji." << std::endl << std::endl << std::endl;

    for (auto i = 0; i < 4; i++) {
        results.push_back((double &&) amountOfCities);
        
        //Algorytm zachłanny.
        sum = 0;
        
        for (auto j = 0; j < numberOfRepetitions; j++) {
            TSPLIB_Parser parser(path);
            s.LoadArrayOfMatrixOfCities(parser.GetArrayOfMatrixCities(), parser.GetDimension(),
                                        parser.GetFileName(), parser.GetGraphType());
            t.TimeStart();
            s.GreedyAlgorithm();
            t.TimeStop();
            sum += t.GetTimeInSeconds();
        }
        sum = sum / numberOfRepetitions;
        results.push_back(sum);
        results.push_back(s.GetTourLength("greedy"));
        
        std::cout << "." << std::endl << std::endl;
        
        //Zwiększenie ilości miast - wczytanie danych z kolejnego pliku.
        
        if (path == "test/atsp/ftv33.atsp") {
            path = "test/atsp/ftv35.atsp";
            amountOfCities = 33;
        }
        else if (path == "test/atsp/ftv35.atsp") {
            path = "test/atsp/ftv38.atsp";
            amountOfCities = 38;
        }
        else if (path == "test/atsp/ftv38.atsp") {
            path = "test/atsp/ftv47.atsp";
            amountOfCities = 47;
        }
    }
    
    file << "Il_miast\tZachłanny\tDługość\n";
    for (int i = 0; i < results.size(); i++) {
        file << results[i] << "\t";
        if (((i + 1) % 3) == 0) {
            file << "\n";
        }
    }
    
    file << std::endl << std::endl << "Czas zakończenia testów - " << t.currentDateTime() << "." << std::endl;
    file.close();
    std::cout << "Test zakończony pomyślnie." << std::endl << std::endl;
}

void AlgorithmTest::TravellingSalesmanProblem_Test_BranchAndBound(int numberOfRepetitions) {
    TravellingSalesmanProblem s;
    TimeMeasurement t;
    std::vector<double> results;
    int amountOfCities = 33;
    double sum = 0;
    std::string path = "test/atsp/ftv33.atsp";
    std::ofstream file;
    file.open("test_atsp_branch_and_bound.txt", std::ios::out);
    file << "Test czasów wykonywania algorytmu podziału i ograniczeń problemu komiwojażera z dnia - " << t.currentDateTime()
         << ".\nDla każdego zestawu danych wyniki uśrednione z " << numberOfRepetitions
         << " losowych instancji." << std::endl << std::endl << std::endl;

    for (auto i = 0; i < 4; i++) {
        results.push_back((double &&) amountOfCities);

        //Algorytm podziału i ograniczeń.
        sum = 0;

        for (auto j = 0; j < numberOfRepetitions; j++) {
            TSPLIB_Parser parser(path);
            s.LoadArrayOfMatrixOfCities(parser.GetArrayOfMatrixCities(), parser.GetDimension(),
                                        parser.GetFileName(), parser.GetGraphType());
            t.TimeStart();
            s.BranchAndBoundAlgorithm();
            t.TimeStop();
            sum += t.GetTimeInSeconds();
        }
        sum = sum / numberOfRepetitions;
        results.push_back(sum);
        results.push_back(s.GetTourLength("branchandbound"));

        std::cout << "." << std::endl << std::endl;

        //Zwiększenie ilości miast - wczytanie danych z kolejnego pliku.

        if (path == "test/atsp/ftv33.atsp") {
            path = "test/atsp/ftv35.atsp";
            amountOfCities = 35;
        }
        else if (path == "test/atsp/ftv35.atsp") {
            path = "test/atsp/ftv38.atsp";
            amountOfCities = 38;
        }
        else if (path == "test/atsp/ftv38.atsp") {
            path = "test/atsp/ftv44.atsp";
            amountOfCities = 44;
        }
    }

    file << "Il_miast\tPodziały_i_ograniczenia\tDługość\n";
    for (int i = 0; i < results.size(); i++) {
        file << results[i] << "\t";
        if (((i + 1) % 3) == 0) {
            file << "\n";
        }
    }

    file << std::endl << std::endl << "Czas zakończenia testów - " << t.currentDateTime() << "." << std::endl;
    file.close();
    std::cout << "Test zakończony pomyślnie." << std::endl << std::endl;
}

