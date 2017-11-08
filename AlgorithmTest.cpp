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

void AlgorithmTest::TravellingSalesmanProblem_Test_BruteForce(int numberOfRepetitions) {
    TravellingSalesmanProblem s;
    TimeMeasurement t;
    std::vector<double> results;
    int firstAmountOfCities = 3;
    int amountOfCities = firstAmountOfCities;
    double sum = 0;
    std::ofstream file;
    file.open("test_tsp_brute_force.txt", std::ios::out);
    file << "Test czasów wykonywania algorytmu zupełnego komiwojażera z dnia - " << t.currentDateTime()
         << ".\nDla każdego zestawu danych wyniki uśrednione z " << numberOfRepetitions
         << " losowych instancji." << std::endl << std::endl << std::endl;

    for (auto i = 0; i < 5; i++) {
        results.push_back((double &&) amountOfCities);

        //Algorytm zupełny.
        sum = 0;

        for (auto j = 0; j < numberOfRepetitions; j++) {
            s.GenerateRandomCities(amountOfCities, 99);
            t.TimeStart();
            s.BruteForceAlgorithm();
            t.TimeStop();
            sum += t.GetTimeInSeconds();
        }
        sum = sum / numberOfRepetitions;
        results.push_back(sum);

        std::cout << "." << std::endl << std::endl;

        //Zwiększenie ilości miast

        if (amountOfCities == 3) {
            amountOfCities = 5;
        } else if (amountOfCities == 5) {
            amountOfCities = 8;
        } else if (amountOfCities == 8) {
            amountOfCities = 11;
        } else if (amountOfCities == 11) {
            amountOfCities = 13;
        }
    }

    file << "Il_miast\tZupełny\n";
    for (int i = 0; i < results.size(); i++) {
        file << results[i] << "\t";
        if (((i + 1) % 2) == 0) {
            file << "\n";
        }
    }

    file << std::endl << std::endl << "Czas zakończenia testów - " << t.currentDateTime() << "." << std::endl;
    file.close();
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

