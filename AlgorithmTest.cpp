//
// Created by mrfarinq on 22.05.17.
//

#include "AlgorithmTest.h"
#include "TimeMeasurement.h"
#include "TravellingSalesmanProblem.h"
#include <cmath>
#include <vector>
#include <fstream>

AlgorithmTest::AlgorithmTest() {

}

AlgorithmTest::~AlgorithmTest() {

}

void AlgorithmTest::TravellingSalesmanProblem_Test_Greedy(int numberOfRepetitions, int maxDistanceBetweenCity) {
    TravellingSalesmanProblem s;
    TimeMeasurement t;
    std::vector<double> results;
    int firstAmountOfCities = 50;
    int amountOfCities = firstAmountOfCities;
    double sum = 0;
    std::ofstream file;
    file.open("test_tsp_greedy.txt", std::ios::out);
    file << "Test czasów wykonywania algorytmu zachłannego problemu komiwojażera z dnia - " << t.currentDateTime()
         << ".\nDla każdego zestawu danych wyniki uśrednione z " << numberOfRepetitions
         << " losowych instancji." << std::endl << std::endl << std::endl;
    
    
    for (auto i = 0; i < 5; i++) {
        results.push_back((double &&) amountOfCities);
        
        //Algorytm zachłanny.
        sum = 0;
        
        for (auto j = 0; j < numberOfRepetitions; j++) {
            s.GenerateRandomCities(amountOfCities, maxDistanceBetweenCity);
            t.TimeStart();
            s.GreedyAlgorithm();
            t.TimeStop();
            sum += t.GetTimeInSeconds();
        }
        sum = sum / numberOfRepetitions;
        results.push_back(sum);
        
        std::cout << "." << std::endl << std::endl;
        
        //Zwiększenie ilości miast.
        
        if (amountOfCities == 50) {
            amountOfCities = 100;
        } else if (amountOfCities == 100) {
            amountOfCities = 200;
        } else if (amountOfCities == 200) {
            amountOfCities = 350;
        } else if (amountOfCities == 350) {
            amountOfCities = 600;
        }
    }
    
    file << "Il_miast\tZachłanny\n";
    for (int i = 0; i < results.size(); i++) {
        file << results[i] << "\t";
        if (((i + 1) % 2) == 0) {
            file << "\n";
        }
    }
    
    file << std::endl << std::endl << "Czas zakończenia testów - " << t.currentDateTime() << "." << std::endl;
    file.close();
}

void AlgorithmTest::TravellingSalesmanProblem_Test_BranchAndBound(int numberOfRepetitions, int maxDistanceBetweenCity) {
    TravellingSalesmanProblem s;
    TimeMeasurement t;
    std::vector<double> results;
    int firstAmountOfCities = 13;
    int amountOfCities = firstAmountOfCities;
    double sum = 0;
    std::ofstream file;
    file.open("test_tsp_brute_force.txt", std::ios::out);
    file << "Test czasów wykonywania algorytmu podziału i ograniczeń komiwojażera z dnia - " << t.currentDateTime()
         << ".\nDla każdego zestawu danych wyniki uśrednione z " << numberOfRepetitions
         << " losowych instancji." << std::endl << std::endl << std::endl;
    
    
    for (auto i = 0; i < 5; i++) {
        results.push_back((double &&) amountOfCities);
        
        //Algorytm podziału i ograniczeń.
        sum = 0;
        
        for (auto j = 0; j < numberOfRepetitions; j++) {
            s.GenerateRandomCities(amountOfCities, maxDistanceBetweenCity);
            t.TimeStart();
            s.BranchAndBoundAlgorithm();
            t.TimeStop();
            sum += t.GetTimeInSeconds();
        }
        sum = sum / numberOfRepetitions;
        results.push_back(sum);
        
        std::cout << "." << std::endl << std::endl;
        
        //Zwiększenie ilości miast.
        
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
    
    file << "Il_miast\tPodziały_i_ograniczenia\n";
    for (int i = 0; i < results.size(); i++) {
        file << results[i] << "\t";
        if (((i + 1) % 2) == 0) {
            file << "\n";
        }
    }
    
    file << std::endl << std::endl << "Czas zakończenia testów - " << t.currentDateTime() << "." << std::endl;
    file.close();
}

